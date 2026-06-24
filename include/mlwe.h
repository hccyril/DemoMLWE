//==================================================================================
// mlwe.h
//
// 基于 OpenFHE 的 Module-LWE（MLWE）模块实现头文件
//
//   模块化容错学习问题（Module Learning With Errors, Module-LWE，简称 MLWE）
// 是介于标准 LWE 与 Ring-LWE 之间的一类格密码学难题。它由 Brakerski、Gentry、
// Vaikuntanathan（BGV）以及 Langlois-Stehlé 等人推广，目前已成为
//   Kyber（ML-KEM）、Dilithium（ML-DSA）等后量子标准化算法的核心困难性假设。
//
//   本模块直接复用 OpenFHE 的底层多项式环原语（NativePoly / ILParams2N /
// DiscreteGaussianGenerator）来“从零”实现一个标准的 MLWE：
//     1) 选取环 R_q = Z_q[x] / (x^n + 1)，其中 n 为 2 的幂，q 为模数；
//     2) 模块秩 k：把密钥/密文视为 R_q 上的 k 维向量；
//     3) 采样矩阵 A ∈ R_q^{k×n_R}（公钥分量）、小噪声 s、e；
//     4) 公钥 b = A·s + e（在 R_q 上计算），私钥即 s。
//
//   MLWE 与 RLWE 的关键差别：
//     - RLWE：A 是单个环元素（k = 1 的特例）。
//     - MLWE：A 是 R_q 上的 k×k（或 k×m）矩阵，秘密与噪声均为 R_q 上的向量。
//   MLWE 在同等安全强度下可在密钥尺寸 / 运算速度 / 安全规约之间取得更平滑的
// 折中，这也是 Kyber 等算法选择它的原因。
//
// 数学约定（贯穿本文件）：
//   - 环 R         = Z[x] / (x^n + 1)
//   - 商环 R_q     = R / qR = Z_q[x] / (x^n + 1)，系数模 q
//   - 多项式乘法   均在 R_q 中进行，即先多项式乘法再对 (x^n + 1, q) 取模
//   - 小噪声分布   χ 为以 0 为中心、标准差 σ 的离散高斯分布（系数级采样）
//
// 文件结构：
//   - include/mlwe.h        : 本文件，类与函数声明 + 公式注释
//   - src/mlwe.cpp          : 实现
//   - src/mlwe_demo.cpp     : 演示程序（密钥生成 / 加密 / 解密 / 正确性验证）
//   - CMakeLists.txt        : 构建脚本
//   - README.md             : 说明文档
//==================================================================================

#ifndef OPENFHE_MLWE_H
#define OPENFHE_MLWE_H

#include <random>
#include <string>
#include <vector>
// OpenFHE 头文件
#include "openfhecore.h"          // 核心数据类型（NativeInteger 等）
#include "lattice/lat-defaults.h" // NativePoly、ILParams2N 等多项式环默认实现
#include "math/nbtheory.h"        // 数论工具：FirstPrime、RootOfUnity 等

namespace mlwe {

// 元素类型统一使用 OpenFHE 的 NativePoly：
//   NativePoly 是单精度（64-bit 或机器字长）多项式环元素实现，
//   内部以系数序列表示一个 R_q 中的元素，乘法自动对 (x^n+1, q) 取模。
// 使用 NativePoly 而非 DCRTPoly：MLWE 教学实现只需单个素模数 q，
// 不需要 RNS（双 CRT）链，NativePoly 更轻量、更直观。
using RingElement = lbcrypto::NativePoly;

//------------------------------------------------------------------------------
// 结构体：MLWE 参数集
//------------------------------------------------------------------------------
// 对应 MLWE 问题的可调参数。以 Kyber-512 风格为参考（n=256, q=3329, k=2），
// 但为便于演示，默认采用更小的安全参数（运行快、可读性好）。
struct MLWEParams {
    usint n;          // 多项式环的维数（次数），x^n+1 中 n，必须为 2 的幂
    usint k;          // 模块的秩（rank），即 R_q 向量的维数
    usint q;          // 模数（必须为素数，且满足 1 mod 2n，使 x^n+1 在 Z_q 上完全分裂）
    usint nu;         // 噪声的离散高斯参数 λ（标准差 σ 的倒数 / 截断半径相关）
    usint base;       // （可选）用于 gadget 分解的基 b，构造 G 矩阵时使用

    MLWEParams(usint n_ = 256, usint k_ = 2, usint q_ = 7681,
               usint nu_ = 8, usint base_ = 0)
        : n(n_), k(k_), q(q_), nu(nu_), base(base_) {}
};

//------------------------------------------------------------------------------
// 类：MLWE 公钥 / 私钥 / 密文
//------------------------------------------------------------------------------
// MLWE 的公钥为 (A, b = A·s + e)；私钥为 s。
// 这里把“公钥”打包为一个结构体方便传递。
struct MLWEPublicKey {
    std::vector<std::vector<RingElement>> A;  // A ∈ R_q^{k×k}：k×k 环元素矩阵
    std::vector<RingElement> b;               // b ∈ R_q^k：k 维环元素向量
};

// 私钥为 k 维环元素向量 s ∈ R_q^k
using MLWESecretKey = std::vector<RingElement>;

// 密文也用一对 (c0, c1) 表示：这是对消息 m ∈ R_q（或其量化版本）的加密结果。
//   c0 = b^T · r + m        （或 Δ·m，取决于是否做缩放）
//   c1 = A^T · r
// 其中 r 为加密时采样的另一组小噪声。该密文格式与 LWE/RLWE 公钥加密一致，
// 便于验证“解密 = c0 - s^T · c1 ≈ m”的正确性。
struct MLWECiphertext {
    RingElement c0;  // 标量环元素（与消息同类型）
    std::vector<RingElement> c1;  // k 维环元素向量
};

//------------------------------------------------------------------------------
// 类：MLWE 上下文（封装参数、环参数、分布采样器）
//------------------------------------------------------------------------------
// 该类持有：
//   1) MLWE 参数集（n, k, q, ...）；
//   2) OpenFHE 的环参数对象 ILParams2N（描述 R_q）；
//   3) 离散高斯采样器（用于秘密与噪声）；
//   4) 均匀采样（用于公钥矩阵 A）。
//
// 把这些状态封装在一个上下文对象里，避免每次采样都要重新构造环参数。
class MLWEContext {
public:
    // 构造函数：根据参数集初始化环参数与分布
    //   - 会校验 q 是否为素数且满足 q ≡ 1 (mod 2n)（这是让 x^n+1 在 Z_q 上
    //     存在 2n 次单位根、从而 OpenFHE 能找到合适生成元的条件）；
    //   - 若不满足则抛出 std::invalid_argument。
    explicit MLWEContext(const MLWEParams& params);

    // 禁止拷贝（内部持有非可拷贝的 OpenFHE 对象）
    MLWEContext(const MLWEContext&) = delete;
    MLWEContext& operator=(const MLWEContext&) = delete;

    //--------------------------------------------------------------
    // 基本访问
    //--------------------------------------------------------------
    const MLWEParams& GetParams() const { return m_params; }

    // 获取环参数（用于构造新的 RingElement）
    const lbcrypto::ILParams2N& GetElementParams() const {
        return *m_elemParams;
    }

    //--------------------------------------------------------------
    // 环元素工厂方法
    //--------------------------------------------------------------
    // 返回一个零元素：系数全为 0 的 R_q 元素
    RingElement MakeElement() const;

    // 返回一个常数元素：所有系数等于常数 c 的 R_q 元素
    //   c ∈ Z_q，对应环中“整数 c”的提升。
    RingElement MakeConstantElement(int64_t c) const;

    //--------------------------------------------------------------
    // 采样器
    //--------------------------------------------------------------
    // 从以 0 为中心、标准差 σ 的离散高斯分布采样一个 R_q 元素
    //   - 系数级独立采样，每个系数取自 χ_σ；
    //   - 用于秘密 s 与噪声 e。
    //   公式：a_i ← χ(i)。
    //   实现细节：本实现使用 std::normal_distribution 做连续高斯采样后四舍五入
    //   到最近整数（“舍入高斯”，rounded Gaussian），这是合法的 χ 分布实例，
    //   且不依赖 OpenFHE 内部分布 API，从而保证跨版本可移植。
    RingElement SampleGaussian() const;

    // 从均匀分布采样一个 R_q 元素（每个系数均匀取自 Z_q）
    //   - 用于公钥矩阵 A 的元素。
    RingElement SampleUniform() const;

private:
    MLWEParams m_params;
    // 环参数：描述 R_q（模数 q、维数 n、2n 次单位根）
    // 注意：使用 shared_ptr 持有，因为 OpenFHE 的 RingElement 内部以
    // shared_ptr<const ILParams2N> 引用环参数。
    std::shared_ptr<lbcrypto::ILParams2N> m_elemParams;

    // C++11 随机数引擎（用于高斯与均匀系数采样）
    // 使用 mutable 因为采样方法是 const 逻辑操作，但引擎内部状态会改变。
    mutable std::mt19937_64 m_rng;
};

//------------------------------------------------------------------------------
// 类：MLWE 方案（密钥生成 / 加密 / 解密 / 算术）
//------------------------------------------------------------------------------
class MLWEScheme {
public:
    explicit MLWEScheme(std::shared_ptr<MLWEContext> ctx)
        : m_ctx(std::move(ctx)) {}

    //--------------------------------------------------------------
    // 密钥生成：KeyGen
    //--------------------------------------------------------------
    // 算法（标准 MLWE 公钥生成）：
    //   输入：参数 (n, k, q, σ)
    //   1. 均匀采样 A ← R_q^{k×k}；                       // 公共矩阵
    //   2. 从小噪声分布采样 s ← χ^k，e ← χ^k；             // 私钥、噪声
    //   3. 计算 b = A·s + e ∈ R_q^k；                      // 在 R_q 中
    //   4. 公钥 pk = (A, b)，私钥 sk = s。
    //
    // 输出：
    //   pk：见 MLWEPublicKey
    //   sk：见 MLWESecretKey
    void KeyGen(MLWEPublicKey& pk, MLWESecretKey& sk) const;

    //--------------------------------------------------------------
    // 加密：Encrypt
    //--------------------------------------------------------------
    // 算法（MLWE-PKE 的标准加密，公钥重随机化风格）：
    //   输入：公钥 (A, b)、消息 m ∈ R_q
    //   1. 采样小噪声 r ← χ^k；
    //   2. c1 = A^T · r        ∈ R_q^k；       // “掩码”部分
    //   3. c0 = b^T · r + m    ∈ R_q；         // 消息部分（含噪声）
    //   4. 输出密文 (c0, c1)。
    //
    // 注：这里为简洁直接把 m 当作 R_q 元素加密（不做缩放）。
    //     解密时得到的将是 m + 小噪声，可作为正确性演示。
    MLWECiphertext Encrypt(const MLWEPublicKey& pk,
                           const RingElement& m) const;

    //--------------------------------------------------------------
    // 解密：Decrypt
    //--------------------------------------------------------------
    // 算法：
    //   输入：密文 (c0, c1)、私钥 s
    //   计算：m' = c0 - s^T · c1 ∈ R_q
    //
    // 正确性分析（解密误差为何“小”）：
    //   m' = c0 - s^T · c1
    //      = (b^T · r + m) - s^T · (A^T · r)
    //      = (A·s + e)^T · r + m - s^T · A^T · r
    //      = s^T · A^T · r + e^T · r + m - s^T · A^T · r   // 展开 b = A·s + e
    //      = m + e^T · r                                    // 线性项相互抵消
    //   其中 e^T · r 是两个小噪声多项式的内积，其系数仍是“小”的（量级约
    //   n · σ^2），远小于 q/2，因此 m' 与 m 在 R_q 中逐系数接近。
    RingElement Decrypt(const MLWECiphertext& ct,
                        const MLWESecretKey& sk) const;

    //--------------------------------------------------------------
    // 工具：环元素算术封装
    //--------------------------------------------------------------
    // 计算向量内积 a^T · b（两个 k 维 R_q 向量）
    RingElement InnerProduct(const std::vector<RingElement>& a,
                             const std::vector<RingElement>& b) const;

    // 计算 c0 - s^T · c1（解密核心运算）
    RingElement DecryptCore(const RingElement& c0,
                            const std::vector<RingElement>& c1,
                            const MLWESecretKey& sk) const;

private:
    std::shared_ptr<MLWEContext> m_ctx;
};

//------------------------------------------------------------------------------
// 工具函数（自由函数，便于演示与调试）
//------------------------------------------------------------------------------
// 把一个 RingElement 的所有系数打印出来（用于检查噪声大小）
std::string ElementToString(const RingElement& e, bool all = false);

// 把一个 RingElement 的系数逐项转为 long long 向量返回
std::vector<int64_t> ElementToVector(const RingElement& e);

// 把 long long 向量打包成一个 RingElement
RingElement VectorToElement(const std::vector<int64_t>& coeffs,
                            const MLWEContext& ctx);

// 计算两个 RingElement 逐系数的最大绝对差（衡量解密误差的常用指标）
int64_t MaxCoeffAbsDiff(const RingElement& a, const RingElement& b);

}  // namespace mlwe

#endif  // OPENFHE_MLWE_H
