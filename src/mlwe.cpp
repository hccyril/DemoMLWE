//==================================================================================
// mlwe.cpp
//
// 基于 OpenFHE 的 Module-LWE（MLWE）实现。
// 本文件实现了 mlwe.h 中声明的全部类与函数。每个函数前的注释都会再次说明
// 它所对应的数学公式与步骤，方便对照阅读。
//
// 依赖的 OpenFHE 原语（均为稳定公开 API）：
//   - lbcrypto::NativePoly            : 多项式环元素类型（R_q 中的元素）
//   - lbcrypto::ILParams2N            : 环参数（模数 q、维数 n、2n 次本原单位根）
//   - lbcrypto::FirstPrime / RootOfUnity: 数论工具
//
// OpenFHE 多项式环元素的运算（均自动在 R_q = Z_q[x]/(x^n+1) 中完成）：
//   - operator+ : 系数逐项相加，模 q
//   - operator* : 多项式乘法（NTT 加速），并对 x^n+1 与 q 取模（负循环卷积）
//   - operator- : 系数逐项相减，模 q
// 这些运算已封装了“模 x^n+1 的负循环卷积”，等价于密码学意义上的环乘法。
//==================================================================================

#include "mlwe.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <sstream>

namespace mlwe {

//------------------------------------------------------------------------------
// MLWEContext 构造：根据参数构造环参数与分布采样器
//------------------------------------------------------------------------------
//
// 数学背景：
//   我们要在 R_q = Z_q[x] / (x^n + 1) 上工作。
//   ILParams2N 封装了三件事：(1) 分圆域阶数 m = 2n； (2) 模数 q；
//   (3) 一个 2n 次本原单位根 M ∈ Z_q（用于 NTT）。
//   构造合法环参数要求 q 为素数且 q ≡ 1 (mod 2n)。
//
MLWEContext::MLWEContext(const MLWEParams& params)
    : m_params(params) {
    // ------------------------------------------------------------
    // 步骤 1：校验环维数 n
    // ------------------------------------------------------------
    // n 必须为 2 的幂（分圆多项式 Φ_{2n}(x)=x^n+1 的标准要求，
    // 也是 OpenFHE 对环维数的强制约束）。
    if (params.n == 0 || (params.n & (params.n - 1)) != 0) {
        throw std::invalid_argument(
            "MLWEContext: n 必须是 2 的幂（cyclotomic 维数）。");
    }

    // ------------------------------------------------------------
    // 步骤 2：处理模数 q
    // ------------------------------------------------------------
    // OpenFHE 的 ILParams2N 需要模数 q 满足：
    //   (1) q 为素数；
    //   (2) q ≡ 1 (mod 2n)，保证 2n 次本原单位根 ω_{2n} 存在于 Z_q 中。
    // 这两点共同保证 x^n+1 在 Z_q 上完全分裂，NTT 变换可正常工作。
    //
    // 这里只做“合法性校验”：若用户给定的 q 不满足条件，直接抛出异常，
    // 提示用户更换参数。这样把参数选择的责任交给上层（避免静默改参数
    // 导致与用户预期不一致）。
    usint n = params.n;
    usint q = params.q;
    if (q % (2 * n) != 1) {
        std::ostringstream oss;
        oss << "MLWEContext: 模数 q=" << q << " 必须满足 q ≡ 1 (mod 2n="
            << (2 * n) << ")，"
            << "以保证 x^n+1 在 Z_q 上完全分裂。请更换 q。";
        throw std::invalid_argument(oss.str());
    }

    // ------------------------------------------------------------
    // 步骤 3：构造 ILParams2N（环参数对象）
    // ------------------------------------------------------------
    // ILParams2N(cyclotomicOrder, modulus, rootOfUnity)
    //   - cyclotomicOrder = 2n（因为 Φ_{2n}=x^n+1，分圆域阶数为 2n）
    //   - modulus         = q
    //   - rootOfUnity     = q 中的 2n 次本原单位根，由 RootOfUnity(2n, q) 计算
    // 注意：modulus / root 必须是 RingElement::Integer（NativeInteger）类型，
    //   OpenFHE 的 RootOfUnity 是模板函数，参数与返回类型一致。
    usint cyclotomicOrder = 2 * n;
    typename RingElement::Integer modulus(q);
    typename RingElement::Integer root =
        lbcrypto::RootOfUnity(cyclotomicOrder, modulus);

    // make_shared 创建 ILParams2N 对象。NativePoly 在构造时会以
    // shared_ptr<const ILParams2N> 的形式引用它，因此对象生命周期由
    // shared_ptr 管理，安全。
    m_elemParams = std::make_shared<lbcrypto::ILParams2N>(
        cyclotomicOrder, modulus, root);

    // ------------------------------------------------------------
    // 步骤 4：初始化随机数引擎
    // ------------------------------------------------------------
    // 使用标准库 mt19937_64 作为底层 PRNG。实际部署中可换成更强的熵源
    //（如 std::random_device 播种）。这里以固定种子+时间种子混合，
    // 既保证可复现性又避免每次结果完全一致。
    std::random_device rd;
    m_rng.seed(static_cast<std::mt19937_64::result_type>(rd()));
}

//------------------------------------------------------------------------------
// MakeElement：构造 R_q 的零元素
//------------------------------------------------------------------------------
RingElement MLWEContext::MakeElement() const {
    // 用环参数构造一个零元素，并显式以系数表示存储。
    // Format::COEFFICIENT 表示系数向量形式（与 EVALUATION/NTT 域相对）。
    return RingElement(m_elemParams, lbcrypto::COEFFICIENT);
}

//------------------------------------------------------------------------------
// MakeConstantElement：构造常数元素 c·1 ∈ R_q
//------------------------------------------------------------------------------
RingElement MLWEContext::MakeConstantElement(int64_t c) const {
    // 步骤：
    //   1. 构造零元素（系数表示）；
    //   2. 把第 0 个系数（即常数项，对应 1∈R_q 的系数）置为 c；
    //   3. OpenFHE 内部会自动对 c 模 q 归一化（含负数）。
    RingElement e(m_elemParams, lbcrypto::COEFFICIENT);
    e.SetValueAtIndex(0, c);
    return e;
}

//------------------------------------------------------------------------------
// SampleGaussian：从离散高斯分布采样一个 R_q 元素
//------------------------------------------------------------------------------
// 数学：对每个系数 i = 0,...,n-1 独立地从中心化离散高斯
//         χ_σ(x) ∝ exp(-x^2 / (2σ^2))   (x ∈ Z)
//       采样，得到 a_i，组成环元素 a(x) = Σ a_i x^i。
//
// 实现：使用“舍入高斯”（rounded Gaussian）——先采样连续高斯 N(0, σ²)，
//   再四舍五入到最近整数。它是合法的 χ 分布实例，且不依赖 OpenFHE 内部分布
//   API，跨版本可移植。σ 取自 m_params.nu。
RingElement MLWEContext::SampleGaussian() const {
    std::normal_distribution<double> dist(0.0, static_cast<double>(m_params.nu));
    RingElement e(m_elemParams, lbcrypto::COEFFICIENT);
    for (usint i = 0; i < m_params.n; ++i) {
        double r = dist(m_rng);
        int64_t v = static_cast<int64_t>(std::llround(r));
        e.SetValueAtIndex(i, v);
    }
    return e;
}

//------------------------------------------------------------------------------
// SampleUniform：从均匀分布采样一个 R_q 元素
//------------------------------------------------------------------------------
// 数学：每个系数独立均匀取自 Z_q = {0,1,...,q-1}。
// 用标准库 uniform_int_distribution 在 [0, q-1] 内采样。
RingElement MLWEContext::SampleUniform() const {
    std::uniform_int_distribution<int64_t> dist(
        0, static_cast<int64_t>(m_params.q - 1));
    RingElement e(m_elemParams, lbcrypto::COEFFICIENT);
    for (usint i = 0; i < m_params.n; ++i) {
        e.SetValueAtIndex(i, dist(m_rng));
    }
    return e;
}

//==============================================================================
// MLWEScheme 实现
//==============================================================================

//------------------------------------------------------------------------------
// KeyGen：密钥生成
//------------------------------------------------------------------------------
// 算法（标准 MLWE 公钥生成）：
//   输入：参数 (n, k, q, σ)
//   1. 均匀采样 A ← R_q^{k×k}；                       // 公共矩阵
//   2. 从小噪声分布采样 s ← χ^k，e ← χ^k；             // 私钥、噪声
//   3. 计算 b = A·s + e ∈ R_q^k；                      // 在 R_q 中
//   4. 公钥 pk = (A, b)，私钥 sk = s。
//
// 说明：本实现中 A 为 k×k 方阵（与 Kyber 等的 k×k 一致）；
//       b_i = Σ_j A_{ij} * s_j + e_i。
void MLWEScheme::KeyGen(MLWEPublicKey& pk, MLWESecretKey& sk) const {
    usint k = m_ctx->GetParams().k;

    // --- 步骤 1：采样公共矩阵 A ∈ R_q^{k×k} ---
    pk.A.assign(k, std::vector<RingElement>(k));
    for (usint i = 0; i < k; ++i) {
        for (usint j = 0; j < k; ++j) {
            pk.A[i][j] = m_ctx->SampleUniform();
        }
    }

    // --- 步骤 2：采样私钥 s ∈ χ^k 与噪声 e ∈ χ^k ---
    sk.assign(k, m_ctx->MakeElement());
    std::vector<RingElement> e(k, m_ctx->MakeElement());
    for (usint i = 0; i < k; ++i) {
        sk[i] = m_ctx->SampleGaussian();
        e[i] = m_ctx->SampleGaussian();
    }

    // --- 步骤 3：计算 b = A·s + e ∈ R_q^k ---
    //   b_i = Σ_{j=0}^{k-1} A_{ij} * s_j + e_i
    // 注意：这里的多项式乘法 * 在 R_q 中进行（自动模 x^n+1 与 q）。
    pk.b.assign(k, m_ctx->MakeElement());
    for (usint i = 0; i < k; ++i) {
        RingElement acc = m_ctx->MakeElement();  // 累加器，初值 0
        for (usint j = 0; j < k; ++j) {
            // 累加 A_{ij} * s_j，每次都是一次 R_q 中的环乘 + 环加
            acc += pk.A[i][j] * sk[j];
        }
        pk.b[i] = acc + e[i];  // 加噪声 e_i
    }
}

//------------------------------------------------------------------------------
// Encrypt：加密
//------------------------------------------------------------------------------
// 算法（MLWE-PKE 的标准加密）：
//   输入：公钥 (A, b)、消息 m ∈ R_q
//   1. 采样小噪声 r ← χ^k；                         // 加密随机性
//   2. c1 = A^T · r        ∈ R_q^k；                // “掩码”部分
//   3. c0 = b^T · r + m    ∈ R_q；                  // 消息部分（含噪声）
//   4. 输出密文 (c0, c1)。
//
// 说明：c1 是 k 维向量，c1_i = Σ_j A_{ji} * r_j（注意是 A 的转置）。
MLWECiphertext MLWEScheme::Encrypt(const MLWEPublicKey& pk,
                                   const RingElement& m) const {
    usint k = m_ctx->GetParams().k;

    // --- 步骤 1：采样加密随机性 r ← χ^k ---
    std::vector<RingElement> r(k, m_ctx->MakeElement());
    for (usint i = 0; i < k; ++i) {
        r[i] = m_ctx->SampleGaussian();
    }

    MLWECiphertext ct;
    ct.c1.assign(k, m_ctx->MakeElement());

    // --- 步骤 2：c1 = A^T · r ---
    //   c1_i = Σ_{j=0}^{k-1} A_{j,i} * r_j     （A 的第 i 列与 r 内积）
    for (usint i = 0; i < k; ++i) {
        RingElement acc = m_ctx->MakeElement();
        for (usint j = 0; j < k; ++j) {
            acc += pk.A[j][i] * r[j];
        }
        ct.c1[i] = acc;
    }

    // --- 步骤 3：c0 = b^T · r + m ---
    //   c0 = Σ_{j=0}^{k-1} b_j * r_j + m
    RingElement c0 = m_ctx->MakeElement();
    for (usint j = 0; j < k; ++j) {
        c0 += pk.b[j] * r[j];
    }
    c0 += m;

    ct.c0 = c0;
    return ct;
}

//------------------------------------------------------------------------------
// Decrypt：解密
//------------------------------------------------------------------------------
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
//   n·σ^2），远小于 q/2，因此 m' 与 m 在 R_q 中逐系数接近。
RingElement MLWEScheme::Decrypt(const MLWECiphertext& ct,
                                const MLWESecretKey& sk) const {
    // 调用核心运算：c0 - s^T · c1
    return DecryptCore(ct.c0, ct.c1, sk);
}

//------------------------------------------------------------------------------
// DecryptCore：解密核心运算 c0 - s^T · c1
//------------------------------------------------------------------------------
//   s^T · c1 = Σ_{i=0}^{k-1} s_i * c1_i   （R_q 中的环乘累加）
RingElement MLWEScheme::DecryptCore(
    const RingElement& c0,
    const std::vector<RingElement>& c1,
    const MLWESecretKey& sk) const {
    usint k = m_ctx->GetParams().k;

    // 计算 s^T · c1
    RingElement inner = m_ctx->MakeElement();
    for (usint i = 0; i < k; ++i) {
        inner += sk[i] * c1[i];
    }
    // m' = c0 - s^T · c1
    return c0 - inner;
}

//------------------------------------------------------------------------------
// InnerProduct：两个 k 维 R_q 向量的内积 a^T · b
//------------------------------------------------------------------------------
RingElement MLWEScheme::InnerProduct(const std::vector<RingElement>& a,
                                     const std::vector<RingElement>& b) const {
    if (a.size() != b.size()) {
        throw std::invalid_argument("InnerProduct: 向量长度不一致。");
    }
    RingElement acc = m_ctx->MakeElement();
    for (size_t i = 0; i < a.size(); ++i) {
        acc += a[i] * b[i];
    }
    return acc;
}

//==============================================================================
// 工具函数实现
//==============================================================================

//------------------------------------------------------------------------------
// ElementToVector：把环元素的系数导出为 int64_t 向量
//------------------------------------------------------------------------------
std::vector<int64_t> ElementToVector(const RingElement& e) {
    // 切换到系数表示，逐项取值。
    RingElement eCoeff = e;
    eCoeff.SetFormat(lbcrypto::COEFFICIENT);
    auto n = eCoeff.GetLength();
    std::vector<int64_t> out(n);
    for (usint i = 0; i < n; ++i) {
        out[i] = eCoeff[i].ConvertToInt();
    }
    // 注意：OpenFHE 的系数返回的是 [0, q) 内的非负值。
    // 为得到“带符号”表示，我们对 > q/2 的值减去 q。
    // 但此处我们返回原始值，由调用方决定如何解读。
    return out;
}

//------------------------------------------------------------------------------
// VectorToElement：把 int64_t 系数向量打包成环元素
//------------------------------------------------------------------------------
RingElement VectorToElement(const std::vector<int64_t>& coeffs,
                            const MLWEContext& ctx) {
    RingElement e = ctx.MakeElement();
    for (size_t i = 0; i < coeffs.size() && i < ctx.GetParams().n; ++i) {
        e.SetValueAtIndex(i, coeffs[i]);
    }
    return e;
}

//------------------------------------------------------------------------------
// ElementToString：把环元素打印成字符串
//------------------------------------------------------------------------------
// all=false：只打印前若干系数 + 汇总统计；
// all=true ：打印全部系数。
std::string ElementToString(const RingElement& e, bool all) {
    std::vector<int64_t> v = ElementToVector(e);
    usint n = v.size();
    usint q = 0;
    // 获取模数（从环参数）
    // 这里通过 e.GetModulus() 拿到 NativeInteger，转 int64_t
    q = static_cast<usint>(e.GetModulus().ConvertToInt());

    // 把 [0,q) 表示转为带符号表示 [−⌊q/2⌋, ⌊q/2⌋)
    auto toSigned = [q](int64_t x) -> int64_t {
        if (x > (int64_t)q / 2) x -= (int64_t)q;
        return x;
    };

    // 统计：绝对值最大、L2 范数等
    int64_t maxAbs = 0;
    double l2 = 0.0;
    for (auto x : v) {
        int64_t s = toSigned(x);
        maxAbs = std::max(maxAbs, std::llabs(s));
        l2 += (double)s * s;
    }
    l2 = std::sqrt(l2);

    std::ostringstream oss;
    oss << "[mod q = " << q << ", n = " << n << "] ";

    if (all) {
        oss << "coeffs (signed): {";
        for (usint i = 0; i < n; ++i) {
            if (i) oss << ", ";
            oss << toSigned(v[i]);
        }
        oss << "}";
    } else {
        oss << "first coeffs (signed): {";
        usint show = std::min<usint>(n, 8);
        for (usint i = 0; i < show; ++i) {
            if (i) oss << ", ";
            oss << toSigned(v[i]);
        }
        oss << (n > show ? ", ..." : "");
        oss << "}";
    }
    oss << " |max|=" << maxAbs << ", ||·||_2=" << l2;
    return oss.str();
}

//------------------------------------------------------------------------------
// MaxCoeffAbsDiff：逐系数绝对差的最大值（衡量解密误差）
//------------------------------------------------------------------------------
// 这里用带符号表示来比较，以正确处理模 q 的环绕。
int64_t MaxCoeffAbsDiff(const RingElement& a, const RingElement& b) {
    std::vector<int64_t> va = ElementToVector(a);
    std::vector<int64_t> vb = ElementToVector(b);
    usint n = std::min(va.size(), vb.size());
    usint q = static_cast<usint>(a.GetModulus().ConvertToInt());

    auto toSigned = [q](int64_t x) -> int64_t {
        if (x > (int64_t)q / 2) x -= (int64_t)q;
        return x;
    };

    int64_t maxDiff = 0;
    for (usint i = 0; i < n; ++i) {
        int64_t d = std::llabs(toSigned(va[i]) - toSigned(vb[i]));
        // 考虑环绕：差值也可能是 q - d（当两个系数分别落在 0 和 q-1 附近时）
        d = std::min(d, (int64_t)q - d);
        maxDiff = std::max(maxDiff, d);
    }
    return maxDiff;
}

}  // namespace mlwe
