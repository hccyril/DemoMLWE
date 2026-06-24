//==================================================================================
// mlwe_demo.cpp
//
// MLWE 演示程序（基于 OpenFHE）。
//
// 本演示完整地走一遍 MLWE 公钥加密方案的流程：
//   1. 构造 MLWE 参数与上下文；
//   2. 密钥生成 KeyGen：得到公钥 (A, b=A·s+e) 与私钥 s；
//   3. 加密 Encrypt：对消息 m 加密得到 (c0, c1)；
//   4. 解密 Decrypt：用私钥恢复 m'，验证 m' ≈ m；
//   5. 误差分析：打印解密误差 |m' − m| 的统计；
//   6. 多组测试：分别加密“常数消息”、“带符号消息”、“随机消息”进行验证；
//   7. 加法同态性：演示 (c0+c0', c1+c1') 解密 ≈ m + m'。
//
// 编译与运行请参见 README.md 与 CMakeLists.txt。
//==================================================================================

#include "mlwe.h"

#include <iostream>
#include <string>
#include <vector>

using namespace mlwe;

//------------------------------------------------------------------------------
// 工具：打印分割线
//------------------------------------------------------------------------------
static void PrintBar(const std::string& title) {
    std::cout << "\n========================================"
                 "========================================\n";
    std::cout << "  " << title << "\n";
    std::cout << "========================================"
                 "========================================\n";
}

//------------------------------------------------------------------------------
// 工具：把消息从带符号 int64 向量打包为 RingElement，并显示
//------------------------------------------------------------------------------
static RingElement PackAndShow(const std::vector<int64_t>& coeffs,
                               const MLWEContext& ctx,
                               const std::string& name) {
    RingElement m = VectorToElement(coeffs, ctx);
    std::cout << "[消息 " << name << "] " << ElementToString(m, true) << "\n";
    return m;
}

//------------------------------------------------------------------------------
// 单次加密-解密正确性验证
//------------------------------------------------------------------------------
// 返回：解密误差的最大绝对差（逐系数，带符号）。
static int64_t RunEncDecap(const MLWEScheme& scheme,
                           const MLWEPublicKey& pk,
                           const MLWESecretKey& sk,
                           const RingElement& m,
                           const std::string& tag) {
    std::cout << "\n----- 测试: " << tag << " -----\n";

    // 加密
    MLWECiphertext ct = scheme.Encrypt(pk, m);

    // 解密
    RingElement m_prime = scheme.Decrypt(ct, sk);

    // 打印解密结果（仅前 8 个系数，避免刷屏）
    std::cout << "  明文 m    : " << ElementToString(m, false) << "\n";
    std::cout << "  解密 m'   : " << ElementToString(m_prime, false) << "\n";

    // 计算并打印逐系数最大误差
    int64_t err = MaxCoeffAbsDiff(m, m_prime);
    std::cout << "  >>> 解密最大逐系数误差 |m' - m|_∞ = " << err << "\n";
    return err;
}

int main() {
    try {
        PrintBar("Module-LWE (MLWE) 演示程序 — 基于 OpenFHE");

        // --------------------------------------------------------------
        // 第 1 步：设置参数并构造上下文
        // --------------------------------------------------------------
        // 这里选用一组“小而可读”的教学参数：
        //   n   = 64     环维数（2 的幂；Kyber 用 256，这里为加速用 64）
        //   k   = 2      模块秩
        //   q   = 7681   模数（素数，且 7681 ≡ 1 (mod 2*64=128) 成立）
        //   nu  = 4      噪声高斯标准差 σ（教学参数）
        //
        // 安全性说明：这组参数仅为演示用，并不达到 NIST 级别安全强度。
        // 真实部署应采用 n=256, k≥2, q≥3329, 并配合中心化二项分布（Kyber 风格）。
        usint n = 64;
        usint k = 2;
        usint q = 7681;   // 7681 = 60*128 + 1, 素数, 满足 q ≡ 1 (mod 2n)
        usint nu = 4;     // σ = 4.0
        MLWEParams params(n, k, q, nu);

        std::cout << "[参数] n(环维数)=" << params.n
                  << ", k(模块秩)=" << params.k
                  << ", q(模数)=" << params.q
                  << ", σ(噪声标准差)=" << params.nu << "\n";
        std::cout << "[环]  R_q = Z_" << params.q
                  << "[x] / (x^" << params.n << " + 1)\n";

        auto ctx = std::make_shared<MLWEContext>(params);

        // --------------------------------------------------------------
        // 第 2 步：密钥生成
        // --------------------------------------------------------------
        PrintBar("第 2 步：密钥生成 KeyGen");
        std::cout << "采样 A ← R_q^{k×k}（均匀），s, e ← χ^k（离散高斯）\n";
        std::cout << "计算 b = A·s + e ∈ R_q^k\n";

        MLWEScheme scheme(ctx);
        MLWEPublicKey pk;
        MLWESecretKey sk;
        scheme.KeyGen(pk, sk);

        std::cout << "公钥 b[0] (前 8 系数): "
                  << ElementToString(pk.b[0], false) << "\n";
        std::cout << "私钥 s[0] (前 8 系数): "
                  << ElementToString(sk[0], false) << "\n";
        std::cout << "（注意：s 的系数很小，集中在 0 附近；b 的系数接近均匀）\n";

        // --------------------------------------------------------------
        // 第 3 步：构造若干测试消息
        // --------------------------------------------------------------
        PrintBar("第 3 步：构造测试消息");

        // 消息 1：全零消息 m = 0 ∈ R_q
        std::vector<int64_t> zeroMsg(params.n, 0);
        RingElement m_zero = PackAndShow(zeroMsg, *ctx, "全零 m=0");

        // 消息 2：常数消息 m = 5·1
        std::vector<int64_t> constMsg(params.n, 0);
        constMsg[0] = 5;
        RingElement m_const = PackAndShow(constMsg, *ctx, "常数 m=5");

        // 消息 3：带符号小消息（前 4 项为 ±1）
        std::vector<int64_t> signMsg(params.n, 0);
        signMsg[0] = 1; signMsg[1] = -1; signMsg[2] = 1; signMsg[3] = -1;
        RingElement m_sign = PackAndShow(signMsg, *ctx, "带符号 m=(1,-1,1,-1,0,...)");

        // 消息 4：随机较大消息（前 4 项在 [0, 100] 内）
        //   说明：当消息系数较大时，可能超过噪声容限导致解密失败——
        //   这正是 MLWE 消息需要被“限制在小范围”或“做模 q 缩放”的原因。
        std::vector<int64_t> randMsg(params.n, 0);
        randMsg[0] = 42; randMsg[1] = 7; randMsg[2] = 99; randMsg[3] = 13;
        RingElement m_rand = PackAndShow(randMsg, *ctx, "随机 m=(42,7,99,13,0,...)");

        // --------------------------------------------------------------
        // 第 4 步：加密-解密正确性验证
        // --------------------------------------------------------------
        PrintBar("第 4 步：加密 Encrypt -> 解密 Decrypt 正确性验证");

        int64_t err_zero  = RunEncDecap(scheme, pk, sk, m_zero,  "全零消息");
        int64_t err_const = RunEncDecap(scheme, pk, sk, m_const, "常数消息");
        int64_t err_sign  = RunEncDecap(scheme, pk, sk, m_sign,  "带符号消息");
        int64_t err_rand  = RunEncDecap(scheme, pk, sk, m_rand,  "随机消息");

        // --------------------------------------------------------------
        // 第 5 步：误差分析汇总
        // --------------------------------------------------------------
        PrintBar("第 5 步：解密误差分析");
        usint q = ctx->GetParams().q;
        std::cout << "模数 q = " << q << "，解密容限 ≈ q/2 = " << (q / 2) << "\n";
        std::cout << "（只要逐系数误差 < q/2，解密即正确）\n\n";
        std::cout << "  全零消息   最大误差: " << err_zero  << "\n";
        std::cout << "  常数消息   最大误差: " << err_const << "\n";
        std::cout << "  带符号消息 最大误差: " << err_sign  << "\n";
        std::cout << "  随机消息   最大误差: " << err_rand  << "\n";

        // 综合判定
        int64_t maxErr = std::max({err_zero, err_const, err_sign, err_rand});
        std::cout << "\n  综合最大误差 = " << maxErr
                  << "  (q/2 = " << (q / 2) << ")\n";
        if (maxErr < (int64_t)(q / 2)) {
            std::cout << "  >>> 结论：所有测试用例解密成功（误差 < q/2）。✓\n";
        } else {
            std::cout << "  >>> 结论：部分用例超出容限（增大 σ 会增加失败概率）。\n";
        }

        // --------------------------------------------------------------
        // 第 6 步：加法同态性演示
        // --------------------------------------------------------------
        // MLWE 密文天然支持加法同态：
        //   若 (c0, c1) 加密 m，(c0', c1') 加密 m'，则
        //   (c0+c0', c1+c1') 加密 m + m'（噪声也随之累加）。
        //
        // 验证：解密 (c0+c0', c1+c1') 应得到 ≈ m + m'。
        PrintBar("第 6 步：加法同态性演示 (m + m') ");

        MLWECiphertext ct1 = scheme.Encrypt(pk, m_const);   // m = 5
        MLWECiphertext ct2 = scheme.Encrypt(pk, m_sign);    // m' = (1,-1,1,-1,..)
        // 明文和（用于对比）
        RingElement m_sum = m_const + m_sign;
        std::cout << "明文和  m + m'      : " << ElementToString(m_sum, false) << "\n";

        // 密文相加（逐元素环加）
        MLWECiphertext ct_sum;
        ct_sum.c0 = ct1.c0 + ct2.c0;
        ct_sum.c1.resize(ctx->GetParams().k, ctx->MakeElement());
        for (usint i = 0; i < ctx->GetParams().k; ++i) {
            ct_sum.c1[i] = ct1.c1[i] + ct2.c1[i];
        }
        // 同态解密
        RingElement m_sum_prime = scheme.Decrypt(ct_sum, sk);
        std::cout << "同态解密 (m+m')'    : "
                  << ElementToString(m_sum_prime, false) << "\n";
        int64_t err_hom = MaxCoeffAbsDiff(m_sum, m_sum_prime);
        std::cout << ">>> 加法同态误差    : " << err_hom
                  << "（约为单次误差的 2 倍，因噪声累加）\n";

        // --------------------------------------------------------------
        // 第 7 步：总结
        // --------------------------------------------------------------
        PrintBar("演示完成");
        std::cout << "本演示展示了 MLWE 的核心操作：\n"
                  << "  • 公钥/私钥生成（基于 R_q 上的均匀矩阵 A 与高斯噪声 s, e）\n"
                  << "  • 加密（采样 r，计算 c0=b^T·r+m, c1=A^T·r）\n"
                  << "  • 解密（m'=c0 - s^T·c1 ≈ m + e^T·r）\n"
                  << "  • 解密误差分析（误差项 e^T·r 为小量）\n"
                  << "  • 加法同态性（密文相加对应明文相加）\n\n"
                  << "MLWE 是 Kyber(ML-KEM) / Dilithium(ML-DSA) 等后量子标准的核心假设。\n";
    } catch (const std::exception& ex) {
        std::cerr << "\n[错误] " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
