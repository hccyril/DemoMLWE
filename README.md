# 基于 OpenFHE 的 Module-LWE（MLWE）模块与演示程序

本项目从零实现了一个 **Module Learning With Errors（Module-LWE，模块化容错学习）**
公钥加密方案，直接复用 [OpenFHE](https://github.com/openfheorg/openfhe-development)
的底层多项式环原语（`NativePoly` / `ILParams2N` / 数论工具），而不依赖
OpenFHE 内置的 BFV/BGV/CKKS 等高级方案。

> 背景说明：OpenFHE 内置了 RLWE（Ring-LWE）类型的方案（BFV/BGV/CKKS 等），
> 但**没有**单独暴露通用的 MLWE 模块。MLWE 是介于标准 LWE 与 RLWE 之间的
> 一类格难题，是 **Kyber（ML-KEM）/ Dilithium（ML-DSA）** 等 NIST 后量子标准
> 算法的核心假设。本项目补齐了这一空白，并提供一个可读、可运行的参考实现。

---

## 目录

1. [MLWE 是什么](#1-mlwe-是什么)
2. [数学原理与公式](#2-数学原理与公式)
3. [项目结构](#3-项目结构)
4. [依赖与环境](#4-依赖与环境)
5. [编译与运行](#5-编译与运行)
6. [输出解读](#6-输出解读)
7. [关键参数说明](#7-关键参数说明)
8. [API 参考](#8-api-参考)
9. [扩展方向](#9-扩展方向)
10. [免责声明](#10-免责声明)

---

## 1. MLWE 是什么

**LWE / RLWE / MLWE 三者关系：**

| 方案 | 公共部分 | 秘密 / 噪声 | 维度 |
|------|----------|-------------|------|
| **LWE**  | 矩阵 `A ∈ Z_q^{m×n}` | 向量 `s, e ∈ Z_q^n` | 标量级 |
| **RLWE** | 环元素 `a ∈ R_q` | 环元素 `s, e ∈ R_q` | 单个环元素（k=1） |
| **MLWE** | 矩阵 `A ∈ R_q^{k×k}` | 向量 `s, e ∈ R_q^k` | k 个环元素 |

其中 `R_q = Z_q[x] / (x^n + 1)` 是分圆多项式商环。

- **RLWE** 是 **MLWE** 在 `k = 1` 时的特例。
- **LWE** 不使用环结构，效率较低；**RLWE** 速度最快但密钥尺寸固定；
- **MLWE** 通过调节 `k` 在 **安全强度 / 密钥尺寸 / 运算效率** 之间取得平滑折中，
  这正是 Kyber 选择它的原因。

**MLWE 的搜索 / 判定困难性：**

- **判定型 MLWE**：给定 `(A, b = A·s + e)`，无法区分 `b` 是由小噪声 `e` 生成
  还是均匀随机。
- **搜索型 MLWE**：给定 `(A, b = A·s + e)`，恢复小秘密 `s` 是困难的。

---

## 2. 数学原理与公式

### 2.1 记号

- `R = Z[x] / (x^n + 1)`，`n` 为 2 的幂（分圆多项式 `Φ_{2n}(x) = x^n + 1`）。
- `R_q = R / qR = Z_q[x] / (x^n + 1)`，模数 `q` 为素数，且 `q ≡ 1 (mod 2n)`
  （保证 `x^n+1` 在 `Z_q` 上完全分裂，NTT 可用）。
- 模块秩 `k`：秘密、噪声均为 `R_q` 上的 `k` 维向量。
- 噪声分布 `χ`：以 0 为中心、标准差 `σ` 的离散高斯（本项目用“舍入高斯”）。

### 2.2 密钥生成 KeyGen

```
输入：参数 (n, k, q, σ)
  1. 均匀采样 A ← R_q^{k×k}                       // 公共矩阵
  2. 从小噪声分布采样 s ← χ^k，e ← χ^k             // 秘密、噪声
  3. 计算 b = A·s + e ∈ R_q^k                      // 在 R_q 中
  4. 公钥 pk = (A, b)，私钥 sk = s
```

其中矩阵—向量乘法 `b_i = Σ_{j=0}^{k-1} A_{ij} · s_j + e_i`，每个 `·` 是
`R_q` 中的多项式乘法（自动模 `x^n+1` 与 `q`）。

### 2.3 加密 Encrypt（MLWE-PKE）

```
输入：公钥 (A, b)、消息 m ∈ R_q
  1. 采样小噪声 r ← χ^k
  2. c1 = A^T · r        ∈ R_q^k     // “掩码”部分
  3. c0 = b^T · r + m    ∈ R_q       // 消息部分（含噪声）
  4. 输出密文 (c0, c1)
```

### 2.4 解密 Decrypt

```
输入：密文 (c0, c1)、私钥 s
计算：m' = c0 - s^T · c1 ∈ R_q
```

### 2.5 正确性分析（解密误差为何“小”）

把 `b = A·s + e` 代入展开：

```
m' = c0 - s^T · c1
   = (b^T · r + m) - s^T · (A^T · r)
   = (A·s + e)^T · r + m - s^T · A^T · r
   = s^T · A^T · r + e^T · r + m - s^T · A^T · r
   = m + e^T · r
```

线性项 `s^T · A^T · r` 完全抵消，剩余 `e^T · r` 是两个小噪声多项式的内积，
其系数量级约 `n · σ^2`，**远小于 `q/2`**，因此 `m'` 与 `m` 在 `R_q` 中逐系数
接近。只要每个系数误差 `< q/2`，就能正确恢复明文。

### 2.6 加法同态性

MLWE 密文天然支持一次（或多次）加法同态：

```
若 (c0, c1) ↔ m，(c0', c1') ↔ m'，则
(c0 + c0', c1 + c1') ↔ m + m'
```

每次加法噪声累加，因此支持的加法次数受 `q / (n·σ^2)` 限制。

---

## 3. 项目结构

```
MLWE_Demo/
├── README.md                  # 本说明文档
├── CMakeLists.txt             # 构建脚本（CMake）
├── include/
│   └── mlwe.h                 # MLWE 模块头文件（含完整公式注释）
└── src/
    ├── mlwe.cpp               # MLWE 模块实现
    └── mlwe_demo.cpp          # 演示程序（main）
```

---

## 4. 依赖与环境

| 组件 | 版本要求 | 说明 |
|------|----------|------|
| **C++ 编译器** | C++17 及以上 | MSVC / GCC / Clang 均可 |
| **CMake** | ≥ 3.16 | 构建系统 |
| **OpenFHE** | v1.x（推荐 1.1+） | 需先编译安装 |

### 安装 OpenFHE（简述）

```bash
# Linux / macOS
git clone https://github.com/openfheorg/openfhe-development.git
cd openfhe-development
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/opt/openfhe
make -j$(nproc)
sudo make install
```

Windows 上可用 MSVC + vcpkg 或手动 CMake 构建，详见
[OpenFHE 官方文档](https://openfhe-development.readthedocs.io/)。

---

## 5. 编译与运行

### 5.1 编译

```bash
cd MLWE_Demo
mkdir build && cd build

# 指定 OpenFHE 安装路径（若未装在系统默认路径）
cmake .. -DOpenFHE_DIR=/opt/openfhe/lib/cmake/OpenFHE

# 构建（Release 模式，启用优化）
cmake --build . --config Release
```

### 5.2 运行

```bash
# Linux/macOS
./mlwe_demo

# Windows
.\Release\mlwe_demo.exe
```

程序无需任何命令行参数，会依次执行：
密钥生成 → 多组加密/解密正确性测试 → 误差分析 → 加法同态性演示。

---

## 6. 输出解读

程序典型输出（节选）如下：

```
========================================
  Module-LWE (MLWE) 演示程序 — 基于 OpenFHE
========================================
[参数] n(环维数)=64, k(模块秩)=2, q(模数)=7681, σ(噪声标准差)=4
[环]  R_q = Z_7681[x] / (x^64 + 1)

========================================
  第 2 步：密钥生成 KeyGen
========================================
采样 A ← R_q^{k×k}（均匀），s, e ← χ^k（离散高斯）
计算 b = A·s + e ∈ R_q^k
...

----------------------------------------
----- 测试: 常数消息 -----
  明文 m    : [mod q = 7681, n = 64] first coeffs (signed): {5, 0, ...} |max|=5, ||·||_2=5
  解密 m'   : [mod q = 7681, n = 64] first coeffs (signed): {5, 0, ...} |max|=5, ||·||_2=5
  >>> 解密最大逐系数误差 |m' - m|_∞ = 12
```

**关键观察：**

- 私钥 `s` 的系数非常小（集中在 0 附近，量级 ≤ 2σ）。
- 公钥 `b` 的系数接近均匀分布（因为 `A·s` 把小秘密“放大”到全 `Z_q`）。
- **解密误差** `|m' - m|_∞` 远小于 `q/2`（7681/2 ≈ 3840），因此解密成功。
- 加法同态的误差约为单次误差的 2 倍（两个噪声 `e^T·r` 相加）。

---

## 7. 关键参数说明

在 `src/mlwe_demo.cpp` 的 `main()` 中可调：

```cpp
usint n  = 64;     // 环维数（2 的幂；Kyber 用 256）
usint k  = 2;      // 模块秩
usint q  = 7681;   // 模数（素数，且需 q ≡ 1 (mod 2n)）
usint nu = 4;      // 噪声高斯标准差 σ（教学参数）
MLWEParams params(n, k, q, nu);
```

### 参数选取约束

| 参数 | 约束 | 备注 |
|------|------|------|
| `n` | 必须为 2 的幂 | 分圆多项式要求 |
| `q` | 素数，且 `q ≡ 1 (mod 2n)` | 否则 NTT 无效；程序会抛异常 |
| `k` | 任意正整数 | 越大越安全但越慢 |
| `nu` (σ) | 正数 | 越大噪声越大，但安全等级也变化 |

### 合法 q 的示例

| n | 合法 q（示例） |
|---|----------------|
| 64 | 7681, 12289 |
| 128 | 12289 |
| 256 | 7681, 12289, 3329 |

> 例如 `n=64, q=7681`：`7681 = 60 × 128 + 1`，素数，满足 `q ≡ 1 (mod 2n=128)`。

### 安全性提醒

本项目使用 **小参数** 仅为 **教学演示** 用，**不达到 NIST 后量子安全强度**。
真实部署应参照 Kyber 参数集（`n=256, k∈{2,3,4}, q=3329`），并采用
**中心化二项分布** 替代高斯分布作为噪声来源（抗侧信道更优）。

---

## 8. API 参考

### 8.1 `mlwe::MLWEParams` — 参数集

```cpp
struct MLWEParams {
    usint n;    // 环维数
    usint k;    // 模块秩
    usint q;    // 模数
    usint nu;   // 噪声标准差 σ
    usint base; // gadget 基（保留，本演示未使用）
    MLWEParams(usint n_=256, usint k_=2, usint q_=7681, usint nu_=8, usint base_=0);
};
```

### 8.2 `mlwe::MLWEContext` — 上下文

```cpp
class MLWEContext {
public:
    explicit MLWEContext(const MLWEParams& params);
    const MLWEParams& GetParams() const;
    const lbcrypto::ILParams2N& GetElementParams() const;

    RingElement MakeElement() const;              // 零元素
    RingElement MakeConstantElement(int64_t c) const;  // 常数元素 c·1
    RingElement SampleGaussian() const;           // 从 χ 采样环元素
    RingElement SampleUniform() const;            // 从 U(Z_q) 采样环元素
};
```

### 8.3 `mlwe::MLWEScheme` — 方案

```cpp
class MLWEScheme {
public:
    explicit MLWEScheme(std::shared_ptr<MLWEContext> ctx);

    void KeyGen(MLWEPublicKey& pk, MLWESecretKey& sk) const;
    MLWECiphertext Encrypt(const MLWEPublicKey& pk, const RingElement& m) const;
    RingElement Decrypt(const MLWECiphertext& ct, const MLWESecretKey& sk) const;

    RingElement InnerProduct(...) const;          // 向量内积 a^T·b
    RingElement DecryptCore(...) const;           // c0 - s^T·c1
};
```

### 8.4 工具函数

```cpp
std::string ElementToString(const RingElement& e, bool all = false);
std::vector<int64_t> ElementToVector(const RingElement& e);
RingElement VectorToElement(const std::vector<int64_t>& coeffs, const MLWEContext& ctx);
int64_t MaxCoeffAbsDiff(const RingElement& a, const RingElement& b);
```

---

## 9. 扩展方向

本项目是 MLWE 的最小可运行实现，可作为以下扩展的基础：

1. **消息缩放与模舍入（modulus rounding）**：实现 `Δ = ⌊q/t⌋` 的 BFV 风格消息
   编码，使消息空间为 `R_t`（`t < q`），解密后做模 `t` 还原精确明文。
2. **NTT 域优化**：把所有元素预先 `SetFormat(EVALUATION)`，避免反复 NTT/INTT。
3. **gadget 分解与 G 矩阵**：实现 `G = I_k ⊗ g^T`，支持更高效的密钥切换。
4. **Regev 风格的密文压缩 / 公钥压缩**：模拟 Kyber 的 bit-packing。
5. **CPA→CCA 转换（Fujisaki-Okamoto）**：把本 MLWE-PKE 升级为 MLWE-KEM。
6. **中心化二项分布（CBD）**：替换高斯噪声，与 Kyber 对齐，抗侧信道更优。

---

## 10. 免责声明

- 本实现 **仅供学习与研究** 使用，**未经密码学安全审计**，**不可用于生产环境**。
- 默认参数为教学小参数，**不提供任何安全保证**。
- 如需生产级 MLWE，请使用经过审计的库，如
  [liboqs](https://github.com/open-quantum-safe/liboqs) 中的 Kyber/Dilithium 实现。

---

## 参考

- Langlois, A., & Stehlé, D. (2015). *Worst-case to average-case reductions for module lattices.*
- Bos, J., et al. (2018). *CRYSTALS – Kyber: a CCA-secure module-lattice-based KEM.*
- [OpenFHE 官方仓库与文档](https://github.com/openfheorg/openfhe-development)
- [OpenFHE API 文档](https://openfhe-development.readthedocs.io/)
- [NIST PQC 标准（ML-KEM / ML-DSA）](https://csrc.nist.gov/projects/post-quantum-cryptography)
