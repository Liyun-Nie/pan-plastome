# 技术说明：四类 pan-plastome 分析

本文介绍四类脚本的**设计、流程与参考文献**。运行命令和输入要求见根目录 [README.zh-CN.md](../README.zh-CN.md)。英文版：[METHODS.md](METHODS.md)。

下文描述以本仓库中的脚本实现为准。

---

## 0. 总览

| 模块 | 科学问题 | 核心推断 |
| --- | --- | --- |
| 有阈值 Hotspot | 正频率位点的高频尾部是否跨物种反复成为 Hotspot？ | 物种分层置换 + BH |
| 无阈值 Hotspot | 完整种内排名是否呈现共同高秩结构？ | 物种内标签置换；`Q_high` / `C_global` / `M_max`；maxT |
| Tajima’s D | CDS / IGS 频率谱是否偏离中性预期？ | Tajima (1989)，单倍体 *n* |
| 结构筛查 v4.1 | 短读长 SNV/indel 流程可能漏检的宏观重排？ | 比对优先的共线性区块筛查 |

**共同原则**

1. 两套 Hotspot 分析都先在**每个物种内部**建立随机零分布，再比较同一位点在不同物种中的表现。  
2. 有阈值和无阈值分析回答的问题不同，应分别解释两套分析的结果。  
3. 结构筛查使用 **single-IR 线性**坐标，分析这些坐标能够表示的重排。  
4. 完整 Hotspot 分析使用 **99,999** 次置换。仓库存放脚本，输入基因组和置换生成的数据保存在本地。

---

## 1. 有阈值 Hotspot（`hotspot_threshold_sensitivity`）

### 1.1 设计

Hotspot 定义为：每个物种**已发生变异**的候选集（`frequency_per_kb > 0`）中，频率最高的 Top *q* 尾部（*q* ∈ {10%, 15%, 20%, 25%, 30%}）。

- CDS 仅计 SNV；IGS 计 SNV + indel + CPX；均换算为每千碱基频率。  
- 排除零频率是这套分析的定义：零表示该物种在该位点没有检出变异。如果把大量零值纳入分位数计算，阈值可能降到零，从而把零频率位点错误地标为 Hotspot。  
- 位点被分为 Hotspot 和非 Hotspot 后，跨物种检验统计**同一位点被反复标为 Hotspot**的次数。

### 1.2 流程

1. 读取 `--workspace` 下的会话排名表（路径见 README）。  
2. 对每个物种、每个 *q*，在正频率子集上用 R `quantile(..., type = 7)`（Hyndman–Fan 第 7 型）计算阈值；`frequency ≥ threshold` 标为 Hotspot。阈值处并列全部保留，实际比例可以略高于名义 Top *q*。  
3. 观察复现次数 \(X_l=\sum_s H_{sl}\)。  
4. **物种分层置换**（主推断）：保持物种 *s* 的正频率候选集 \(E_s\) 与 Hotspot 数 \(h_s\)，在 \(E_s\) 中无放回抽取标签；物种间不交换；重复 *B* 次；上尾  
   \(P_l=(1+b_l)/(1+B)\)（Phipson–Smyth 形式）。正式 *B* = 99,999。  
5. 在每个「区域类型 × 阈值」族内做 Benjamini–Hochberg 校正。  
入口：`scripts/hotspot_threshold_sensitivity/hotspot_threshold_sensitivity.R`。

### 1.3 参考文献

1. Hyndman, R.J. & Fan, Y. Sample quantiles in statistical packages. *Am. Stat.* **50**, 361–365 (1996).  
2. Strasser, H. & Weber, C. On the asymptotic theory of permutation statistics. *Math. Methods Stat.* **8**, 220–250 (1999).  
3. Phipson, B. & Smyth, G.K. Permutation P-values should never be zero. *Stat. Appl. Genet. Mol. Biol.* **9**, Article 39 (2010).  
4. Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate. *J. R. Stat. Soc. B* **57**, 289–300 (1995).  
---

## 2. 无阈值 Hotspot（`hotspot_threshold_free_analysis`）

### 2.1 设计

不预先规定 Top 25%（或其它分位）才叫 Hotspot，而使用每个物种**全部可分析位点**的完整频率排名（含真实零频率）。缺失保持为缺失，不改写成零。

1. 种内按频率升序赋秩（并列用平均秩，零频率并列同样处理）。  
2. Rankit 百分位 \(u=(r-0.5)/n\)，再 \(z=\Phi^{-1}(u)\)。  
3. 单位点跨物种得分 \(T_l=\sum_s z_{sl}/\sqrt{K_l}\)（\(K_l\) 为覆盖物种数）。除以 \(\sqrt{K_l}\) 用于统一尺度。P 值由物种内置换得到。  
4. 三个整体统计量描述互补的模式：  
   - \(Q_{\mathrm{high}}=\sum_l[\max(T_l,0)]^2\)：主要整体统计量，看一批位点是否共同偏向高频端；  
   - \(C_{\mathrm{global}}\)：按重叠位点加权的种内百分位 Pearson 相关；  
   - \(M_{\mathrm{max}}=\max_l T_l\)：是否至少存在一个极端共同高秩位点。  
   三者使用同一套种内排名。  
5. 零模型：在每个物种内打乱“得分属于哪个位点”（保留频率多重集、零值并列与缺失结构）。  
6. 逐位点显著性：Westfall–Young 单步 maxT（对照置换最大值分布）。边际 P 的 BH 校正仅作辅助对照。

### 2.2 流程

1. 读取与有阈值相同的排名表。SHA-256 与行数门控来自 YAML（SHA 可回退到 helpers 中的常量）。  
2. 种内秩 → \(u\) → \(z\)；计算 \(T_l\)、\(Q_{\mathrm{high}}\)、\(C_{\mathrm{global}}\)、\(M_{\mathrm{max}}\)。  
3. 主置换 99,999 次。稳健性（较少次数；零值广度、正频率强度、长度残差、中心化秩、属级下采样）仅当 YAML 写明 `run_robustness: true` 时运行。  
4. 写出全局表、逐位点表、贡献分解、成对一致性等。

入口：`scripts/hotspot_threshold_free_analysis/hotspot_threshold_free_analysis.R`。

### 2.3 参考文献

1. Bliss, C.I. *Statistics in Biology* Vol. 1 (McGraw-Hill, 1967).（Rankit）  
2. Beasley, T.M., Erickson, S. & Allison, D.B. Rank-based inverse normal transformations… *Behav. Genet.* **39**, 580–595 (2009).  
3. Westfall, P.H. & Young, S.S. *Resampling-based Multiple Testing* (Wiley, 1993).  
4. Phipson & Smyth (2010)；Strasser & Weber (1999)。  
5. Benjamini & Hochberg (1995) — 辅助 BH；单位点主校正为 maxT。  
6. Spearman, C. The proof and measurement of association between two things. *Am. J. Psychol.* **15**, 72–101 (1904). — 秩相关为 \(C_{\mathrm{global}}\) 的背景；本实现计算全物种内百分位在重叠位点上的相关。

`Q_high`、`C_global` 和 `M_max` 是本分析定义的名称。

---

## 3. Tajima’s D（`selection_test_results`）

### 3.1 设计

比较成对差异均值 \(\pi\) 与基于分离位点数的 Watterson 估计 \(\theta_W=S/a_1\)。在中性、恒定群体大小下期望相等；\(D=(\pi-\theta_W)/\sqrt{\widehat{\mathrm{Var}}}\) 的偏离可以来自选择，也可以来自群体历史。

单倍体质体实现：

- \(n\) = accession 数，即 `{物种}/sample.txt` 的非空行数；  
- \(S\) = 策展变异表中该基因 / 区段的唯一 SNP 坐标数；  
- \(\pi=\sum_p 2j(p)\,(n-j(p))/(n(n-1))\)，\(j(p)\) 为该位点携带交替等位基因的唯一样本数（未出现在变异表中的样本按参考等位基因计）；  
- \(\theta\) 的生物学含义为 \(2N\mu\)，但 **D 的计算公式本身不含倍性系数**，只依赖 \(n\) 与 \(S\)；  
- Tajima (1989) 在无重组模型下推导方差；本实现对质体数据使用 Table 2 临界值；  
- 显著性：按 SNPGenie 方式查 Tajima (1989) Table 2。临界值不对称。样本量按区间匹配（例如 \(n\ge 50\) 使用 \(n=50\) 行）。\(S<3\) 或 \(n<4\) 的基因 / 区段跳过。

\(D<0\) 表示低频变异偏多，可能与纯化选择或群体扩张有关；\(D>0\) 表示中间频率变异偏多，可能与平衡选择或群体收缩有关。分析时结合生物学与群体历史证据解释这些可能性。

CDS：每个物种 × 基因一个 *D*。IGS：每个间隔区一个 *D*，并另报 IGS / CDS / 全部 SNP 合并后的全基因组口径 *D*（后者用 `--genome-regions` 中的 single-IR 长度归一）。合并表采用全基因组汇总口径，CDS 表采用按基因汇总口径。

### 3.2 流程

1. 读入变异总表、长度 / 区段表、样本目录。  
2. 按基因或 IGS 聚合分离位点，计算 \(S\)、\(\pi\)、\(D\) 与显著性档。  
3. 写出汇总表；可选 `visualize_tajima_d.py` 作 CDS 图。

入口：`selection_test_tajima.py`、`selection_test_tajima_igs.py`。两组必需参数不同，见 README。

### 3.3 参考文献

1. Tajima, F. Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics* **123**, 585–595 (1989).  
2. SNPGenie `Tajima_D.R`（Chase W. Nelson）— Table 2 查找与系数链的实现对照。

---

## 4. 结构筛查 v4.1（`synteny_analysis`）

### 4.1 设计

短读长 SNV/indel 流程可能漏掉较大的倒位和易位。本模块先进行全长序列比对，再筛查共线性变化，作为短读长分析的补充：

1. 每个样本相对物种参考做全长比对。主路径：minimap2 `asm5`。交叉核对：代表样本 blastn。  
2. 双轴共线性区块合并（`collinear_block_merge.py`）：同链相邻片段在查询轴与参考轴的间隙均落入 \([-\texttt{overlap_tol}, \texttt{max_gap}]\) 时合并。在主参数下检测 **区块间** 倒位 / 易位。  
3. 基于 PAF 缓存的参数网格、方法一致性表、PAF 图、链内局部异常，以及代表样本 nucmer 复核。  
4. 主证据链使用比对区块。`window_rearrangement_test.py` 和 `audit_window_skips.py` 提供其它模块导入的辅助函数。

坐标系：single-IR **线性** FASTA，与变异检测参考一致；检测范围为这些线性坐标能够表示的事件。

**SyRI。** 当 `PATH` 中可以找到 `syri` 时，`syri_representative_screen.py` 会运行 SyRI；否则改用 nucmer 和 `show-coords`。在 **nucmer 回退路径**中，反向比对块记为倒位证据，多个不共线的比对簇记为重排候选。

常用主参数：`min_len=500`、`max_gap=200`、`overlap_tol=50`、`translocation_jump=10000`（具体值在本地 YAML 中设置）。事件汇总报告 single-IR 线性比对中检测到的结构模式。

### 4.2 流程（概念）

```
FASTA + GenBank 参考
  → synteny_block_analysis.py
  → aggregate_synteny_results.py
  → （可选）敏感性 / BLAST 一致性 / PAF 图 / 局部异常
  → synteny_report_orchestrator.py
  → （可选）syri_representative_screen.py   # 有 SyRI 则用，否则 nucmer
```

命令及哪些脚本接受 `--config` 见 README。部分后续脚本默认读取本克隆下的结果目录，除非改路径参数。

### 4.3 参考文献

1. Li, H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* **34**, 3094–3100 (2018).  
2. Camacho, C. *et al.* BLAST+: architecture and applications. *BMC Bioinformatics* **10**, 421 (2009).  
3. Kurtz, S. *et al.* Versatile and open software for comparing large genomes. *Genome Biol.* **5**, R12 (2004).（MUMmer / nucmer）  
4. Goel, M. *et al.* SyRI: finding genomic rearrangements and local sequence differences from whole-genome assemblies. *Genome Biol.* **20**, 277 (2019). — 可选复核工具。

合并规则与事件定义以 `collinear_block_merge.py`、`synteny_block_analysis.py` 的实现为准。

---

## 5. 与 cpopvar 的关系

| 对象 | 位置 | 说明 |
| --- | --- | --- |
| 下游 Hotspot 分析（分层置换 / 连续秩 / maxT） | **本仓库** | 由上文所列脚本实现 |
| 排名表生成（预处理 / 归一化） | cpopvar 会话输出 | 作为本仓库脚本的输入 |

---

## 6. 复现注意

1. 基因组、PAF 缓存、结果 TSV 和 YAML 文件保存在本地。  
2. R 端 Hotspot 读取扁平 YAML 键；无阈值分析读取 `expected_*` 校验键。结构筛查的 Python 端读取可由 `yaml.safe_load` 解析为映射的 YAML 文件。  
3. 完整 Hotspot 分析使用 99,999 次置换；辅助测试和短置换用于代码检查与开发。  
4. 微型测试覆盖 Hotspot helpers 和 `test_collinear_block_merge.py`；输入校验使用 cpopvar 风格树中的本地排名表。

---

## 7. 文档维护

- 用法与命令：根 `README.zh-CN.md` / `README.md`  
- 方法细节：本文 `METHODS.zh-CN.md` / `METHODS.md`  
- 维护者：Liyun Nie \<nieliyun18@163.com\> · ORCID 0000-0002-5288-0041
