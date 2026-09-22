# pan-plastome

叶绿体 pan-plastome 研究配套脚本，包括**超变区（Hotspot）分析**、**Tajima’s D** 和**基于序列比对的结构筛查**。

| | |
| --- | --- |
| 维护者 | Liyun Nie \<nieliyun18@163.com\> |
| ORCID | [0000-0002-5288-0041](https://orcid.org/0000-0002-5288-0041) |
| 许可证 | MIT（2026） |
| 相关 R 包 | [Liyun-Nie/cpopvar](https://github.com/Liyun-Nie/cpopvar) |

这些脚本可用于了解分析方法，也可配合用户**自备**的排名表和基因组组装重新运行分析。R 包 `cpopvar` 负责变异处理、排名表生成和交互式工作流；本仓库提供下游 Hotspot 分析，包括物种分层置换，以及无阈值连续秩与 Westfall–Young maxT。

本仓库发布分析脚本。结果表、YAML 配置文件、质体 FASTA/GenBank 文件和序列数据由用户在本地准备。数据获取方式见论文的数据可用性说明，方法细节见 [`docs/METHODS.zh-CN.md`](docs/METHODS.zh-CN.md)。

**English**：[README.md](README.md) · [docs/METHODS.md](docs/METHODS.md)

---

## 目录

```
pan-plastome/
├── README.md
├── README.zh-CN.md
├── LICENSE
├── docs/
│   ├── METHODS.md
│   └── METHODS.zh-CN.md
└── scripts/
    ├── hotspot_threshold_sensitivity/
    ├── hotspot_threshold_free_analysis/
    ├── selection_test_results/
    └── synteny_analysis/
```

| 目录 | 用途 |
| --- | --- |
| `hotspot_threshold_sensitivity/` | 有阈值 Hotspot（正频率 Top 10–30%）；物种分层置换 |
| `hotspot_threshold_free_analysis/` | 无阈值连续秩；`Q_high` / `C_global` / `M_max`；maxT |
| `selection_test_results/` | CDS 与 IGS Tajima’s D（可选作图） |
| `synteny_analysis/` | 比对优先的共线性区块筛查（v4.1） |

---

## 运行环境

| 模块 | 运行时 | 常用依赖 |
| --- | --- | --- |
| Hotspot（两套） | R ≥ 4.1 | 辅助测试一般只需 base R |
| Tajima | Python 3 | `pandas`、`numpy`；作图需 `matplotlib`、`seaborn`、`scipy` |
| 结构筛查 | Python 3 | `pandas`、`PyYAML`、`biopython`；比对器 `minimap2`；可选 `blastn`、`nucmer` / `show-coords`、SyRI |

完整分析使用原始输入数据和本地配置文件。下文的小型测试用于检查部分函数；完整分析使用 99,999 次置换。

---

## 1. 有阈值 Hotspot

**用途：** 在每个物种内，把*非零*变异频率中最高的 10%、15%、20%、25% 和 30% 标为 Hotspot。CDS 只统计 SNV；IGS 统计 SNV、indel 和 CPX，并换算为每千碱基频率。随后检验同一位点在多个物种中反复成为 Hotspot 的次数是否高于随机预期。随机化过程会保留每个物种可参与分析的位点及 Hotspot 数量。

**入口：** `scripts/hotspot_threshold_sensitivity/hotspot_threshold_sensitivity.R`  
**辅助 / 测试。** `hotspot_sensitivity_helpers.R`、`test_hotspot_threshold_sensitivity.R`

### 需要自备的输入

将 `--workspace` 设置为 **cpopvar 风格的树根**。脚本会拼接：

```
<workspace>/app_data/sessions/<session_id>/results/plots/M03_hotspot/M03_hotspot_standard/M03_all_genes_ranked_with_thresholds.csv
<workspace>/app_data/sessions/<session_id>/results/plots/M04_poigs_hotspot_engine/M04_igs_hotspot_standard/M04_all_poigs_with_thresholds.csv
```

脚本会在 `--workspace` 后自行补上 `app_data/sessions/<session_id>`。上述两个结果目录中的候选 Hotspot CSV 文件也需要同时提供。请把 YAML 中的 `session_id` 设置为自己的会话目录名。

**配置：** 用 `--config` 传入扁平 YAML。嵌套块会被忽略。若未指定文件，脚本会找本克隆旁的 `config/hotspot_threshold_sensitivity.yml`；文件不存在时当作空配置，改用代码内默认值。实际使用仍建议自备 YAML，以便记录 `session_id` 与置换次数。示例键：

```yaml
session_id: your_session_id
n_perm: 99999
seed: 20251104
alpha: 0.05
quantile_type: 7
```

命令行可覆盖（文档以连字符为准，`--n_perm` 也能被解析）：`--workspace`、`--output`、`--config`、`--n-perm`、`--seed`、`--skip-perm`。

### 运行

```bash
# 在仓库根目录
Rscript scripts/hotspot_threshold_sensitivity/hotspot_threshold_sensitivity.R \
  --config /path/to/hotspot_threshold_sensitivity.yml \
  --workspace /path/to/cpopvar_tree_root \
  --output /path/to/output_dir

# 辅助测试（合成数据检查；仅当本克隆下有排名表时才做 25% 回归）
Rscript scripts/hotspot_threshold_sensitivity/test_hotspot_threshold_sensitivity.R
```

`--skip-perm` 会运行预处理和汇总步骤，并跳过置换阶段。

---

## 2. 无阈值 Hotspot

**用途：** 这套分析不设置 Top-% 界线，而是使用完整的种内排名。真实的零频率会保留，缺失值仍按缺失处理。脚本把排名转换为 Rankit 逆正态得分，再在每个物种内部置换位点标签，检验不同物种是否共享高排名位点。`Q_high` 是主要整体统计量，`C_global` 和 `M_max` 用于描述互补的模式；单位点检验采用 Westfall–Young 单步 maxT（`p_maxT`）。只有至少覆盖两个物种的位点才会进入分析。

**入口：** `scripts/hotspot_threshold_free_analysis/hotspot_threshold_free_analysis.R`  
**辅助 / 测试：** `hotspot_threshold_free_helpers.R`（会载入有阈值 helpers）、`test_hotspot_threshold_free_analysis.R`

`--workspace` 含义与第 1 节相同，排名表路径也相同。

### 完整运行所需配置

通过 `--config` 提供扁平 YAML 文件。完整运行会从该文件读取数据集专用的校验键。

配置文件应包括运行参数，以及根据用户自己的排名表计算出的校验值：

```yaml
session_id: your_session_id
run_mode: official
n_perm: 99999
n_perm_dev: 9999
n_perm_robust: 9999
seed: 12345
alpha: 0.05
chunk_size: 250
k_min_cross_species: 2
min_overlap_cglobal: 3
run_robustness: true
```

此外还需提供脚本使用的 `expected_*` 校验键，包括 CDS 和 IGS 的行数、非零记录数、SHA-256 值及 Top 25% 汇总计数。这些值随数据集变化，请根据自己的排名表计算。

**命令行可覆盖**（以连字符为准，下划线别名也可）：`--workspace`、`--output`、`--config`、`--n-perm`、`--n-perm-robust`、`--seed`、`--mode`、`--chunk-size`、`--skip-perm`、`--skip-robustness`。稳健性置换仅在 YAML 写明 `run_robustness: true` 且未加 `--skip-robustness` 时运行。`--mode dev` 使用 `n_perm_dev`，除非同时指定 `--n-perm`。

```bash
Rscript scripts/hotspot_threshold_free_analysis/hotspot_threshold_free_analysis.R \
  --config /path/to/hotspot_threshold_free.yml \
  --workspace /path/to/cpopvar_tree_root \
  --output /path/to/output_dir

Rscript scripts/hotspot_threshold_free_analysis/test_hotspot_threshold_free_analysis.R
```

辅助测试覆盖秩与置换代数。当克隆目录下存在排名表时，测试还会运行 SHA / 25% 回归检查。完整分析使用 99,999 次置换。

---

## 3. Tajima’s D（CDS / IGS）

**用途：** 检验 CDS 和 IGS 的位点频率谱是否偏离中性预期。质体按单倍体处理，因此 **n 是样本数**，即各物种 `sample.txt` 中的非空行数。计算公式依据 Tajima (1989)，显著性按照 Table 2 的临界值判定，查表方式与 SNPGenie 一致。

| 脚本 | 角色 | 必需参数 |
| --- | --- | --- |
| `selection_test_tajima.py` | 按基因的 CDS Tajima’s D | `--variant-data`、`--sample-dir`、`--cds-lengths` |
| `selection_test_tajima_igs.py` | 按 IGS 的 D，并给出合并的全基因组口径 D | `--variant-data`、`--sample-dir`、`--region-info`、`--cds-lengths`、`--group-info`、`--genome-regions` |
| `visualize_tajima_d.py` | 可选 CDS 出图 | `--tajima-results`、`--variant-data`、`--sample-dir`、`--species-order` |

`--sample-dir` 是含 `{物种}/sample.txt` 的目录（一般为 single-IR FASTA 树）。CDS 脚本的 `--gene-function-map` 默认使用本目录随仓发布的 `gene_function_map.csv`。

```bash
cd scripts/selection_test_results

python3 selection_test_tajima.py \
  --variant-data /path/to/all_combined_data.csv \
  --sample-dir /path/to/single_IR_fastas \
  --cds-lengths /path/to/output_cds_lengths.csv \
  --output-dir /path/to/tajima_d

python3 selection_test_tajima_igs.py \
  --variant-data /path/to/all_combined_data.csv \
  --sample-dir /path/to/single_IR_fastas \
  --region-info /path/to/output_gene_intergenic_intron_pos_length.csv \
  --cds-lengths /path/to/output_cds_lengths.csv \
  --group-info /path/to/group_info.csv \
  --genome-regions /path/to/species_genome_regions.csv \
  --output-dir /path/to/tajima_d_igs

python3 visualize_tajima_d.py \
  --tajima-results /path/to/all_species_tajima_d_summary.tsv \
  --variant-data /path/to/all_combined_data.csv \
  --sample-dir /path/to/single_IR_fastas \
  --species-order /path/to/species_label_order.csv \
  --output-dir /path/to/tajima_d/figures
```

可选参数：`--min-S`（默认 3）、`--min-n`（默认 4）；IGS 的 `--skip-plots`；作图脚本的 `--hotspot-genes`。分离位点数 *S* < 3 或样本量 *n* < 4 的基因 / 区段会跳过。

---

## 4. 结构筛查（v4.1）

**用途：** 通过全长序列比对寻找**区块间倒位或易位候选**。分析使用与变异检测一致的 **single-IR 线性** FASTA 坐标，并报告这些比对能够表示的结构模式。

**主链：** minimap2（`asm5`）→ 双轴共线性区块合并（`collinear_block_merge.py`）→ 区块间事件。常用区块参数：`min_len=500`、`max_gap=200`、`overlap_tol=50`、`translocation_jump=10000`。BLAST 可用于交叉核对。若 `PATH` 中没有 SyRI，`syri_representative_screen.py` 会改用 nucmer 和 `show-coords`。

**主入口：** `scripts/synteny_analysis/synteny_block_analysis.py`  
主入口为 `synteny_block_analysis.py`。其它模块会从 `window_rearrangement_test.py` 和 `audit_window_skips.py` 导入辅助函数。

### 必须在本地创建的内容

Python 脚本会通过 `yaml.safe_load` 读取 `--config`。用户在本地创建 YAML 文件。相对于本克隆的默认路径为：

| 参数 | 默认本地路径 |
| --- | --- |
| `--config` | `config/synteny_screening.yml` |
| `--sample-dir` | `data/00-rawfa/single_IR` |
| `--ref-dir` | `data/ref/single_IR` |

把 FASTA 与 GenBank 放到上述目录只是图方便，仍需一份能解析成字典的 YAML。`synteny_block_analysis.py` 在缺少 `block` / `minimap2` / `blastn` 键时会填入代码默认值；`aggregate_synteny_results.py` 目前会读取 `results_dirs.win200` 做归档窗检对照，因此传入的 YAML 必须包含该键（可以指向本地空目录）。

```bash
cd scripts/synteny_analysis

python3 synteny_block_analysis.py \
  --config /path/to/synteny_screening.yml \
  --aligner minimap2 --save-paf-all \
  --sample-dir /path/to/00-rawfa/single_IR \
  --ref-dir /path/to/ref/single_IR \
  --output-dir /path/to/blocks

python3 aggregate_synteny_results.py \
  --config /path/to/synteny_screening.yml \
  --blocks-dir /path/to/blocks \
  --output-dir /path/to/synteny_out

python3 synteny_report_orchestrator.py \
  --config /path/to/synteny_screening.yml \
  --output-dir /path/to/synteny_out

python3 test_collinear_block_merge.py
```

**接受** `--config`（默认 `config/synteny_screening.yml`）的脚本：`synteny_block_analysis.py`、`aggregate_synteny_results.py`、`synteny_report_orchestrator.py`、`synteny_block_sensitivity.py`、`run_phase1b_viz.py`、`paf_local_anomaly.py`、`syri_representative_screen.py` 以及作图脚本。`synteny_block_sensitivity.py` 还必须加 `--from-cache`，并提供 PAF 缓存。

**不接受** `--config` 的脚本：`build_method_agreement.py`（`--summary`、`--rep-list`、`--blast-blocks-dir`、`--out`）、`audit_synteny_results.py`（仅 `--json`）。二者默认读取本克隆 `results/synteny_analysis/` 下的结果，除非改路径参数。

`test_collinear_block_merge.py` 用于测试区块合并逻辑。

---

## 本地保存的输入与输出

- Hotspot、Tajima 和结构筛查的**结果表**，包括 PAF 缓存  
- YAML **配置文件**  
- 质体 FASTA / GenBank 面板与原始变异总表  

公开仓包含上文列出的四类脚本。

---

## 引用

使用本脚本请引用 pan-plastome 稿件；使用 R 包请引用 [cpopvar](https://github.com/Liyun-Nie/cpopvar)。公式与设计说明：[`docs/METHODS.zh-CN.md`](docs/METHODS.zh-CN.md) / [`docs/METHODS.md`](docs/METHODS.md)。
