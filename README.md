# pan-plastome

Scripts for three parts of the pan-plastome study: **hypervariable-region (Hotspot) analysis**, **Tajima’s D**, and **alignment-based structural screening**.

| | |
| --- | --- |
| Maintainer | Liyun Nie \<nieliyun18@163.com\> |
| ORCID | [0000-0002-5288-0041](https://orcid.org/0000-0002-5288-0041) |
| License | MIT (2026) |
| Related R package | [Liyun-Nie/cpopvar](https://github.com/Liyun-Nie/cpopvar) |

Use these scripts to understand the analysis methods or to repeat them with **your own** ranked tables and genome assemblies. The `cpopvar` R package handles variant processing, ranked-table generation, and the interactive workflow. This repository contains the downstream Hotspot analyses: species-stratified permutation and threshold-free rank analysis with Westfall–Young maxT.

This repository distributes the analysis scripts. Users supply the result tables used as inputs, YAML configuration files, plastome FASTA/GenBank files, and sequence datasets locally. See the manuscript data-availability statement and [`docs/METHODS.md`](docs/METHODS.md).

**中文说明**：[README.zh-CN.md](README.zh-CN.md) · [docs/METHODS.zh-CN.md](docs/METHODS.zh-CN.md)

---

## Layout

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

| Directory | What it does |
| --- | --- |
| `hotspot_threshold_sensitivity/` | Threshold Hotspot (Top 10–30% of *positive* within-species frequencies); species-stratified permutation |
| `hotspot_threshold_free_analysis/` | Threshold-free ranks; `Q_high` / `C_global` / `M_max`; Westfall–Young maxT |
| `selection_test_results/` | CDS and IGS Tajima’s D (optional figures) |
| `synteny_analysis/` | Alignment-first collinear-block screen (v4.1) |

---

## Requirements

| Suite | Runtime | Typical extras |
| --- | --- | --- |
| Hotspot (both) | R ≥ 4.1 | Helper tests use base R |
| Tajima | Python 3 | `pandas`, `numpy`; figures need `matplotlib`, `seaborn`, `scipy` |
| Synteny | Python 3 | `pandas`, `PyYAML`, `biopython`; aligner `minimap2`; optional `blastn`, `nucmer` / `show-coords`, SyRI |

The complete analyses use the original input data and local configuration files. The small tests below check individual functions; the full analyses use 99,999 permutations.

---

## 1. Threshold-based Hotspot

**Purpose:** Within each species, the script marks the top 10%, 15%, 20%, 25%, and 30% of *nonzero* variant frequencies as Hotspots. CDS frequencies include SNVs only; IGS frequencies include SNVs, indels, and CPX variants and are reported per kb. The test then asks whether the same loci are marked as Hotspots in multiple species more often than expected by chance. The randomisation preserves each species’ eligible loci and number of Hotspots.

**Entry:** `scripts/hotspot_threshold_sensitivity/hotspot_threshold_sensitivity.R`  
**Helpers / tests:** `hotspot_sensitivity_helpers.R`, `test_hotspot_threshold_sensitivity.R`

### Inputs you must supply

Set `--workspace` to the **cpopvar-style tree root**. The script joins:

```
<workspace>/app_data/sessions/<session_id>/results/plots/M03_hotspot/M03_hotspot_standard/M03_all_genes_ranked_with_thresholds.csv
<workspace>/app_data/sessions/<session_id>/results/plots/M04_poigs_hotspot_engine/M04_igs_hotspot_standard/M04_all_poigs_with_thresholds.csv
```

The script adds `app_data/sessions/<session_id>` to `--workspace`. The candidate Hotspot CSV files in the same two result folders are also required. Set `session_id` in your YAML file to match your session directory.

**Config:** Pass a flat YAML with `--config`. Nested blocks are ignored. If the file is omitted, the script looks for `config/hotspot_threshold_sensitivity.yml` next to this clone and treats a missing file as empty (built-in defaults apply). YAML is still the practical way to record `session_id` and permutation settings. Example keys:

```yaml
session_id: your_session_id
n_perm: 99999
seed: 20251104
alpha: 0.05
quantile_type: 7
```

CLI overrides (hyphenated form is the documented spelling; `--n_perm` is also accepted): `--workspace`, `--output`, `--config`, `--n-perm`, `--seed`, `--skip-perm`.

### Run

```bash
# From the clone root
Rscript scripts/hotspot_threshold_sensitivity/hotspot_threshold_sensitivity.R \
  --config /path/to/hotspot_threshold_sensitivity.yml \
  --workspace /path/to/cpopvar_tree_root \
  --output /path/to/output_dir

# Helper tests (synthetic checks; optional 25% regression only if ranked tables sit under this clone)
Rscript scripts/hotspot_threshold_sensitivity/test_hotspot_threshold_sensitivity.R
```

`--skip-perm` runs the preprocessing and summary steps without the permutation stage.

---

## 2. Threshold-free Hotspot

**Purpose:** This analysis uses the full within-species ranking instead of choosing a Top-% cutoff. True zero values are kept, while missing values remain missing. The ranks are converted to Rankit inverse-normal scores, and labels are permuted within each species to test whether high-ranking loci are shared across species. `Q_high` is the main overall statistic; `C_global` and `M_max` describe complementary patterns. Locus-level testing uses Westfall–Young single-step maxT (`p_maxT`). Only loci covered in at least two species are included.

**Entry:** `scripts/hotspot_threshold_free_analysis/hotspot_threshold_free_analysis.R`  
**Helpers / tests:** `hotspot_threshold_free_helpers.R` (sources the threshold helpers), `test_hotspot_threshold_free_analysis.R`

`--workspace` has the same meaning as in §1. Ranked-table paths are identical.

### Configuration for a full run

Provide a flat YAML file with `--config`. The full run reads dataset-specific validation keys from this file.

The configuration should define the run settings and the validation values calculated from your own ranked tables:

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

Also provide the `expected_*` validation keys used by the script: row counts, nonzero counts, SHA-256 values, and Top 25% summary counts for CDS and IGS. Calculate these dataset-specific values from your ranked tables.

**CLI overrides** (hyphenated form documented; underscore aliases are accepted): `--workspace`, `--output`, `--config`, `--n-perm`, `--n-perm-robust`, `--seed`, `--mode`, `--chunk-size`, `--skip-perm`, `--skip-robustness`. Robustness permutations run when YAML has `run_robustness: true`; `--skip-robustness` disables this stage. `--mode dev` uses `n_perm_dev`, while `--n-perm` sets an explicit value.

```bash
Rscript scripts/hotspot_threshold_free_analysis/hotspot_threshold_free_analysis.R \
  --config /path/to/hotspot_threshold_free.yml \
  --workspace /path/to/cpopvar_tree_root \
  --output /path/to/output_dir

Rscript scripts/hotspot_threshold_free_analysis/test_hotspot_threshold_free_analysis.R
```

Helper tests cover the rank and permutation algebra. SHA / 25% regression checks run when the ranked tables are available under the clone. Full analysis runs use 99,999 draws.

---

## 3. Tajima’s D (CDS / IGS)

**Purpose:** Test whether the site-frequency spectra of CDS and IGS regions differ from neutral expectations. Plastid genomes are treated as haploid, so **n is the number of accessions**, counted as the non-empty lines in each species’ `sample.txt`. The calculation follows Tajima (1989), and significance is assigned from the Table 2 critical values using an SNPGenie-style lookup.

| Script | Role | Required arguments |
| --- | --- | --- |
| `selection_test_tajima.py` | Per-gene CDS Tajima’s D | `--variant-data`, `--sample-dir`, `--cds-lengths` |
| `selection_test_tajima_igs.py` | Per-IGS D plus pooled genome-wide D | `--variant-data`, `--sample-dir`, `--region-info`, `--cds-lengths`, `--group-info`, `--genome-regions` |
| `visualize_tajima_d.py` | Optional CDS figures | `--tajima-results`, `--variant-data`, `--sample-dir`, `--species-order` |

`--sample-dir` is a directory of `{Species}/sample.txt` folders (typically single-IR FASTA trees). `--gene-function-map` on the CDS script defaults to the shipped `gene_function_map.csv` in this folder.

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

Optional flags: `--min-S` (default 3), `--min-n` (default 4); IGS `--skip-plots`; visualiser `--hotspot-genes`. Genes with *S* < 3 or *n* < 4 are skipped.

---

## 4. Structural screening (v4.1)

**Purpose:** Use whole-assembly alignment to find **between-block** inversion or translocation candidates. The analysis uses **single-IR linear** FASTA coordinates, matching the variant-calling references, and reports structural patterns represented in these alignments.

**Primary chain:** minimap2 (`asm5`) → dual-axis collinear block merge (`collinear_block_merge.py`) → inter-block events. Typical block settings: `min_len=500`, `max_gap=200`, `overlap_tol=50`, `translocation_jump=10000`. BLAST is an optional cross-check. If SyRI is not on `PATH`, `syri_representative_screen.py` falls back to nucmer + `show-coords`.

**Main entry:** `scripts/synteny_analysis/synteny_block_analysis.py`  
The primary entry point is `synteny_block_analysis.py`. Other modules import utility functions from `window_rearrangement_test.py` and `audit_window_skips.py`.

### What you must create locally

Python scripts call `yaml.safe_load` on `--config`. There is **no shipped YAML**. Defaults (relative to this clone) are:

| Flag | Default local path |
| --- | --- |
| `--config` | `config/synteny_screening.yml` |
| `--sample-dir` | `data/00-rawfa/single_IR` |
| `--ref-dir` | `data/ref/single_IR` |

Placing FASTA and GenBank under those folders is convenient, but you still need a valid YAML dictionary. `synteny_block_analysis.py` fills `block` / `minimap2` / `blastn` if those keys are absent; `aggregate_synteny_results.py` currently reads `results_dirs.win200` for an archived window-test cross-check, so that key must exist in the YAML you pass (it may point at an empty local directory).

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

Scripts that **do** take `--config` (defaulting to `config/synteny_screening.yml`): `synteny_block_analysis.py`, `aggregate_synteny_results.py`, `synteny_report_orchestrator.py`, `synteny_block_sensitivity.py`, `run_phase1b_viz.py`, `paf_local_anomaly.py`, `syri_representative_screen.py`, and the plotting helpers. `synteny_block_sensitivity.py` additionally requires `--from-cache` and a PAF cache.

Scripts that use direct path options instead of `--config`: `build_method_agreement.py` (`--summary`, `--rep-list`, `--blast-blocks-dir`, `--out`) and `audit_synteny_results.py` (`--json` only). Their default result location is `results/synteny_analysis/` under the clone, and their path flags can select another location.

`test_collinear_block_merge.py` is a unit test for the block-merging logic.

---

## Inputs and outputs stored locally

- Hotspot, Tajima, and synteny **result tables**, including PAF caches  
- YAML **configuration files**  
- Plastome FASTA / GenBank panels and raw variant tables  

The public repository contains the four script groups listed above.

---

## Citation

Please cite the pan-plastome manuscript when using these scripts, and cite [cpopvar](https://github.com/Liyun-Nie/cpopvar) when using the R package. Formulas and design notes: [`docs/METHODS.md`](docs/METHODS.md) / [`docs/METHODS.zh-CN.md`](docs/METHODS.zh-CN.md).
