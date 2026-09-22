# Technical methods: four pan-plastome analyses

This document explains the **design, workflow, and references** for the four groups of scripts. Commands and input requirements are in the root [README.md](../README.md). Chinese version: [METHODS.zh-CN.md](METHODS.zh-CN.md).

The descriptions below follow the implementation in this repository.

---

## 0. Overview

| Module | Scientific question | Core inference |
| --- | --- | --- |
| Threshold Hotspot | Do high-frequency tails of *nonzero* loci recur as Hotspots across species? | Species-stratified permutation + BH |
| Threshold-free Hotspot | Is there shared high-rank structure in the full within-species rankings? | Species-wise label permutation; `Q_high` / `C_global` / `M_max`; maxT |
| Tajima’s D | Do CDS / IGS site-frequency spectra depart from neutrality? | Tajima (1989), haploid *n* |
| Structural screen v4.1 | Macro rearrangements that short-read SNV/indel pipelines may miss? | Alignment-first collinear-block screen |

**Shared principles**

1. Both Hotspot analyses create the randomised null distribution **within each species** and then compare the same loci across species.  
2. The threshold and threshold-free analyses answer different questions and should be interpreted separately.  
3. Structural screening uses **single-IR linear** coordinates and covers rearrangements represented in those coordinates.  
4. The full Hotspot analyses use **99,999** permutation draws. The repository stores the scripts, while input genomes and generated permutation data remain local.

---

## 1. Threshold-based Hotspot (`hotspot_threshold_sensitivity`)

### 1.1 Design

A Hotspot is the top *q* tail (*q* ∈ {10%, 15%, 20%, 25%, 30%}) of **positive** within-species frequencies (`frequency_per_kb > 0`).

- CDS: SNV only. IGS: SNV + indel + CPX. All rates are per kb.  
- Zeros are excluded because they mean that no variation was detected at that locus in that species. Including many zeros could reduce the cutoff to zero and incorrectly label zero-frequency loci as Hotspots. This choice is specific to this analysis.  
- After loci are labelled as Hotspots or non-Hotspots, the cross-species test measures how often the **same locus is repeatedly labelled as a Hotspot**.

### 1.2 Workflow

1. Read session ranked CDS/IGS tables under `--workspace` (see README).  
2. For each species and each *q*, compute the threshold on the positive subset with R `quantile(..., type = 7)` (Hyndman–Fan type 7). Mark `frequency ≥ threshold`. Ties at the cutoff are kept, so the realised fraction can exceed the nominal Top *q*.  
3. Observed recurrence \(X_l=\sum_s H_{sl}\).  
4. **Species-stratified permutation** (primary): keep each species’ eligible set \(E_s\) and Hotspot count \(h_s\); redraw labels without replacement inside \(E_s\); no mixing across species; *B* replicates; upper-tail  
   \(P_l=(1+b_l)/(1+B)\) (Phipson–Smyth form). Official *B* = 99,999.  
5. Benjamini–Hochberg within each region-type × cutoff family.  
Entry: `scripts/hotspot_threshold_sensitivity/hotspot_threshold_sensitivity.R`.

### 1.3 References

1. Hyndman, R.J. & Fan, Y. Sample quantiles in statistical packages. *Am. Stat.* **50**, 361–365 (1996).  
2. Strasser, H. & Weber, C. On the asymptotic theory of permutation statistics. *Math. Methods Stat.* **8**, 220–250 (1999).  
3. Phipson, B. & Smyth, G.K. Permutation P-values should never be zero. *Stat. Appl. Genet. Mol. Biol.* **9**, Article 39 (2010).  
4. Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate. *J. R. Stat. Soc. B* **57**, 289–300 (1995).  
---

## 2. Threshold-free Hotspot (`hotspot_threshold_free_analysis`)

### 2.1 Design

This analysis uses the **full** within-species ranking of analysable loci, including true zeros. Missing coverage remains missing.

1. Rank frequencies ascending within species (average ranks for ties, including zero ties).  
2. Rankit percentile \(u=(r-0.5)/n\), then \(z=\Phi^{-1}(u)\).  
3. Locus score \(T_l=\sum_s z_{sl}/\sqrt{K_l}\) (*K_l* = species coverage). Division by \(\sqrt{K_l}\) standardises the scale. P-values come from the within-species permutations.  
4. The three global statistics describe complementary patterns:  
   - \(Q_{\mathrm{high}}=\sum_l[\max(T_l,0)]^2\) — primary omnibus for joint high-rank structure;  
   - \(C_{\mathrm{global}}\) — overlap-weighted mean Pearson correlation of species-internal percentiles \(u\);  
   - \(M_{\mathrm{max}}=\max_l T_l\) — whether at least one locus is extreme.  
   All three use the same within-species ranks.  
5. Null: permute locus labels **within each species** (preserve the frequency multiset, zero ties, and missingness).  
6. Locus-level inference: Westfall–Young single-step maxT against the permutation max-\(T\) null. BH on marginal locus P-values is an auxiliary contrast.

### 2.2 Workflow

1. Same ranked tables as the threshold suite. SHA-256 and row-count gates are read from YAML (SHA can fall back to helper constants).  
2. Ranks → \(u\) → \(z\); compute \(T_l\), \(Q_{\mathrm{high}}\), \(C_{\mathrm{global}}\), \(M_{\mathrm{max}}\).  
3. Primary 99,999 permutations. Robustness (fewer draws; zero prevalence, positive intensity, length residual, centred ranks, genus down-sampling) runs only if YAML sets `run_robustness: true`.  
4. Write global / locus / contribution / pairwise tables and figures.

Entry: `scripts/hotspot_threshold_free_analysis/hotspot_threshold_free_analysis.R`.

### 2.3 References

1. Bliss, C.I. *Statistics in Biology* Vol. 1 (McGraw-Hill, 1967). (Rankit)  
2. Beasley, T.M., Erickson, S. & Allison, D.B. Rank-based inverse normal transformations… *Behav. Genet.* **39**, 580–595 (2009).  
3. Westfall, P.H. & Young, S.S. *Resampling-based Multiple Testing* (Wiley, 1993).  
4. Phipson & Smyth (2010); Strasser & Weber (1999).  
5. Benjamini & Hochberg (1995) — auxiliary BH; primary locus control is maxT.  
6. Spearman, C. The proof and measurement of association between two things. *Am. J. Psychol.* **15**, 72–101 (1904). — rank association underpinning \(C_{\mathrm{global}}\), which uses species-internal percentiles rather than a re-ranked overlapping subset.

The names `Q_high`, `C_global`, and `M_max` are defined for this analysis.

---

## 3. Tajima’s D (`selection_test_results`)

### 3.1 Design

Compare pairwise diversity \(\pi\) with Watterson’s \(\theta_W=S/a_1\). Under neutrality and constant size the expectations match; \(D=(\pi-\theta_W)/\sqrt{\widehat{\mathrm{Var}}}\) can reflect selection **or** demography.

Haploid plastome implementation:

- *n* = number of accessions, counted as non-empty lines in `{Species}/sample.txt`;  
- *S* = number of unique SNP positions in the curated variant table;  
- \(\pi=\sum_p 2j(p)\,(n-j(p))/(n(n-1))\), where \(j(p)\) is the number of distinct samples carrying an alternate allele at site *p* (samples absent from the variant table are treated as reference at that site);  
- biological \(\theta\) is \(2N\mu\), but **D’s formula has no ploidy factor**—only *n* and *S*;  
- Tajima (1989) derives the variance under a no-recombination model. The implementation uses the Table 2 critical values for plastid data;  
- significance: SNPGenie-style lookup in Tajima (1989) Table 2. Limits are asymmetric. Sample sizes are matched by interval (e.g. *n* ≥ 50 uses the *n* = 50 row). Genes / regions with *S* < 3 or *n* < 4 are skipped.

Negative *D* indicates an excess of rare variants and may be associated with purifying selection or population expansion. Positive *D* indicates an excess of intermediate-frequency variants and may be associated with balancing selection or population contraction. These alternatives are interpreted together with biological and demographic evidence.

CDS: one *D* per species × gene. IGS: per-region *D*, plus pooled genome-wide *D* for IGS, CDS, and all SNPs (the last scaled by single-IR length from `--genome-regions`). Those pooled tables are a different aggregation from the per-gene CDS table.

### 3.2 Workflow

1. Read the variant table, length / region tables, and sample directory.  
2. Aggregate segregating sites; compute *S*, \(\pi\), *D*, and significance bands.  
3. Write summaries; optional `visualize_tajima_d.py` for CDS figures.

Entries: `selection_test_tajima.py`, `selection_test_tajima_igs.py`. Required flags differ (README).

### 3.3 References

1. Tajima, F. Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics* **123**, 585–595 (1989).  
2. SNPGenie `Tajima_D.R` (Chase W. Nelson) — implementation cross-check for the Table 2 lookup and coefficient chain.

---

## 4. Structural screen v4.1 (`synteny_analysis`)

### 4.1 Design

Short-read SNV/indel pipelines may miss large inversions and translocations. This suite therefore starts with whole-assembly alignment and then screens for changes in collinearity:

1. Align each sample to the species reference. Primary: minimap2 `asm5`. Cross-check: blastn on representatives.  
2. Dual-axis collinear block merge (`collinear_block_merge.py`): merge adjacent same-strand segments when both query and reference gaps lie in \([-\texttt{overlap_tol}, \texttt{max_gap}]\). Detect **inter-block** inversion / translocation under primary parameters.  
3. Sensitivity grids from cached PAF, method-agreement tables, PAF visualisation, within-chain local anomalies, and representative nucmer checks.  
4. The primary evidence chain uses alignment blocks. `window_rearrangement_test.py` and `audit_window_skips.py` provide helper functions imported by other modules.

Coordinates: single-IR **linear** FASTA, matching the variant-calling references. Detection covers events represented in these linear coordinates.

**SyRI.** `syri_representative_screen.py` runs SyRI when `syri` is available on `PATH`. Otherwise it uses nucmer and `show-coords`. The **nucmer fallback** records reverse-strand blocks as inversion evidence and multiple non-collinear clusters as rearrangement candidates.

Typical primary block settings: `min_len=500`, `max_gap=200`, `overlap_tol=50`, `translocation_jump=10000` (exact values are set in the local YAML). Event summaries report the structural patterns detected in the single-IR linear alignments.

### 4.2 Workflow (conceptual)

```
FASTA + GenBank ref
  → synteny_block_analysis.py
  → aggregate_synteny_results.py
  → (optional) sensitivity / BLAST agreement / PAF figures / local anomaly
  → synteny_report_orchestrator.py
  → (optional) syri_representative_screen.py   # SyRI if installed, else nucmer
```

Commands and which scripts accept `--config` are listed in the README. Several follow-on scripts assume result files under the clone unless paths are overridden.

### 4.3 References

1. Li, H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* **34**, 3094–3100 (2018).  
2. Camacho, C. *et al.* BLAST+: architecture and applications. *BMC Bioinformatics* **10**, 421 (2009).  
3. Kurtz, S. *et al.* Versatile and open software for comparing large genomes. *Genome Biol.* **5**, R12 (2004). (MUMmer / nucmer)  
4. Goel, M. *et al.* SyRI: finding genomic rearrangements and local sequence differences from whole-genome assemblies. *Genome Biol.* **20**, 277 (2019). — optional validation tool.

Merge rules and event definitions are those implemented in `collinear_block_merge.py` and `synteny_block_analysis.py`.

---

## 5. Relationship to cpopvar

| Object | Where | Note |
| --- | --- | --- |
| Downstream Hotspot analyses (stratified permutation / ranks / maxT) | **This repository** | Implemented by the scripts described above |
| Ranked-table generation | cpopvar session outputs | Used as input by the scripts in this repository |

---

## 6. Reproducibility notes

1. Genomes, PAF caches, result TSVs, and YAML files are stored locally.  
2. R Hotspot readers use flat YAML keys. Threshold-free analysis reads the `expected_*` validation keys. Synteny Python readers require a YAML file that `yaml.safe_load` parses as a mapping.  
3. Full Hotspot analyses use 99,999 draws; helper tests and shorter runs support code checks and development.  
4. Mini tests cover the Hotspot helpers and `test_collinear_block_merge.py`. Input-validation checks use local ranked tables under a cpopvar-style tree.

---

## 7. Maintenance

- Usage: root `README.md` / `README.zh-CN.md`  
- Methods: this file / `METHODS.zh-CN.md`  
- Maintainer: Liyun Nie \<nieliyun18@163.com\> · ORCID 0000-0002-5288-0041
