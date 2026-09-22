# Automatic tests for hotspot threshold sensitivity helpers.
# Run: Rscript scripts/hotspot_threshold_sensitivity/test_hotspot_threshold_sensitivity.R

this_file <- tryCatch({
  ofile <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(sub("^--file=", "", ofile[[1]]))
}, error = function(e) normalizePath("test_hotspot_threshold_sensitivity.R"))

mod_dir <- dirname(this_file)
scripts_root <- dirname(mod_dir)
pan_dir <- dirname(scripts_root)
workspace_root <- normalizePath(pan_dir)
scripts_dir <- mod_dir
source(file.path(mod_dir, "hotspot_sensitivity_helpers.R"), local = FALSE)

n_pass <- 0L
n_fail <- 0L
failures <- character()

check <- function(desc, cond) {
  ok <- isTRUE(cond)
  if (ok) {
    n_pass <<- n_pass + 1L
    cat("PASS  ", desc, "\n", sep = "")
  } else {
    n_fail <<- n_fail + 1L
    failures <<- c(failures, desc)
    cat("FAIL  ", desc, "\n", sep = "")
  }
  invisible(ok)
}

almost_equal <- function(a, b, tol = 1e-10) {
  is.finite(a) && is.finite(b) && abs(a - b) <= tol
}

cat("=== hotspot threshold sensitivity tests ===\n")

# 1) quantile mapping
check("top 10% -> 0.90", identical(quantile_prob_for_top_percent(10), 0.90))
check("top 15% -> 0.85", identical(quantile_prob_for_top_percent(15), 0.85))
check("top 20% -> 0.80", identical(quantile_prob_for_top_percent(20), 0.80))
check("top 25% -> 0.75", identical(quantile_prob_for_top_percent(25), 0.75))
check("top 30% -> 0.70", identical(quantile_prob_for_top_percent(30), 0.70))
err <- try(quantile_prob_for_top_percent(12), silent = TRUE)
check("unsupported cutoff errors", inherits(err, "try-error"))

# 2) nonzero universe: zeros must not change the threshold
df_nz <- data.frame(
  species = rep("A", 6),
  locus_id = letters[1:6],
  locus_display = letters[1:6],
  region_type = "CDS",
  frequency_per_kb = c(0, 0, 1, 2, 3, 4),
  stringsAsFactors = FALSE
)
got <- identify_hotspots_nonzero_quantile(df_nz, 0.75, 7L)
thr_nonzero_only <- as.numeric(stats::quantile(c(1, 2, 3, 4), 0.75, type = 7, names = FALSE))
thr_with_zeros <- as.numeric(stats::quantile(c(0, 0, 1, 2, 3, 4), 0.75, type = 7, names = FALSE))
check("nonzero threshold uses only freq>0", almost_equal(unique(got$threshold), thr_nonzero_only))
check("nonzero threshold differs from all-row quantile", !almost_equal(thr_nonzero_only, thr_with_zeros))
check("zeros are never selected when threshold>0", !any(got$is_hotspot[df_nz$frequency_per_kb == 0]))

# 3) ties retained by >=
df_ties <- data.frame(
  species = rep("A", 5),
  locus_id = letters[1:5],
  locus_display = letters[1:5],
  region_type = "CDS",
  frequency_per_kb = c(1, 2, 3, 3, 3),
  stringsAsFactors = FALSE
)
got_ties <- identify_hotspots_nonzero_quantile(df_ties, 0.75, 7L)
thr_ties <- as.numeric(stats::quantile(c(1, 2, 3, 3, 3), 0.75, type = 7, names = FALSE))
check("ties threshold is 3", almost_equal(unique(got_ties$threshold), 3))
check("ties: all three 3s are selected", sum(got_ties$is_hotspot) == 3L)
check("ties inflate actual fraction above nominal 0.25",
      unique(got_ties$actual_selected_fraction) > unique(got_ties$nominal_top_fraction))
check("R type 7 used for ties fixture", almost_equal(thr_ties, 3))

# 4) BH family is region_type x cutoff, not pooled
# Use unequal families so pooled adjustment is observably different.
p <- c(0.001, 0.04, 0.02, 0.80)
family <- c("CDS|25", "CDS|25", "IGS|25", "IGS|25")
adj_by <- bh_adjust_by_family(p, family)
adj_pool <- stats::p.adjust(p, method = "BH")
check("BH within family differs from pooled BH", !isTRUE(all.equal(adj_by, adj_pool)))
check("BH CDS family first p stays smallest", adj_by[1] < adj_by[2])
check("BH families of size 2: 0.001 -> 0.002", almost_equal(adj_by[1], 0.002))
check("BH families of size 2: 0.04 -> 0.04", almost_equal(adj_by[2], 0.04))
check("IGS family is adjusted independently",
      almost_equal(adj_by[3], 0.04) && almost_equal(adj_by[4], 0.80))

# 5) toy hypergeometric, hand-checked
# N=20, K=4, n=8, x=3; P(X>=3)=P(3)+P(4)
comb <- function(a, b) {
  if (b < 0 || b > a) return(0)
  choose(a, b)
}
N <- 20; K <- 4; n <- 8; x <- 3
p_hand <- (comb(K, 3) * comb(N - K, n - 3) + comb(K, 4) * comb(N - K, n - 4)) / comb(N, n)
p_fun <- hypergeometric_upper_tail(x, K, N, n)
p_r <- stats::phyper(x - 1, K, N - K, n, lower.tail = FALSE)
check("toy hypergeom matches choose() expansion", almost_equal(p_fun, p_hand, 1e-12))
check("toy hypergeom matches phyper upper tail", almost_equal(p_fun, p_r, 1e-12))
check("x=0 yields p=1", almost_equal(hypergeometric_upper_tail(0, K, N, n), 1))

# 6) pooled hypergeometric sensitivity must use only nonzero pairs
toy_nz_hyper <- data.frame(
  species = c("S1", "S1", "S1", "S2", "S2", "S2"),
  locus_id = c("L1", "L2", "L3", "L1", "L2", "L3"),
  locus_display = c("L1", "L2", "L3", "L1", "L2", "L3"),
  region_type = "CDS",
  frequency_per_kb = c(4, 2, 0, 3, 5, 0),
  is_hotspot = c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)
toy_nz_tests <- locus_test_table(toy_nz_hyper)
check("nonzero hypergeometric N excludes zero pairs",
      all(toy_nz_tests$N_nonzero == 4L))
check("nonzero hypergeometric K is locus positive coverage",
      all(toy_nz_tests$K_nonzero[toy_nz_tests$locus_id %in% c("L1", "L2")] == 2L))
check("zero-positive locus stays in prespecified BH family with K=0 and P=1",
      with(toy_nz_tests[toy_nz_tests$locus_id == "L3", ],
           K_nonzero == 0L && p_hyper_nonzero_raw == 1))
check("nonzero hypergeometric expected x uses nonzero universe",
      all(vapply(toy_nz_tests$expected_x_nonzero[toy_nz_tests$K_nonzero > 0],
                 almost_equal,
                 logical(1), b = 1)))
check("pooled hypergeometric is explicitly supplementary",
      all(toy_nz_tests$hypergeom_role == "coarse_sensitivity_only"))

# 7) fixed-seed stratified permutation reproducibility
toy <- data.frame(
  species = c(rep("S1", 4), rep("S2", 4)),
  locus_id = rep(c("L1", "L2", "L3", "L4"), 2),
  locus_display = rep(c("L1", "L2", "L3", "L4"), 2),
  region_type = "CDS",
  frequency_per_kb = c(4, 3, 2, 1, 5, 1, 1, 1),
  stringsAsFactors = FALSE
)
toy$is_hotspot <- c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
perm_a <- stratified_permutation_pvalues(toy, n_perm = 200L, seed = 20251104L)
perm_b <- stratified_permutation_pvalues(toy, n_perm = 200L, seed = 20251104L)
perm_c <- stratified_permutation_pvalues(toy, n_perm = 200L, seed = 99L)
check("permutation same seed is identical", isTRUE(all.equal(perm_a$p_perm_raw, perm_b$p_perm_raw)))
check("permutation different seed can differ", !isTRUE(all.equal(perm_a$p_perm_raw, perm_c$p_perm_raw)))
check("permutation p in (0,1]", all(perm_a$p_perm_raw > 0 & perm_a$p_perm_raw <= 1))
check("widely shared locus has smaller perm p than a private one",
      perm_a$p_perm_raw[perm_a$locus_id == "L1"] < perm_a$p_perm_raw[perm_a$locus_id == "L4"])
perm_full <- stratified_permutation_pvalues(
  toy,
  n_perm = 200L,
  seed = 20251104L,
  return_null_distribution = TRUE
)
check("full permutation return contains p-values and null distribution",
      identical(sort(names(perm_full)), c("null_distribution", "pvalues")))
check("full permutation p-values match compact return",
      isTRUE(all.equal(perm_full$pvalues, perm_a)))
null_sums <- tapply(
  perm_full$null_distribution$n_permutations,
  perm_full$null_distribution$locus_id,
  sum
)
null_probs <- tapply(
  perm_full$null_distribution$probability,
  perm_full$null_distribution$locus_id,
  sum
)
check("each locus null distribution contains every permutation",
      all(null_sums == 200L))
check("each locus null probabilities sum to one",
      all(abs(null_probs - 1) < 1e-12))
check("null recurrence counts stay within species range",
      all(perm_full$null_distribution$x_null >= 0L &
            perm_full$null_distribution$x_null <= 2L))

# 8) 25% baseline regression against frozen session + Table 12 root cause
session_m03 <- file.path(workspace_root, "app_data", "sessions", "workflow_test_2025_11_4",
                         "results", "plots", "M03_hotspot", "M03_hotspot_standard",
                         "M03_all_genes_ranked_with_thresholds.csv")
session_m04 <- file.path(workspace_root, "app_data", "sessions", "workflow_test_2025_11_4",
                         "results", "plots", "M04_poigs_hotspot_engine", "M04_igs_hotspot_standard",
                         "M04_all_poigs_with_thresholds.csv")
table12 <- file.path(workspace_root, "manuscript", "pan-plastome", "version",
                     "first_revision_0525",
                     "Supplementary Table 12. Variant frequency for each IGS based on  SNV + indel + CPX.tsv")

if (file.exists(session_m03) && file.exists(session_m04)) {
  cds <- standardize_frequency_table(read_csv_utf8(session_m03), "CDS", "gene")
  igs <- standardize_frequency_table(read_csv_utf8(session_m04), "IGS", "poiGS_ID")
  cds25 <- identify_hotspots_nonzero_quantile(cds, 0.75, 7L)
  igs25 <- identify_hotspots_nonzero_quantile(igs, 0.75, 7L)
  cds_sh <- summarize_sharing_classes(sharing_by_locus(cds25)$n_species)
  igs_sh <- summarize_sharing_classes(sharing_by_locus(igs25)$n_species)
  check("session CDS 25% calls=109", sum(cds25$is_hotspot) == 109L)
  check("session CDS 25% unique=50", cds_sh$n_unique == 50L)
  check("session CDS private/shared2/shared3+=24/6/20",
        cds_sh$private == 24L && cds_sh$shared2 == 6L && cds_sh$shared3plus == 20L)
  check("session CDS stored flags match recompute",
        identical(as.logical(cds25$is_hotspot), as.logical(cds$stored_is_hotspot)))
  check("session IGS 25% calls=183", sum(igs25$is_hotspot) == 183L)
  check("session IGS 25% unique=63", igs_sh$n_unique == 63L)
  check("session IGS private/shared2/shared3+=21/10/32",
        igs_sh$private == 21L && igs_sh$shared2 == 10L && igs_sh$shared3plus == 32L)
  check("session IGS stored flags match recompute",
        identical(as.logical(igs25$is_hotspot), as.logical(igs$stored_is_hotspot)))
  gb <- igs25[igs25$species == "Gossypium_barbadense", ]
  check("GB rpl16-rps3 is below session type7 threshold",
        isTRUE(gb$frequency_per_kb[gb$locus_display == "rpl16-rps3"] < unique(gb$threshold)[[1]]))
  check("GB rpl16-rps3 is not a session hotspot",
        isTRUE(!gb$is_hotspot[gb$locus_display == "rpl16-rps3"]))
  check("approved manuscript baseline matches session IGS",
        sum(igs25$is_hotspot) == MANUSCRIPT_BASELINE_25$IGS$n_calls)
} else {
  cat("SKIP optional 25% session regression (ranked tables not in this clone)\n")
}

if (file.exists(table12)) {
  t12 <- read_tsv_utf8(table12)
  flag_col <- grep("hotspot", names(t12), ignore.case = TRUE, value = TRUE)[[1]]
  n_true <- sum(as_logical_flag(t12[[flag_col]]))
  t12_std <- data.frame(
    species = species_key(t12$Species),
    locus_id = as.character(t12$IGS),
    locus_display = as.character(t12$IGS),
    region_type = "IGS",
    frequency_per_kb = as.numeric(t12[[grep("frequency", names(t12), ignore.case = TRUE, value = TRUE)[[1]]]]),
    stringsAsFactors = FALSE
  )
  t12_re <- identify_hotspots_nonzero_quantile(t12_std, 0.75, 7L)
  extra <- t12_std$species == "Gossypium_barbadense" & t12_std$locus_display == "rpl16-rps3"
  check("Table12 stored flags = 184", n_true == 184L)
  check("Table12 type7 recompute = 183", sum(t12_re$is_hotspot) == 183L)
  check("Table12 extra TRUE is GB rpl16-rps3",
        isTRUE(as_logical_flag(t12[[flag_col]][extra])) &&
          isTRUE(!t12_re$is_hotspot[extra]))
} else {
  cat("SKIP optional Table12 183/184 root-cause test (file not in this clone)\n")
}

# GB 44 vs 43 nonzero fixture: both type7 thresholds exclude 5.988; snapping to 5.988 includes it
gb_nz44 <- c(
  17.2413793103448, 14.0845070422535, 9.92063492063492, 9.75609756097561,
  8.48765432098766, 7.8740157480315, 7.8125, 7.75193798449612, 7.40740740740741,
  6.75349734684033, 6.5359477124183, 5.98802395209581, 5.91016548463357,
  5.81395348837209, 5.6980056980057, 5.68181818181818, 5.55555555555556,
  4.98504486540379, 4.92610837438424, 4.83091787439614, 4.55927051671733,
  4.49640287769784, 4.4762757385855, 4.29876410531972, 4.27350427350427,
  4.04448938321537, 3.9032006245121, 3.80228136882129, 3.37268128161889,
  3.25379609544469, 3.19148936170213, 3.18725099601594, 2.63852242744063,
  2.57731958762887, 2.53646163601776, 2.3696682464455, 2.27963525835866,
  2.21238938053097, 1.91623323295921, 1.77304964539007, 1.68406871000337,
  1.42450142450142, 0.929368029739777, 0.824780058651026
)
thr44 <- as.numeric(stats::quantile(gb_nz44, 0.75, type = 7, names = FALSE))
thr43 <- as.numeric(stats::quantile(gb_nz44[gb_nz44 != 0.929368029739777], 0.75, type = 7, names = FALSE))
boundary <- 5.98802395209581
check("GB n=44 type7 threshold ~6.125", almost_equal(thr44, 6.12500489217643, 1e-8))
check("GB n=43 type7 threshold ~6.262", almost_equal(thr43, 6.261985832, 1e-8))
check("both type7 thresholds exclude rpl16-rps3", boundary < thr44 && boundary < thr43)
check("snapping threshold to 5.988 would include the extra call", boundary >= boundary)

cat("\n", n_pass, " passed, ", n_fail, " failed\n", sep = "")
if (n_fail > 0) {
  cat("Failures:\n- ", paste(failures, collapse = "\n- "), "\n", sep = "")
  quit(status = 1)
}
quit(status = 0)
