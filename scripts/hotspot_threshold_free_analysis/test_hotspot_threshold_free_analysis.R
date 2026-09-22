# Automatic tests for threshold-free hotspot concordance helpers.
# Run: Rscript scripts/hotspot_threshold_free_analysis/test_hotspot_threshold_free_analysis.R

this_file <- tryCatch({
  ofile <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(sub("^--file=", "", ofile[[1]]))
}, error = function(e) normalizePath("test_hotspot_threshold_free_analysis.R"))

mod_dir <- dirname(this_file)
scripts_root <- dirname(mod_dir)
pan_dir <- dirname(scripts_root)
workspace_root <- normalizePath(pan_dir)
scripts_dir <- mod_dir
source(file.path(scripts_root, "hotspot_threshold_sensitivity",
                 "hotspot_sensitivity_helpers.R"), local = FALSE)
source(file.path(mod_dir, "hotspot_threshold_free_helpers.R"), local = FALSE)

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

almost_equal_vec <- function(a, b, tol = 1e-10) {
  length(a) == length(b) && all(is.finite(a) & is.finite(b)) && max(abs(a - b)) <= tol
}

cat("=== hotspot threshold-free analysis tests ===\n")

helper_src <- readLines(file.path(scripts_dir, "hotspot_threshold_free_helpers.R"),
                        warn = FALSE, encoding = "UTF-8")
analysis_src <- readLines(file.path(scripts_dir, "hotspot_threshold_free_analysis.R"),
                          warn = FALSE, encoding = "UTF-8")
check("helpers do not call perform_sharing_statistical_test",
      !any(grepl("perform_sharing_statistical_test\\s*\\(", helper_src)))
check("analysis does not call perform_sharing_statistical_test",
      !any(grepl("perform_sharing_statistical_test\\s*\\(", analysis_src)))
check("automated reruns preserve the PM-reviewed detailed report",
      !any(grepl(
        "writeLines\\s*\\([^\\n]*RESULTS_REPORT_CN\\.md",
        analysis_src
      )) &&
        any(grepl("RESULTS_SUMMARY_AUTO_CN\\.md", analysis_src)))
check("helpers do not attach cpopvar",
      !any(grepl("library\\s*\\(\\s*cpopvar|source\\(.*cpopvar/R", helper_src)))
check("helpers do not use identify_hotspots",
      !any(grepl("identify_hotspots", helper_src)))
check("helpers do not use num_species >= 3",
      !any(grepl("num_species\\s*>=\\s*3", helper_src)))
check("cpopvar is not attached", isTRUE(assert_no_cpopvar_attached()))

# 1) zero ties: midrank of three zeros then 1,2
x_ties <- c(0, 0, 0, 1, 2)
sc_ties <- midrank_normal_scores(x_ties)
check("zero ties midranks are 2,2,2,4,5",
      almost_equal_vec(sc_ties$rank, c(2, 2, 2, 4, 5)))
check("zero ties n=5", identical(sc_ties$n, 5L))
u_hand <- (c(2, 2, 2, 4, 5) - 0.5) / 5
check("zero ties u=(rank-0.5)/n", almost_equal_vec(sc_ties$u, u_hand))
check("zero ties z=qnorm(u)", almost_equal_vec(sc_ties$z, stats::qnorm(u_hand)))

# 2) missing vs zero: NA is excluded from n, not ranked as zero
x_miss <- c(0, 0, NA, 1, 2)
sc_miss <- midrank_normal_scores(x_miss)
sc_zero_as_missing <- midrank_normal_scores(c(0, 0, 0, 1, 2))
check("missing is not ranked", is.na(sc_miss$rank[3]) && is.na(sc_miss$u[3]) && is.na(sc_miss$z[3]))
check("missing n excludes NA", identical(sc_miss$n, 4L))
check("two zeros with missing: midranks 1.5,1.5,NA,3,4",
      almost_equal_vec(sc_miss$rank[c(1, 2, 4, 5)], c(1.5, 1.5, 3, 4)))
check("treating missing as zero would change ranks",
      !almost_equal(sc_miss$rank[1], sc_zero_as_missing$rank[1]))

# 3) monotone transform invariance of midranks and z
set.seed(42)
x_raw <- c(0, 0.2, 1.5, 3, 8, 8)
sc_a <- midrank_normal_scores(x_raw)
sc_b <- midrank_normal_scores(exp(x_raw))
sc_c <- midrank_normal_scores(x_raw^3 + 4)
check("monotone exp preserves midrank", almost_equal_vec(sc_a$rank, sc_b$rank))
check("monotone cubic preserves z", almost_equal_vec(sc_a$z, sc_c$z))
check("monotone exp preserves u", almost_equal_vec(sc_a$u, sc_b$u))

# 4) Q/T hand calculation
z_hand <- matrix(c(1, 0, 1, 0), nrow = 2, byrow = TRUE,
                 dimnames = list(c("S1", "S2"), c("L1", "L2")))
th <- locus_T(z_hand, k_min = 2L)
check("hand T_L1 = sqrt(2)", almost_equal(th$T[["L1"]], sqrt(2)))
check("hand T_L2 = 0", almost_equal(th$T[["L2"]], 0))
check("hand Q_high = 2", almost_equal(Q_high(th$T), 2))
check("hand M_max = sqrt(2)", almost_equal(M_max(th$T), sqrt(2)))

# 4b) Q contribution decomposition and descriptive species counts
loc_hand <- data.frame(
  region_type = c("CDS", "CDS"),
  locus_canonical = c("L1", "L2"),
  locus_display = c("L1", "L2"),
  K_l = c(2L, 2L),
  included_in_cross_species = c(TRUE, TRUE),
  T_l = c(sqrt(2), 0),
  mean_percentile_u = c(0.75, 0.5),
  p_marginal = c(0.01, 1),
  p_bh = c(0.02, 1),
  p_maxT = c(0.03, 1),
  stringsAsFactors = FALSE
)
q_hand <- q_contribution_table(loc_hand, z_hand, alpha = 0.05)
check("Q contribution components sum to Q_high",
      almost_equal(sum(q_hand$Q_component), Q_high(loc_hand$T_l)))
check("Q contribution fractions sum to one",
      almost_equal(sum(q_hand$Q_contribution_fraction), 1))
check("L1 contributes 100 percent in hand example",
      almost_equal(q_hand$Q_contribution_percent[q_hand$locus_canonical == "L1"], 100))
check("L2 nonpositive contribution is zero",
      almost_equal(q_hand$Q_component[q_hand$locus_canonical == "L2"], 0))
check("Q table identifies M_max locus",
      identical(q_hand$locus_canonical[q_hand$is_M_max_locus], "L1"))
check("Q table counts above-median species",
      identical(q_hand$n_species_above_median[q_hand$locus_canonical == "L1"], 2L))
check("Q table keeps above-median species names",
      identical(q_hand$species_above_median[q_hand$locus_canonical == "L1"], "S1,S2"))
check("Q table carries maxT significance",
      isTRUE(q_hand$maxT_significant[q_hand$locus_canonical == "L1"]))

# 5) K=1 is audit-only
z_k <- matrix(c(1, 0, NA, 1, 0, 2), nrow = 2, byrow = TRUE,
              dimnames = list(c("S1", "S2"), c("A", "B", "C")))
tk <- locus_T(z_k, k_min = 2L)
check("K=1 locus C has K=1", identical(unname(tk$K[["C"]]), 1L))
check("K=1 locus T is NA", is.na(tk$T[["C"]]))
check("Q ignores K=1", almost_equal(Q_high(tk$T), Q_high(tk$T[c("A", "B")])))

# 6) permutation preserves species margins and missingness
z_margin <- matrix(c(1, 0, NA, 2, 0.5, 0, NA, 3, 1.2, NA),
                   nrow = 2, byrow = TRUE,
                   dimnames = list(c("S1", "S2"), paste0("L", 1:5)))
u_margin <- z_margin
perm_m <- permute_omnibus(z_margin, u_margin, n_perm = 25L, seed = 99L,
                          chunk_size = 10L, k_min = 2L, compute_C = FALSE)
ok1 <- is.finite(z_margin[1, ])
ok2 <- is.finite(z_margin[2, ])
set.seed(99L)
one <- z_margin
one[1, ok1] <- z_margin[1, ok1][sample.int(sum(ok1))]
check("missing cells stay missing after a manual perm",
      all(is.na(one[1, !ok1])) && all(is.na(z_margin[1, !ok1])))
check("manual perm preserves species-1 value multiset",
      almost_equal_vec(sort(one[1, ok1]), sort(z_margin[1, ok1])))
check("S2 missing pattern is species-specific",
      !identical(ok1, ok2))

# 7) fixed seed reproducibility
toy_z <- matrix(c(2, 1, 0, 0,
                  2.2, 0.4, 0.3, 0.1,
                  1.8, 0.5, 0.2, 0.0),
                nrow = 3, byrow = TRUE,
                dimnames = list(paste0("S", 1:3), paste0("L", 1:4)))
toy_u <- score_matrices_by_species(toy_z)$u
p1 <- permute_omnibus(toy_z, toy_u, n_perm = 80L, seed = 20260915L, chunk_size = 20L)
p2 <- permute_omnibus(toy_z, toy_u, n_perm = 80L, seed = 20260915L, chunk_size = 20L)
p3 <- permute_omnibus(toy_z, toy_u, n_perm = 80L, seed = 7L, chunk_size = 20L)
check("same seed same Q p", almost_equal(p1$p_Q, p2$p_Q))
check("same seed same M p", almost_equal(p1$p_M, p2$p_M))
check("same seed same C p", almost_equal(p1$p_C, p2$p_C))
check("same seed same pairwise permutation p",
      almost_equal_vec(p1$p_pair, p2$p_pair))
check("three toy species yield three pairwise tests",
      identical(length(p1$p_pair), 3L))
check("pairwise permutation p values are in (0,1]",
      all(p1$p_pair > 0 & p1$p_pair <= 1))
check("same seed same T p", almost_equal_vec(p1$p_T, p2$p_T))
check("different seed can change Q null draws",
      !isTRUE(all.equal(p1$Q_null, p3$Q_null)))
check("empirical p in (0,1]", all(p1$p_Q > 0 && p1$p_Q <= 1))
check("compute_C=FALSE skips pairwise inference",
      length(perm_m$p_pair) == 0L && length(perm_m$pair_keys) == 0L)

# 8) toy planted signal vs no-signal
sig <- matrix(0, 4, 8, dimnames = list(paste0("S", 1:4), paste0("L", 1:8)))
sig[, "L1"] <- c(10, 12, 11, 9)
set.seed(1)
sig[, -1] <- matrix(runif(4 * 7, 0, 1), 4, 7)
sig_sc <- score_matrices_by_species(sig)
sig_perm <- permute_omnibus(sig_sc$z, sig_sc$u, n_perm = 199L, seed = 123L, chunk_size = 50L)
check("planted shared high locus has the largest T",
      names(which.max(sig_perm$T_obs)) == "L1")
check("planted signal Q p is at the plus-one floor or very small",
      sig_perm$p_Q <= 5 / 200)
check("planted L1 maxT p is among the smallest",
      isTRUE(almost_equal(unname(sig_perm$p_maxT[["L1"]]), min(sig_perm$p_maxT, na.rm = TRUE))))

set.seed(2)
nosig <- matrix(runif(4 * 8, 0, 1), 4, 8,
                dimnames = list(paste0("S", 1:4), paste0("L", 1:8)))
nosig_sc <- score_matrices_by_species(nosig)
nosig_perm <- permute_omnibus(nosig_sc$z, nosig_sc$u, n_perm = 199L, seed = 123L, chunk_size = 50L)
check("no-signal Q p is not at the plus-one floor",
      nosig_perm$p_Q > 1 / 200)
check("signal Q p is smaller than no-signal Q p",
      sig_perm$p_Q < nosig_perm$p_Q)

# 9) BH and maxT
p_raw <- c(0.001, 0.04, 0.20, 0.80)
bh <- stats::p.adjust(p_raw, method = "BH")
fam <- bh_adjust_by_family(p_raw, rep("CDS", 4))
check("BH helper matches p.adjust on one family", almost_equal_vec(bh, fam))
check("maxT p is monotone in T",
      all(diff(sig_perm$p_maxT[order(-sig_perm$T_obs)]) >= -1e-12))
check("BH p exists only for K>=2 loci", all(is.finite(sig_perm$p_T_bh)))

# 10) 16 genus downsampling panels
combos <- genus_downsample_combos()
check("exactly 16 genus combinations", identical(length(combos), 16L))
keep_sets <- vapply(combos, function(x) paste(sort(x$kept_representatives), collapse = "|"), character(1))
check("16 combinations are unique", identical(length(unique(keep_sets)), 16L))
check("each combo keeps one Glycine/Oryza/Solanum/Gossypium",
      all(vapply(combos, function(x) {
        kept <- x$kept_representatives
        length(kept) == 4L &&
          any(grepl("^Glycine_", kept)) &&
          any(grepl("^Oryza_", kept)) &&
          any(grepl("^Solanum_", kept)) &&
          any(grepl("^Gossypium_", kept))
      }, logical(1))))
check("each combo drops four congeners",
      all(vapply(combos, function(x) length(x$dropped_species) == 4L, logical(1))))

# 11) canonical IGS direction and prefix
check("poiGS prefix stripped before canonical sort",
      identical(canonical_igs_id("poiGS_rbcL-atpB"), "atpB-rbcL"))
check("same-name IGS kept",
      identical(canonical_igs_id("poiGS_ndhA-ndhA"), "ndhA-ndhA"))
check("already sorted IGS unchanged after strip",
      identical(canonical_igs_id("poiGS_atpB-rbcL"), "atpB-rbcL"))

dup <- data.frame(
  species = c("A", "A"),
  locus_canonical = c("g1", "g1"),
  locus_raw = c("g1", "g1b"),
  stringsAsFactors = FALSE
)
conflicts <- audit_key_conflicts(dup)
check("duplicate species x locus is detected", nrow(conflicts) == 2L)
err <- try(stop_on_key_conflicts(dup, "CDS"), silent = TRUE)
check("duplicate keys halt formal analysis", inherits(err, "try-error"))

# 12) plus-one p and Monte Carlo SE
check("plus-one p with 0 exceedances is 1/(B+1)",
      almost_equal(empirical_p_plus_one(0, 99), 1 / 100))
check("MC SE formula",
      almost_equal(monte_carlo_se(0.5, 99), sqrt(0.5 * 0.5 / 100)))

# 13) frozen input gate and 25% baseline must remain 109/183
session_m03 <- file.path(workspace_root, "app_data", "sessions", "workflow_test_2025_11_4",
                         "results", "plots", "M03_hotspot", "M03_hotspot_standard",
                         "M03_all_genes_ranked_with_thresholds.csv")
session_m04 <- file.path(workspace_root, "app_data", "sessions", "workflow_test_2025_11_4",
                         "results", "plots", "M04_poigs_hotspot_engine", "M04_igs_hotspot_standard",
                         "M04_all_poigs_with_thresholds.csv")
if (file.exists(session_m03) && file.exists(session_m04)) {
  cds_raw <- read_csv_utf8(session_m03)
  igs_raw <- read_csv_utf8(session_m04)
  cds_std <- standardize_rank_table(cds_raw, "CDS")
  igs_std <- standardize_rank_table(igs_raw, "IGS")
  cds_sha <- file_sha256(session_m03)
  igs_sha <- file_sha256(session_m04)
  check("CDS sha256 matches freeze record", identical(cds_sha, FROZEN_SHA256$CDS))
  check("IGS sha256 matches freeze record", identical(igs_sha, FROZEN_SHA256$IGS))
  check("CDS long rows=1323", identical(nrow(cds_std), 1323L))
  check("CDS nonzero=410", identical(as.integer(sum(cds_std$frequency_per_kb > 0)), 410L))
  check("IGS long rows=1288", identical(nrow(igs_std), 1288L))
  check("IGS nonzero=701", identical(as.integer(sum(igs_std$frequency_per_kb > 0)), 701L))
  cds_mats <- build_species_locus_matrices(cds_std)
  igs_mats <- build_species_locus_matrices(igs_std)
  check("CDS analyzable cells=1323", identical(as.integer(cds_mats$n_analyzable), 1323L))
  check("IGS analyzable cells=1288", identical(as.integer(igs_mats$n_analyzable), 1288L))
  check("CDS missing cells are not filled with zero",
        cds_mats$n_missing > 0L && all(!is.finite(cds_mats$frequency[!cds_mats$analyzable])))
  check("no CDS key conflicts", nrow(audit_key_conflicts(cds_std)) == 0L)
  check("no IGS key conflicts", nrow(audit_key_conflicts(igs_std)) == 0L)
  check("all IGS ids have poiGS_ prefix", all(grepl("^poiGS_", igs_std$locus_raw)))
  igs_stripped <- unique(strip_poigs_prefix(igs_std$locus_raw))
  igs_reversed <- vapply(
    strsplit(igs_stripped, "-", fixed = TRUE),
    function(parts) paste(rev(parts), collapse = "-"),
    character(1)
  )
  check("IGS direction aliases do not collapse distinct loci",
        length(unique(igs_std$locus_canonical)) == length(igs_stripped))
  check("IGS table has no coexisting A-B and B-A aliases",
        !any(igs_reversed %in% igs_stripped & igs_reversed != igs_stripped))

  cds_hs <- cds_std
  cds_hs$locus_id <- cds_hs$locus_canonical
  igs_hs <- igs_std
  igs_hs$locus_id <- igs_hs$locus_canonical
  cds25 <- identify_hotspots_nonzero_quantile(cds_hs, 0.75, 7L)
  igs25 <- identify_hotspots_nonzero_quantile(igs_hs, 0.75, 7L)
  cds_sh <- summarize_sharing_classes(sharing_by_locus(cds25)$n_species)
  igs_sh <- summarize_sharing_classes(sharing_by_locus(igs25)$n_species)
  check("five-tier CDS 25% calls remain 109", sum(cds25$is_hotspot) == 109L)
  check("five-tier IGS 25% calls remain 183", sum(igs25$is_hotspot) == 183L)
  check("five-tier CDS 25% unique remain 50", cds_sh$n_unique == 50L)
  check("five-tier IGS 25% unique remain 63", igs_sh$n_unique == 63L)
  check("five-tier CDS private/shared2/shared3+ remain 24/6/20",
        cds_sh$private == 24L && cds_sh$shared2 == 6L && cds_sh$shared3plus == 20L)
  check("five-tier IGS private/shared2/shared3+ remain 21/10/32",
        igs_sh$private == 21L && igs_sh$shared2 == 10L && igs_sh$shared3plus == 32L)
} else {
  cat("SKIP optional frozen SHA / 25% regression (ranked tables not in this clone)\n")
}

cat("\n", n_pass, " passed, ", n_fail, " failed\n", sep = "")
if (n_fail > 0) {
  cat("Failures:\n- ", paste(failures, collapse = "\n- "), "\n", sep = "")
  quit(status = 1)
}
quit(status = 0)
