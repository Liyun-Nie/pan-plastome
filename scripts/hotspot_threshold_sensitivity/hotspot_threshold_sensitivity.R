# Hotspot threshold sensitivity (top 10/15/20/25/30%).
# Computational input is frozen session workflow_test_2025_11_4 only.
# New annotation outputs and Supplementary Table 12 are diagnostic, not inputs.

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(
    workspace = NULL,
    output = NULL,
    n_perm = NULL,
    seed = NULL,
    skip_perm = FALSE,
    config = NULL
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--skip-perm", "--skip_perm")) {
      out$skip_perm <- TRUE
      i <- i + 1L
    } else if (startsWith(key, "--") && i < length(args)) {
      val <- args[[i + 1L]]
      name <- sub("^--", "", key)
      name <- gsub("-", "_", name)
      out[[name]] <- val
      i <- i + 2L
    } else {
      stop("Unrecognized argument: ", key)
    }
  }
  out
}

cli <- parse_args(args)

this_file <- tryCatch({
  ofile <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(ofile)) {
    normalizePath(sub("^--file=", "", ofile[[1]]))
  } else {
    normalizePath("hotspot_threshold_sensitivity.R")
  }
}, error = function(e) {
  normalizePath(file.path("scripts", "hotspot_threshold_sensitivity",
                          "hotspot_threshold_sensitivity.R"))
})

mod_dir <- dirname(this_file)
scripts_root <- dirname(mod_dir)
pan_dir <- dirname(scripts_root)
workspace_root <- if (!is.null(cli$workspace)) {
  normalizePath(cli$workspace)
} else {
  # Standalone clone: session tree is supplied via --workspace
  normalizePath(pan_dir)
}

helper_path <- file.path(mod_dir, "hotspot_sensitivity_helpers.R")
source(helper_path, local = FALSE)

config_path <- if (!is.null(cli$config)) {
  cli$config
} else {
  file.path(pan_dir, "config", "hotspot_threshold_sensitivity.yml")
}
cfg <- read_simple_yaml_flat(config_path)

session_id <- if (!is.null(cfg$session_id)) cfg$session_id else "workflow_test_2025_11_4"
n_perm <- as.integer(if (!is.null(cli$n_perm)) cli$n_perm else if (!is.null(cfg$n_perm)) cfg$n_perm else 99999L)
seed <- as.integer(if (!is.null(cli$seed)) cli$seed else if (!is.null(cfg$seed)) cfg$seed else 20251104L)
alpha <- if (!is.null(cfg$alpha)) as.numeric(cfg$alpha) else 0.05
quantile_type <- if (!is.null(cfg$quantile_type)) as.integer(cfg$quantile_type) else 7L
skip_perm <- isTRUE(cli$skip_perm)

output_dir <- if (!is.null(cli$output)) {
  cli$output
} else {
  file.path(pan_dir, "results", "hotspot_threshold_sensitivity")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

session_root <- file.path(workspace_root, "app_data", "sessions", session_id)
m03_dir <- file.path(session_root, "results", "plots", "M03_hotspot", "M03_hotspot_standard")
m04_dir <- file.path(session_root, "results", "plots", "M04_poigs_hotspot_engine", "M04_igs_hotspot_standard")

paths <- list(
  m03_all = file.path(m03_dir, "M03_all_genes_ranked_with_thresholds.csv"),
  m03_hot = file.path(m03_dir, "M03_candidate_hotspots.csv"),
  m03_summary = file.path(m03_dir, "M03_candidate_hotspots_summary.csv"),
  m04_all = file.path(m04_dir, "M04_all_poigs_with_thresholds.csv"),
  m04_hot = file.path(m04_dir, "M04_candidate_poigs_hotspots.csv"),
  m04_freq = file.path(m04_dir, "M04_poigs_frequencies.csv"),
  m04_map = file.path(m04_dir, "M04_poigs_mapping.csv"),
  session_log = file.path(session_root, "logs", "session_2025-11-04_10-16-24.log"),
  p03_norm = file.path(session_root, "results", "processed_data", "P03_normalized", "normalized_frequencies.csv"),
  table12 = file.path(workspace_root, "manuscript", "pan-plastome", "version",
                      "first_revision_0525",
                      "Supplementary Table 12. Variant frequency for each IGS based on  SNV + indel + CPX.tsv"),
  new_m03 = file.path(workspace_root, "app_data", "sessions", "workflow_test_2026_01_05_3",
                      "results", "plots", "M03_hotspot", "M03_hotspot",
                      "M03_all_genes_ranked_with_thresholds.csv"),
  new_m04 = file.path(workspace_root, "app_data", "sessions", "workflow_test_2026_01_05_3",
                      "results", "plots", "M04_poigs_hotspot_engine", "M04_igs_hotspot",
                      "M04_all_poigs_with_thresholds.csv"),
  yml_ref = file.path(workspace_root, "project_management", "full_refactor_test.yml")
)

cat("Workspace: ", workspace_root, "\n", sep = "")
cat("Session: ", session_id, "\n", sep = "")
cat("Output: ", output_dir, "\n", sep = "")
cat("n_perm: ", n_perm, " seed: ", seed, " skip_perm: ", skip_perm, "\n", sep = "")

required_inputs <- c("m03_all", "m03_hot", "m04_all", "m04_hot")
missing_required <- vapply(required_inputs, function(k) !file.exists(paths[[k]]), logical(1))
if (any(missing_required)) {
  stop("Missing required session inputs:\n",
       paste(unlist(paths[required_inputs][missing_required]), collapse = "\n"))
}

manifest_all <- rbind(
  file_manifest_row(paths$m03_all, "cds_full_ranked", "PRIMARY computational input; CDS SNV frequencies + stored 25% flags"),
  file_manifest_row(paths$m03_hot, "cds_hotspot_25", "Session 25% CDS hotspot calls"),
  file_manifest_row(paths$m03_summary, "cds_hotspot_25_summary", "Session summary (109/50)"),
  file_manifest_row(paths$m04_all, "igs_full_ranked", "PRIMARY computational input; IGS aggregated SNV+INDEL+complex frequencies + stored 25% flags"),
  file_manifest_row(paths$m04_hot, "igs_hotspot_25", "Session 25% IGS hotspot calls"),
  file_manifest_row(paths$m04_freq, "igs_frequencies_upstream", "Upstream aggregated IGS frequency table; not used once ranked table is present"),
  file_manifest_row(paths$m04_map, "igs_mapping_later_mtime", "DIAGNOSTIC only; mtime is later than Nov 4 2025 session outputs"),
  file_manifest_row(paths$session_log, "session_log", "Generation log; frequency_percentile 0.75"),
  file_manifest_row(paths$p03_norm, "p03_normalized_later_mtime", "DIAGNOSTIC only; mtime later than Nov 4; NOT a computational input"),
  file_manifest_row(paths$table12, "supplementary_table12", "DIAGNOSTIC only; stored Is_Hotspot flags"),
  file_manifest_row(paths$new_m03, "new_annotation_m03", "DIAGNOSTIC only; new annotation CDS ranked table"),
  file_manifest_row(paths$new_m04, "new_annotation_m04", "DIAGNOSTIC only; new annotation IGS ranked table"),
  file_manifest_row(paths$yml_ref, "current_yml_reference", "Current YAML still records frequency_percentile 0.75; session log is the frozen generation record")
)
manifest <- manifest_all[
  manifest_all$role %in% c("cds_full_ranked", "igs_full_ranked"),
  ,
  drop = FALSE
]
manifest$path <- c(
  file.path(
    "app_data", "sessions", session_id, "results", "plots", "M03_hotspot",
    "M03_hotspot_standard", "M03_all_genes_ranked_with_thresholds.csv"
  ),
  file.path(
    "app_data", "sessions", session_id, "results", "plots",
    "M04_poigs_hotspot_engine", "M04_igs_hotspot_standard",
    "M04_all_poigs_with_thresholds.csv"
  )
)
write_tsv_utf8(manifest, file.path(output_dir, "00_input_manifest.tsv"))

m03_raw <- read_csv_utf8(paths$m03_all)
m04_raw <- read_csv_utf8(paths$m04_all)
m03_hot_raw <- read_csv_utf8(paths$m03_hot)
m04_hot_raw <- read_csv_utf8(paths$m04_hot)

cds <- standardize_frequency_table(m03_raw, "CDS", "gene")
igs <- standardize_frequency_table(m04_raw, "IGS", "poiGS_ID")

if (!all(tolower(m03_raw$region_type) == "cds")) {
  stop("M03 ranked table is not CDS-only")
}
if (!all(tolower(m03_raw$var_type) %in% c("snp", "snv"))) {
  stop("M03 ranked table is not SNV/snp-only")
}
if (!all(toupper(m04_raw$region_type) == "IGS")) {
  stop("M04 ranked table is not IGS-only")
}

run_one_cutoff <- function(freq_df, top_percent) {
  q <- quantile_prob_for_top_percent(top_percent)
  annotated <- identify_hotspots_nonzero_quantile(
    freq_df, quantile_prob = q, quantile_type = quantile_type
  )
  tests <- locus_test_table(annotated, alpha = alpha)
  tests$top_percent <- top_percent
  tests$quantile_prob <- q
  tests$family <- paste(tests$region_type, top_percent, sep = "|")
  tests$p_hyper_nonzero_bh <- bh_adjust_by_family(
    tests$p_hyper_nonzero_raw, tests$family
  )
  tests$sig_hyper_nonzero_bh <- tests$p_hyper_nonzero_bh < alpha
  if (skip_perm) {
    tests$p_perm_raw <- NA_real_
    tests$p_perm_bh <- NA_real_
    tests$sig_perm_bh <- NA
    tests$n_perm <- 0L
    tests$perm_seed <- seed
    tests$perm_eligible <- "nonzero_frequency"
    tests$methods_agree_bh <- NA
    null_distribution <- data.frame(
      region_type = character(),
      top_percent = numeric(),
      quantile_prob = numeric(),
      locus_id = character(),
      locus_display = character(),
      x_obs = integer(),
      x_null = integer(),
      n_permutations = integer(),
      probability = numeric(),
      p_perm_raw = numeric(),
      p_perm_bh = numeric(),
      sig_perm_bh = logical(),
      stringsAsFactors = FALSE
    )
  } else {
    perm_result <- stratified_permutation_pvalues(
      annotated,
      n_perm = n_perm,
      seed = seed,
      eligible_nonzero = TRUE,
      return_null_distribution = TRUE
    )
    perm <- perm_result$pvalues
    tests <- merge(tests, perm, by = "locus_id", all.x = TRUE, sort = FALSE)
    tests$p_perm_bh <- bh_adjust_by_family(tests$p_perm_raw, tests$family)
    tests$sig_perm_bh <- tests$p_perm_bh < alpha
    tests$methods_agree_bh <- (
      tests$sig_hyper_nonzero_bh == tests$sig_perm_bh
    )
    null_distribution <- merge(
      perm_result$null_distribution,
      tests[, c(
        "locus_id", "locus_display", "x", "p_perm_raw",
        "p_perm_bh", "sig_perm_bh"
      )],
      by = "locus_id",
      all.x = TRUE,
      sort = FALSE
    )
    names(null_distribution)[names(null_distribution) == "x"] <- "x_obs"
    null_distribution$region_type <- unique(annotated$region_type)[[1]]
    null_distribution$top_percent <- top_percent
    null_distribution$quantile_prob <- q
    null_distribution <- null_distribution[, c(
      "region_type", "top_percent", "quantile_prob",
      "locus_id", "locus_display", "x_obs", "x_null",
      "n_permutations", "probability", "p_perm_raw",
      "p_perm_bh", "sig_perm_bh"
    )]
  }
  tests$primary_inference <- "species_stratified_permutation"
  tests$tested_all_analyzable_loci <- TRUE
  list(
    annotated = annotated,
    summary = overall_summary_row(annotated, top_percent, q),
    per_species = per_species_summary(annotated, top_percent, q),
    tests = tests,
    calls = annotated[annotated$is_hotspot, , drop = FALSE],
    null_distribution = null_distribution
  )
}

cutoffs <- HOTSPOT_CUTOFF_MAP$top_percent
region_tables <- list(CDS = cds, IGS = igs)

all_summary <- list()
all_species <- list()
all_calls <- list()
all_tests <- list()
all_null_distributions <- list()
k <- 1L
for (rt in names(region_tables)) {
  for (tp in cutoffs) {
    cat("Running ", rt, " top ", tp, "%\n", sep = "")
    res <- run_one_cutoff(region_tables[[rt]], tp)
    all_summary[[k]] <- res$summary
    all_species[[k]] <- res$per_species
    call_df <- res$calls
    call_df$top_percent <- tp
    all_calls[[k]] <- call_df[, c("region_type", "top_percent", "species", "locus_id",
                                  "locus_display", "frequency_per_kb", "threshold",
                                  "n_nonzero", "n_hotspot_species",
                                  "actual_selected_fraction", "nominal_top_fraction")]
    all_tests[[k]] <- res$tests
    all_null_distributions[[k]] <- res$null_distribution
    k <- k + 1L
  }
}

summary_df <- do.call(rbind, all_summary)
species_df <- do.call(rbind, all_species)
calls_df <- do.call(rbind, all_calls)
tests_df <- do.call(rbind, all_tests)
null_distributions_df <- do.call(rbind, all_null_distributions)
tests_df <- tests_df[order(tests_df$region_type, tests_df$top_percent,
                           -tests_df$x, tests_df$p_perm_raw, tests_df$locus_id), ]
if (nrow(null_distributions_df) > 0L) {
  null_distributions_df <- null_distributions_df[
    order(
      null_distributions_df$region_type,
      null_distributions_df$top_percent,
      null_distributions_df$locus_id,
      null_distributions_df$x_null
    ),
  ]
}

# 25% regression vs session stored flags and manuscript targets
cds25 <- identify_hotspots_nonzero_quantile(cds, 0.75, quantile_type)
igs25 <- identify_hotspots_nonzero_quantile(igs, 0.75, quantile_type)

count_row <- function(annotated, source_label, expected) {
  sh <- sharing_by_locus(annotated)
  cl <- summarize_sharing_classes(sh$n_species)
  data.frame(
    source = source_label,
    region_type = unique(annotated$region_type)[[1]],
    n_calls = sum(annotated$is_hotspot),
    n_unique = cl$n_unique,
    n_private = cl$private,
    n_shared2 = cl$shared2,
    n_shared3plus = cl$shared3plus,
    expect_n_calls = expected$n_calls,
    expect_n_unique = expected$n_unique,
    expect_private = expected$private,
    expect_shared2 = expected$shared2,
    expect_shared3plus = expected$shared3plus,
    matches_expected = (
      sum(annotated$is_hotspot) == expected$n_calls &&
        cl$n_unique == expected$n_unique &&
        cl$private == expected$private &&
        cl$shared2 == expected$shared2 &&
        cl$shared3plus == expected$shared3plus
    ),
    stringsAsFactors = FALSE
  )
}

cds_stored <- cds
cds_stored$is_hotspot <- cds$stored_is_hotspot
igs_stored <- igs
igs_stored$is_hotspot <- igs$stored_is_hotspot

recon <- rbind(
  count_row(cds25, "session_recomputed_type7", MANUSCRIPT_BASELINE_25$CDS),
  count_row(cds_stored, "session_stored_flag", MANUSCRIPT_BASELINE_25$CDS),
  count_row(igs25, "session_recomputed_type7", MANUSCRIPT_BASELINE_25$IGS),
  count_row(igs_stored, "session_stored_flag", MANUSCRIPT_BASELINE_25$IGS)
)

# candidate files
m03_hot_std <- standardize_frequency_table(m03_hot_raw, "CDS", "gene")
m03_hot_std$is_hotspot <- TRUE
m04_hot_std <- standardize_frequency_table(m04_hot_raw, "IGS", "poiGS_ID")
m04_hot_std$is_hotspot <- TRUE
recon <- rbind(
  recon,
  count_row(m03_hot_std, "session_candidate_hotspots_file", MANUSCRIPT_BASELINE_25$CDS),
  count_row(m04_hot_std, "session_candidate_hotspots_file", MANUSCRIPT_BASELINE_25$IGS)
)

table12_row <- NULL
gb_detail <- NULL
if (file.exists(paths$table12)) {
  t12 <- read_tsv_utf8(paths$table12)
  names(t12) <- gsub("\\s+", "_", names(t12))
  t12$species <- species_key(t12$Species)
  t12$locus_display <- as.character(t12$IGS)
  freq_col <- grep("frequency", names(t12), ignore.case = TRUE, value = TRUE)[[1]]
  flag_col <- grep("hotspot", names(t12), ignore.case = TRUE, value = TRUE)[[1]]
  thr_col <- grep("threshold", names(t12), ignore.case = TRUE, value = TRUE)[[1]]
  t12$frequency_per_kb <- as.numeric(t12[[freq_col]])
  t12$stored_is_hotspot <- as_logical_flag(t12[[flag_col]])
  t12$stored_threshold <- as.numeric(t12[[thr_col]])
  t12_std <- data.frame(
    species = t12$species,
    locus_id = t12$locus_display,
    locus_display = t12$locus_display,
    region_type = "IGS",
    frequency_per_kb = t12$frequency_per_kb,
    stored_is_hotspot = t12$stored_is_hotspot,
    stringsAsFactors = FALSE
  )
  t12_re <- identify_hotspots_nonzero_quantile(t12_std, 0.75, quantile_type)
  t12_flag <- t12_std
  t12_flag$is_hotspot <- t12_std$stored_is_hotspot
  recon <- rbind(
    recon,
    count_row(t12_flag, "table12_stored_Is_Hotspot", MANUSCRIPT_BASELINE_25$IGS),
    count_row(t12_re, "table12_recomputed_type7", MANUSCRIPT_BASELINE_25$IGS)
  )

  gb_sess <- igs25[igs25$species == "Gossypium_barbadense", ]
  gb_t12 <- t12_re[t12_re$species == "Gossypium_barbadense", ]
  gb_t12_flag <- t12_flag[t12_flag$species == "Gossypium_barbadense", ]
  gb_detail <- data.frame(
    item = c(
      "session_GB_rows", "session_GB_nonzero", "session_GB_type7_threshold",
      "session_GB_hotspot_calls", "session_rpl16-rps3_freq", "session_rpl16-rps3_is_hotspot",
      "table12_GB_rows", "table12_GB_nonzero", "table12_GB_type7_recalc_threshold",
      "table12_GB_stored_threshold", "table12_GB_stored_hotspot_calls",
      "table12_GB_recalc_hotspot_calls", "table12_rpl16-rps3_stored_flag"
    ),
    value = c(
      nrow(gb_sess),
      sum(gb_sess$frequency_per_kb > 0),
      unique(gb_sess$threshold)[[1]],
      sum(gb_sess$is_hotspot),
      gb_sess$frequency_per_kb[gb_sess$locus_display == "rpl16-rps3"],
      gb_sess$is_hotspot[gb_sess$locus_display == "rpl16-rps3"],
      nrow(gb_t12),
      sum(gb_t12$frequency_per_kb > 0),
      unique(gb_t12$threshold)[[1]],
      unique(t12$stored_threshold[t12$species == "Gossypium_barbadense"])[[1]],
      sum(gb_t12_flag$is_hotspot),
      sum(gb_t12$is_hotspot),
      t12_flag$is_hotspot[t12_flag$species == "Gossypium_barbadense" & t12_flag$locus_display == "rpl16-rps3"]
    ),
    stringsAsFactors = FALSE
  )
}

if (file.exists(paths$new_m03) && file.exists(paths$new_m04)) {
  new_cds <- standardize_frequency_table(read_csv_utf8(paths$new_m03), "CDS", "gene")
  new_igs <- standardize_frequency_table(read_csv_utf8(paths$new_m04), "IGS", "poiGS_ID")
  new_cds25 <- identify_hotspots_nonzero_quantile(new_cds, 0.75, quantile_type)
  new_igs25 <- identify_hotspots_nonzero_quantile(new_igs, 0.75, quantile_type)
  recon <- rbind(
    recon,
    count_row(new_cds25, "new_annotation_recomputed_type7", MANUSCRIPT_BASELINE_25$CDS),
    count_row(new_igs25, "new_annotation_recomputed_type7", MANUSCRIPT_BASELINE_25$IGS)
  )
}

cds_gate <- isTRUE(recon$matches_expected[recon$source == "session_recomputed_type7" & recon$region_type == "CDS"])
igs_gate <- isTRUE(recon$matches_expected[recon$source == "session_recomputed_type7" & recon$region_type == "IGS"])
igs_session_internal <- (
  sum(igs25$is_hotspot) == SESSION_COMPUTED_25$IGS$n_calls &&
    summarize_sharing_classes(sharing_by_locus(igs25)$n_species)$n_unique == SESSION_COMPUTED_25$IGS$n_unique
)
formal_conclusion_allowed <- isTRUE(cds_gate && igs_gate)

disagree <- tests_df[!is.na(tests_df$methods_agree_bh) & !tests_df$methods_agree_bh, , drop = FALSE]
disagree_ge3 <- disagree[!is.na(disagree$shared_ge3) & disagree$shared_ge3, , drop = FALSE]

write_tsv_utf8(summary_df, file.path(output_dir, "01_overall_summary.tsv"))
write_tsv_utf8(species_df, file.path(output_dir, "02_per_species_thresholds.tsv"))
write_tsv_utf8(calls_df, file.path(output_dir, "03_hotspot_calls.tsv"))
calls25_df <- calls_df[calls_df$top_percent == 25, , drop = FALSE]
calls25_attachment_df <- data.frame(
  region_type = calls25_df$region_type,
  top_percent = calls25_df$top_percent,
  species = calls25_df$species,
  locus_id = calls25_df$locus_id,
  locus_display = calls25_df$locus_display,
  frequency_per_kb = calls25_df$frequency_per_kb,
  species_threshold = calls25_df$threshold,
  n_positive_candidate_loci_in_species = calls25_df$n_nonzero,
  n_hotspots_in_species = calls25_df$n_hotspot_species,
  actual_selected_fraction = calls25_df$actual_selected_fraction,
  nominal_top_fraction = calls25_df$nominal_top_fraction,
  stringsAsFactors = FALSE
)
write_tsv_utf8(
  calls25_attachment_df,
  file.path(output_dir, "03a_hotspot_calls_25pct.tsv")
)
primary_tests_df <- data.frame(
  region_type = tests_df$region_type,
  top_percent = tests_df$top_percent,
  quantile_prob = tests_df$quantile_prob,
  locus_id = tests_df$locus_id,
  locus_display = tests_df$locus_display,
  n_positive_frequency_species = tests_df$K_nonzero,
  observed_hotspot_species = tests_df$x,
  hotspot_species_list = tests_df$species_list,
  p_perm_raw = tests_df$p_perm_raw,
  p_perm_bh = tests_df$p_perm_bh,
  sig_perm_bh = tests_df$sig_perm_bh,
  n_perm = tests_df$n_perm,
  perm_seed = tests_df$perm_seed,
  perm_eligible = tests_df$perm_eligible,
  bh_family = tests_df$family,
  primary_inference = tests_df$primary_inference,
  tested_all_analyzable_loci = tests_df$tested_all_analyzable_loci,
  stringsAsFactors = FALSE
)
write_tsv_utf8(
  primary_tests_df,
  file.path(output_dir, "04_locus_sharing_tests.tsv")
)

hypergeom_dir <- file.path(output_dir, "supplementary_nonzero_hypergeometric")
dir.create(hypergeom_dir, recursive = TRUE, showWarnings = FALSE)
hypergeom_tests_df <- tests_df[, c(
  "region_type", "top_percent", "quantile_prob",
  "locus_id", "locus_display", "N_nonzero", "K_nonzero",
  "n_hotspot_calls", "x", "expected_x_nonzero",
  "enrichment_ratio_nonzero", "p_hyper_nonzero_raw",
  "p_hyper_nonzero_bh", "sig_hyper_nonzero_bh",
  "species_list", "hypergeom_role",
  "p_perm_raw", "p_perm_bh", "sig_perm_bh",
  "methods_agree_bh"
)]
write_tsv_utf8(
  hypergeom_tests_df,
  file.path(hypergeom_dir, "01_nonzero_pair_hypergeometric_tests.tsv")
)
writeLines(
  c(
    "# 非零组合超几何补充检验",
    "",
    "本子目录**不是**主要统计推断附件，也不应与上级目录的核心结果表并列提供给合作者。",
    "",
    "超几何模型把所有物种的正频率species×locus组合合并为一个抽样总体，不能保持每个物种自己的候选locus数和Hotspot数。它只用于判断结论对这种粗略总体假设是否敏感。与物种分层置换不一致时，**以物种分层置换为准**。",
    "",
    "核心统计附件仍为上级目录：",
    "",
    "- `../03a_hotspot_calls_25pct.tsv`：稿件25% Hotspot名单",
    "- `../04_locus_sharing_tests.tsv`：主要推断（物种分层置换）",
    "",
    "## 本目录文件",
    "",
    "- `01_nonzero_pair_hypergeometric_tests.tsv`：非零组合超几何原始及BH校正P值，并附带同一locus的置换P值以便对照",
    "",
    "## `01_nonzero_pair_hypergeometric_tests.tsv`列说明",
    "",
    "| 列名 | 含义 |",
    "| --- | --- |",
    "| `region_type` | `CDS`或`IGS` |",
    "| `top_percent` | 10、15、20、25或30 |",
    "| `quantile_prob` | 对应分位概率 |",
    "| `locus_id` / `locus_display` | 标准化标识和显示名称 |",
    "| `N_nonzero` | 当前区域类型、当前阈值下全部正频率species×locus组合数 |",
    "| `K_nonzero` | 该locus具有正频率的物种数 |",
    "| `n_hotspot_calls` | 当前阈值下全部Hotspot调用数 |",
    "| `x` | 该locus观察Hotspot次数 |",
    "| `expected_x_nonzero` | 超几何零模型下的期望次数 |",
    "| `enrichment_ratio_nonzero` | 观察次数除以期望次数 |",
    "| `p_hyper_nonzero_raw` | 超几何上尾原始P值 |",
    "| `p_hyper_nonzero_bh` | BH校正后的超几何P值 |",
    "| `sig_hyper_nonzero_bh` | 超几何BH P<0.05时为TRUE |",
    "| `species_list` | 实际将该locus判定为Hotspot的物种名单 |",
    "| `hypergeom_role` | 固定为`coarse_sensitivity_only` |",
    "| `p_perm_raw` / `p_perm_bh` / `sig_perm_bh` | 同一locus的主要物种分层置换结果 |",
    "| `methods_agree_bh` | 超几何BH显著性与置换BH显著性是否一致 |",
    "",
    "25% IGS下，超几何BH可能与置换BH在边界locus上不一致。正式表述以置换结果为准。"
  ),
  file.path(hypergeom_dir, "README_CN.md"),
  useBytes = FALSE
)
write_tsv_utf8(
  null_distributions_df,
  file.path(output_dir, "05_locus_permutation_null_distributions.tsv")
)

plot_threshold_null_panel <- function(null_df, test_row) {
  sub <- null_df[
    null_df$region_type == test_row$region_type &
      null_df$top_percent == test_row$top_percent &
      null_df$locus_id == test_row$locus_id,
    ,
    drop = FALSE
  ]
  sub <- sub[order(sub$x_null), , drop = FALSE]
  mids <- barplot(
    sub$probability,
    names.arg = sub$x_null,
    col = "gray88",
    border = "gray55",
    xlab = "Cross-species hotspot recurrence count",
    ylab = "Permutation probability",
    main = paste0(
      test_row$region_type, " Top ", test_row$top_percent,
      "%: ", test_row$locus_display
    ),
    cex.names = 0.75
  )
  obs_idx <- match(as.integer(test_row$x), sub$x_null)
  if (is.finite(obs_idx)) {
    abline(v = mids[obs_idx], col = "red", lwd = 2)
  }
  legend(
    "topright",
    legend = c(
      paste0("Permutation null (B=", format(n_perm, big.mark = ","), ")"),
      paste0("Observed recurrence = ", test_row$x)
    ),
    fill = c("gray88", NA),
    border = c("gray55", NA),
    lty = c(NA, 1),
    lwd = c(NA, 2),
    col = c("gray55", "red"),
    bty = "n",
    cex = 0.72
  )
  mtext(
    paste0(
      "raw P=", signif(test_row$p_perm_raw, 5),
      "; BH P=", signif(test_row$p_perm_bh, 5)
    ),
    side = 3,
    line = 0.2,
    cex = 0.72
  )
}

if (!skip_perm && nrow(null_distributions_df) > 0L) {
  focal <- (
    tests_df$region_type == "CDS" &
      tests_df$top_percent == 25 &
      tests_df$locus_display == "ycf1"
  ) | (
    tests_df$region_type == "IGS" &
      tests_df$top_percent == 25 &
      tests_df$locus_display == "rpoC2-rps2"
  )
  plot_rows <- tests_df[
    (!is.na(tests_df$sig_perm_bh) & tests_df$sig_perm_bh) | focal,
    ,
    drop = FALSE
  ]
  plot_rows$plot_priority <- ifelse(focal[
    match(
      paste(plot_rows$region_type, plot_rows$top_percent, plot_rows$locus_id),
      paste(tests_df$region_type, tests_df$top_percent, tests_df$locus_id)
    )
  ], 0L, 1L)
  plot_rows <- plot_rows[
    order(
      plot_rows$plot_priority,
      plot_rows$region_type,
      plot_rows$top_percent,
      plot_rows$p_perm_bh,
      plot_rows$locus_id
    ),
    ,
    drop = FALSE
  ]
  pdf(
    file.path(output_dir, "threshold_locus_observed_vs_null_comparison.pdf"),
    width = 11,
    height = 8.5
  )
  par(mfrow = c(2, 2), mar = c(4.2, 4.4, 3.3, 1.1))
  for (i in seq_len(nrow(plot_rows))) {
    plot_threshold_null_panel(null_distributions_df, plot_rows[i, ])
  }
  dev.off()
}

summary_report_lines <- function(region_type) {
  sub <- summary_df[summary_df$region_type == region_type, , drop = FALSE]
  vapply(seq_len(nrow(sub)), function(i) {
    r <- sub[i, ]
    paste0(
      "- Top ", r$top_percent, "%：", r$n_hotspot_calls,
      " calls，", r$n_unique_loci, "个locus；private/shared-by-2/shared-by-3+ = ",
      r$n_private, "/", r$n_shared2, "/", r$n_shared3plus,
      "；最大共享", r$max_shared_species, "个物种。"
    )
  }, character(1))
}

significant_report_lines <- function(region_type, column) {
  sub <- tests_df[
    tests_df$region_type == region_type &
      !is.na(tests_df[[column]]) & tests_df[[column]],
    ,
    drop = FALSE
  ]
  if (!nrow(sub)) {
    return("- 无。")
  }
  by_cutoff <- split(sub, sub$top_percent)
  vapply(names(by_cutoff), function(tp) {
    z <- by_cutoff[[tp]]
    paste0("- Top ", tp, "%：", paste(z$locus_display, collapse = "、"), "。")
  }, character(1))
}

key_result_line <- function(region_type, locus_display, top_percent) {
  z <- tests_df[
    tests_df$region_type == region_type &
      tests_df$locus_display == locus_display &
      tests_df$top_percent == top_percent,
    ,
    drop = FALSE
  ]
  if (!nrow(z)) {
    return(NULL)
  }
  z <- z[1, ]
  paste0(
    "- `", z$locus_display, "`：x=", z$x,
    "；分层置换raw P=", signif(z$p_perm_raw, 6),
    "，BH P=", signif(z$p_perm_bh, 6),
    "；非零pair超几何raw P=", signif(z$p_hyper_nonzero_raw, 6),
    "，BH P=", signif(z$p_hyper_nonzero_bh, 6), "。"
  )
}

perm_sig_cutoffs <- function(region_type) {
  sub <- tests_df[
    tests_df$region_type == region_type &
      !is.na(tests_df$sig_perm_bh) & tests_df$sig_perm_bh,
    ,
    drop = FALSE
  ]
  sort(unique(as.integer(sub$top_percent)))
}

cutoff_significance_sentence <- function() {
  cds_sig <- perm_sig_cutoffs("CDS")
  igs_sig <- perm_sig_cutoffs("IGS")
  none_25_30 <- !any(c(25L, 30L) %in% c(cds_sig, igs_sig))
  has_10_20 <- any(c(10L, 15L, 20L) %in% c(cds_sig, igs_sig))
  lead <- paste0(
    "阈值从Top 10%放宽到Top 30%时，Hotspot调用数和共享locus数平稳增加，没有出现天然分界点；",
    "25%是操作性选择，不是数据自动产生的唯一阈值。"
  )
  if (has_10_20 && none_25_30) {
    return(paste0(
      lead,
      "主要物种分层置换在Top 10%至Top 20%发现少数显著复现locus，",
      "但Top 25%和Top 30%均无多重校正后显著的单个locus。"
    ))
  }
  fmt <- function(v) {
    if (!length(v)) {
      return("无")
    }
    paste0("Top ", paste(v, collapse = "%、Top "), "%")
  }
  paste0(
    lead,
    "主要物种分层置换在CDS的", fmt(cds_sig),
    "以及IGS的", fmt(igs_sig),
    "发现多重校正后显著的单个locus。"
  )
}

lookup_test_row <- function(region_type, locus_display, top_percent) {
  z <- tests_df[
    tests_df$region_type == region_type &
      tests_df$locus_display == locus_display &
      tests_df$top_percent == top_percent,
    ,
    drop = FALSE
  ]
  if (!nrow(z)) {
    return(NULL)
  }
  z[1, ]
}

key_25_narrative <- function() {
  y <- lookup_test_row("CDS", "ycf1", 25)
  r <- lookup_test_row("IGS", "rpoC2-rps2", 25)
  y_txt <- if (is.null(y)) {
    "`ycf1`未在25% CDS结果中找到。"
  } else if (isTRUE(y$sig_perm_bh)) {
    paste0(
      "`ycf1`在", y$x, "个物种中成为25% CDS Hotspot，物种分层置换BH P=",
      signif(y$p_perm_bh, 6), "，达到多重校正后显著。"
    )
  } else {
    paste0(
      "`ycf1`在", y$x, "个物种中成为25% CDS Hotspot，",
      "但其共享程度在保持各物种候选规模和Hotspot数的随机情形中并不罕见，",
      "不能写成显著富集。"
    )
  }
  r_txt <- if (is.null(r)) {
    "`rpoC2-rps2`未在25% IGS结果中找到。"
  } else if (isTRUE(r$sig_perm_bh)) {
    paste0(
      "`rpoC2-rps2`在", r$x, "个物种中成为25% IGS Hotspot，物种分层置换BH P=",
      signif(r$p_perm_bh, 6), "，达到多重校正后显著。"
    )
  } else if (is.finite(r$p_perm_bh) && r$p_perm_bh > 0.05 && r$p_perm_bh < 0.10) {
    paste0(
      "`rpoC2-rps2`在", r$x, "个物种中成为25% IGS Hotspot，",
      "其物种分层置换BH P=", signif(r$p_perm_bh, 6),
      "，只能表述为接近显著但证据不足。"
    )
  } else {
    paste0(
      "`rpoC2-rps2`在", r$x, "个物种中成为25% IGS Hotspot，物种分层置换BH P=",
      signif(r$p_perm_bh, 6), "，不能写成显著富集。"
    )
  }
  cds25 <- any(
    tests_df$region_type == "CDS" &
      tests_df$top_percent == 25 &
      !is.na(tests_df$sig_perm_bh) & tests_df$sig_perm_bh
  )
  igs25 <- any(
    tests_df$region_type == "IGS" &
      tests_df$top_percent == 25 &
      !is.na(tests_df$sig_perm_bh) & tests_df$sig_perm_bh
  )
  wrap <- if (!cds25 && !igs25) {
    paste0(
      "25%下两类区域均无物种分层置换BH<0.05的单locus。",
      "该结果只表示二元Top 25%复现未通过当前检验，不否定其他阈值的显著结果，",
      "也不否定无阈值连续秩检验检测到的分布式共同结构。"
    )
  } else if (cds25 && igs25) {
    paste0(
      "25%下CDS和IGS均有物种分层置换BH<0.05的单locus。",
      "该结果只针对二元Top 25%复现，不否定其他阈值的显著结果，",
      "也不否定无阈值连续秩检验检测到的分布式共同结构。"
    )
  } else if (igs25) {
    paste0(
      "25%下IGS有物种分层置换BH<0.05的单locus，CDS则无。",
      "该结果只针对二元Top 25%复现，不否定其他阈值的显著结果，",
      "也不否定无阈值连续秩检验检测到的分布式共同结构。"
    )
  } else {
    paste0(
      "25%下CDS有物种分层置换BH<0.05的单locus，IGS则无。",
      "该结果只针对二元Top 25%复现，不否定其他阈值的显著结果，",
      "也不否定无阈值连续秩检验检测到的分布式共同结构。"
    )
  }
  list(combo = paste0(y_txt, r_txt), wrap = wrap)
}

narrative_25 <- key_25_narrative()

report_lines <- c(
  "# 有阈值Hotspot敏感性分析中文结果报告",
  "",
  paste0("**冻结输入**：`", session_id, "`"),
  paste0("**主要统计检验**：物种分层置换（每个区域类型和阈值", n_perm, "次）"),
  "**补充敏感性检验**：非零物种与locus组合超几何检验",
  "**最终状态**：分析和冻结数据回归均通过。",
  "",
  "## 一、先看结论",
  "",
  "本分析先在每个物种内部的正频率locus中定义变异频率最高的Top 10%、15%、20%、25%和30%为Hotspot，再检验同一locus跨物种反复成为Hotspot的次数是否高于随机预期。",
  "",
  cutoff_significance_sentence(),
  "",
  "## 二、分析对象及为什么排除零频率",
  "",
  "CDS仅使用单核苷酸变异频率；IGS使用单核苷酸变异、插入缺失和复杂变异的合计频率。所有频率均按区域长度换算为每千碱基变异数。",
  "",
  "有阈值分析回答的是：在某个物种已经发生变异的locus中，哪些属于最高频尾部？因此，每个物种仅把`frequency_per_kb>0`的locus纳入阈值候选集。零频率表示该物种在该locus未检测到变异，按当前定义不可能成为Hotspot；如果把大量零值纳入分位数且阈值降至零，`frequency>=threshold`可能把零频率locus错误纳入Hotspot。",
  "",
  "这种排除不是所有阈值分析都必须遵守的统计定理，而是本研究将Hotspot定义为“正频率候选集中的高频尾部”所决定的条件分析。无阈值分析回答完整频率景观问题，因此保留真实零值；两者检验目标不同。",
  "",
  "## 三、每个物种如何计算Hotspot阈值",
  "",
  "设物种`s`有`m_s`个正频率候选locus，其升序频率为`y_(1),...,y_(m_s)`。名义Top比例`q`对应分位概率`p=1-q`。这里使用R函数`quantile(..., type=7)`，即R默认的第7型样本分位数估计方法；“第7型”是Hyndman–Fan对九种样本分位数算法的编号，不是“第七个分位点”。该方法计算`h=1+(m_s-1)p`；令`j=floor(h)`和`γ=h-j`，则阈值为`c_s(p)=(1-γ)y_(j)+γy_(j+1)`，也就是在相邻两个有序观测值之间作线性插值。",
  "",
  "若正频率locus `l`满足`y_sl>=c_s(p)`，则定义`H_sl(q)=1`，否则为0。阈值处并列值全部保留，因此实际Hotspot比例可能略高于名义Top比例。Top 10%、15%、20%、25%和30%分别对应0.90、0.85、0.80、0.75和0.70分位数。",
  "",
  "`private`、`shared-by-2`和`shared-by-3+`分别表示某个locus在1个、2个或至少3个物种中成为Hotspot。这些类别只描述观察共享程度，不参与显著性定义。",
  "",
  "## 四、跨物种复现统计量与物种分层置换",
  "",
  "对于locus `l`，观察到的跨物种Hotspot复现次数为`X_l=Σ_s H_sl`。`X_l=6`表示该locus在6个物种中被定义为当前阈值下的Hotspot。",
  "",
  "零假设为：在保持每个物种自己的正频率候选集和实际Hotspot数量后，Hotspot标签可在该物种候选locus之间随机交换，locus身份不产生额外的跨物种共同复现。每次置换执行以下步骤：",
  "",
  "1. 对17个物种分别处理，物种之间不交换数据；",
  "2. 保持物种`s`的正频率候选集`E_s`及实际Hotspot数`h_s`不变；",
  "3. 在`E_s`中无放回随机抽取恰好`h_s`个locus并赋予Hotspot标签；",
  "4. 完成所有物种后，重新计算每个locus的置换复现次数`X_l^(b)`；",
  paste0("5. 重复", n_perm, "次，形成每个locus自己的离散置换零分布；"),
  "6. 以`P_l=(1+b_l)/(1+B)`计算上尾经验P值，其中`b_l`为置换复现次数不小于观察`X_l`的次数，`B`为置换总数。",
  "",
  "因此，本分析不是检验随机频率能否超过观察阈值，而是在阈值已经把数据二值化后，检验相同locus跨物种反复获得Hotspot标签的次数是否异常。原始经验P值随后在每个“区域类型×阈值”检验族内进行Benjamini–Hochberg错误发现率校正。",
  "",
  "## 五、观察值与置换零分布可视化",
  "",
  "`05_locus_permutation_null_distributions.tsv`保存每个区域类型、阈值和locus在0至17次随机复现上的置换计数与概率，可完整重建所有逐locus零分布。",
  "",
  "`threshold_locus_observed_vs_null_comparison.pdf`绘制所有物种分层置换BH显著结果，并额外纳入25%阈值下用于审查的`ycf1`和`rpoC2-rps2`。灰色柱表示置换零分布，红线表示真实观察复现次数，图中同时标注原始及BH校正P值。单个真实数据集只产生一个观察复现次数，因此比较的是“一个观察值与一个置换零分布”，不存在独立的观察分布。",
  "",
  "## 六、五档描述性结果",
  "",
  "### CDS",
  "",
  summary_report_lines("CDS"),
  "",
  "### IGS",
  "",
  summary_report_lines("IGS"),
  "",
  "25%冻结基线为CDS 109 calls/50个locus和IGS 183 calls/63个locus。旧IGS 184来自`Gossypium barbadense rpl16-rps3`的边界标记错误；冻结session重新计算及存储标记均确认最终基线为183。",
  "",
  "这里的调用次数是物种与locus组合数。同一locus在5个物种中入选会贡献5次调用，但只算1个不同locus。阈值放宽必然产生更多调用和共享，因此共享数量增加本身不能证明富集。",
  "",
  "## 七、物种分层置换结果（主要推断）",
  "",
  "### CDS中BH显著的locus",
  "",
  significant_report_lines("CDS", "sig_perm_bh"),
  "",
  "### IGS中BH显著的locus",
  "",
  significant_report_lines("IGS", "sig_perm_bh"),
  "",
  "### 25%关键结果",
  "",
  key_result_line("CDS", "ycf1", 25),
  key_result_line("IGS", "rpoC2-rps2", 25),
  "",
  narrative_25$combo,
  "",
  narrative_25$wrap,
  "",
  "## 八、非零组合超几何检验（补充敏感性）",
  "",
  "为回应合并总体定义的敏感性，另外计算非零物种与locus组合超几何检验。设`N_+`为全部正频率组合数，`K_l+`为目标locus具有正频率的物种数，`n`为当前阈值的全部Hotspot调用数，`x_l`为目标locus观察Hotspot次数，则`X_l~Hypergeometric(N_+,K_l+,n)`，上尾P值为`Pr(X_l>=x_l)`。",
  "",
  "该模型排除了零频率组合，但把不同物种合并为一个抽样池，不能保持各物种不同的候选数和Hotspot数。因此，它仅用于判断结论对粗略总体假设是否敏感，不承担主要推断；与物种分层置换不一致时，以后者为准。",
  "",
  "### CDS中BH显著的locus",
  "",
  significant_report_lines("CDS", "sig_hyper_nonzero_bh"),
  "",
  "### IGS中BH显著的locus",
  "",
  significant_report_lines("IGS", "sig_hyper_nonzero_bh"),
  "",
  "凡非零pair超几何与分层置换结论不一致，以分层置换为主。核心逐locus置换结果保存在`04_locus_sharing_tests.tsv`；超几何结果已单独放入`supplementary_nonzero_hypergeometric/01_nonzero_pair_hypergeometric_tests.tsv`，不与核心附件并列。",
  "",
  "## 九、与无阈值检验的关系",
  "",
  "两套分析都在物种内部建立零模型，再按locus身份跨物种汇总。有阈值分析把正频率locus转换为Hotspot与非Hotspot，检验特定界线上的重复越界；无阈值分析保留全部可分析locus及其完整物种内排名，检验多个locus共同形成的连续高排名结构。",
  "",
  "所以，Top 25%没有显著单locus与无阈值`Q_high`显著并不矛盾。前者要求某一个locus在特定二元界线上足够极端；后者能够累积多个locus方向一致但不一定每次越过25%界线的信息。",
  "",
  "## 十、统计依据与应用先例",
  "",
  "第7类样本分位数、条件置换、有限置换P值和错误发现率校正均有成熟统计依据。保持每个样本特异事件负担后检验基因组位置跨样本复发的思想也已用于GISTIC及泛癌拷贝数变异研究；这些研究支持分层内置换的设计原则，但不代表使用了与本研究完全相同的Hotspot定义和统计公式。",
  "",
  "1. Hyndman RJ, Fan Y. Sample Quantiles in Statistical Packages. *The American Statistician*. 1996;50:361–365. doi:10.1080/00031305.1996.10473566.",
  "2. Strasser H, Weber C. On the asymptotic theory of permutation statistics. *Mathematical Methods of Statistics*. 1999;8:220–250.",
  "3. Phipson B, Smyth GK. Permutation P-values Should Never Be Zero. *Statistical Applications in Genetics and Molecular Biology*. 2010;9:Article 39. doi:10.2202/1544-6115.1585.",
  "4. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate. *Journal of the Royal Statistical Society Series B*. 1995;57:289–300. doi:10.1111/j.2517-6161.1995.tb02031.x.",
  "5. Beroukhim R, Getz G, Nghiemphu L, et al. Assessing the significance of chromosomal aberrations in cancer. *PNAS*. 2007;104:20007–20012. doi:10.1073/pnas.0710052104.",
  "6. Beroukhim R, Mermel CH, Porter D, et al. The landscape of somatic copy-number alteration across human cancers. *Nature*. 2010;463:899–905. doi:10.1038/nature08822.",
  "7. Zack TI, Schumacher SE, Carter SL, et al. Pan-cancer patterns of somatic copy number alteration. *Nature Genetics*. 2013;45:1134–1140. doi:10.1038/ng.2760.",
  "",
  "## 十一、解释边界",
  "",
  "- 25%及其他百分比是操作性阈值，不是自然断点。",
  "- core/private等类别是观察后描述，不能用于循环证明共享显著。",
  "- 分层置换检验的是当前17个物种中的条件随机复现，不等同于17次完全独立的进化重复。",
  "- 本分析不检验染色体物理聚集，也不识别长度、功能约束或系统发育等具体因果机制。",
  "- 完整逐locus结果以`04_locus_sharing_tests.tsv`为准，完整置换零分布以`05_locus_permutation_null_distributions.tsv`为准，两套分析的联合技术流程见`../HOTSPOT_COMBINED_TECHNICAL_FLOWCHART_XMIND_CN.md`。"
)
writeLines(
  report_lines[!is.na(report_lines)],
  file.path(output_dir, "RESULTS_REPORT_CN.md"),
  useBytes = FALSE
)

meta <- c(
  paste0("session_id=", session_id),
  paste0("n_perm=", n_perm),
  paste0("seed=", seed),
  paste0("skip_perm=", skip_perm),
  paste0("alpha=", alpha),
  paste0("quantile_type=", quantile_type),
  paste0("cds_gate_109=", cds_gate),
  paste0("igs_approved_gate_183=", igs_gate),
  paste0("igs_session_internal_183=", igs_session_internal),
  paste0("formal_conclusion_allowed=", formal_conclusion_allowed),
  "primary_inference=species_stratified_permutation",
  "hypergeometric_universe=nonzero_pairs",
  "hypergeometric_role=coarse_sensitivity_only",
  paste0("n_method_disagreements=", nrow(disagree)),
  paste0("n_method_disagreements_shared_ge3=", nrow(disagree_ge3)),
  paste0("generated_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
)
writeLines(meta, file.path(output_dir, "run_metadata.txt"))

cat("\n=== DONE ===\n")
cat("CDS 25% manuscript gate: ", cds_gate, "\n", sep = "")
cat("IGS 25% approved 183 gate: ", igs_gate, "\n", sep = "")
cat("formal_conclusion_allowed: ", formal_conclusion_allowed, "\n", sep = "")
print(summary_df)
if (!formal_conclusion_allowed) {
  cat("\nBLOCKER: the approved 25% baseline was not reproduced. Diagnostic written; no manuscript text updated.\n")
}
