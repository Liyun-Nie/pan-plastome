# Threshold-free hotspot concordance analysis.
# Primary input: frozen session workflow_test_2025_11_4 ranked CDS/IGS tables.
# Does not modify cpopvar/, the archived release snapshot, or the historical session.

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(
    workspace = NULL,
    output = NULL,
    n_perm = NULL,
    n_perm_robust = NULL,
    seed = NULL,
    mode = NULL,
    chunk_size = NULL,
    config = NULL,
    skip_robustness = FALSE,
    skip_perm = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--skip-robustness", "--skip_robustness")) {
      out$skip_robustness <- TRUE
      i <- i + 1L
    } else if (key %in% c("--skip-perm", "--skip_perm")) {
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
    normalizePath("hotspot_threshold_free_analysis.R")
  }
}, error = function(e) {
  normalizePath(file.path("scripts", "hotspot_threshold_free_analysis",
                          "hotspot_threshold_free_analysis.R"))
})

mod_dir <- dirname(this_file)
scripts_root <- dirname(mod_dir)
pan_dir <- dirname(scripts_root)
workspace_root <- if (!is.null(cli$workspace)) {
  normalizePath(cli$workspace)
} else {
  normalizePath(pan_dir)
}

source(file.path(scripts_root, "hotspot_threshold_sensitivity",
                 "hotspot_sensitivity_helpers.R"), local = FALSE)
source(file.path(mod_dir, "hotspot_threshold_free_helpers.R"), local = FALSE)
assert_no_cpopvar_attached()

config_path <- if (!is.null(cli$config)) {
  cli$config
} else {
  file.path(pan_dir, "config", "hotspot_threshold_free.yml")
}
cfg <- read_simple_yaml_flat(config_path)

session_id <- if (!is.null(cfg$session_id)) cfg$session_id else "workflow_test_2025_11_4"
run_mode <- if (!is.null(cli$mode)) cli$mode else if (!is.null(cfg$run_mode)) cfg$run_mode else "official"
n_perm_official <- as.integer(if (!is.null(cfg$n_perm)) cfg$n_perm else 99999L)
n_perm_dev <- as.integer(if (!is.null(cfg$n_perm_dev)) cfg$n_perm_dev else 9999L)
n_perm <- as.integer(if (!is.null(cli$n_perm)) {
  cli$n_perm
} else if (identical(run_mode, "dev")) {
  n_perm_dev
} else {
  n_perm_official
})
n_perm_robust <- as.integer(if (!is.null(cli$n_perm_robust)) {
  cli$n_perm_robust
} else if (!is.null(cfg$n_perm_robust)) {
  cfg$n_perm_robust
} else {
  n_perm_dev
})
seed <- as.integer(if (!is.null(cli$seed)) cli$seed else if (!is.null(cfg$seed)) cfg$seed else 20260915L)
alpha <- if (!is.null(cfg$alpha)) as.numeric(cfg$alpha) else 0.05
chunk_size <- as.integer(if (!is.null(cli$chunk_size)) cli$chunk_size else if (!is.null(cfg$chunk_size)) cfg$chunk_size else 250L)
k_min <- as.integer(if (!is.null(cfg$k_min_cross_species)) cfg$k_min_cross_species else 2L)
min_overlap <- as.integer(if (!is.null(cfg$min_overlap_cglobal)) cfg$min_overlap_cglobal else 3L)
run_robustness <- isTRUE(cfg$run_robustness) && !isTRUE(cli$skip_robustness)
skip_perm <- isTRUE(cli$skip_perm)

output_dir <- if (!is.null(cli$output)) {
  cli$output
} else {
  file.path(pan_dir, "results", "hotspot_threshold_free_analysis")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

session_root <- file.path(workspace_root, "app_data", "sessions", session_id)
m03_path <- file.path(session_root, "results", "plots", "M03_hotspot", "M03_hotspot_standard",
                      "M03_all_genes_ranked_with_thresholds.csv")
m04_path <- file.path(session_root, "results", "plots", "M04_poigs_hotspot_engine",
                      "M04_igs_hotspot_standard", "M04_all_poigs_with_thresholds.csv")

cat("Workspace: ", workspace_root, "\n", sep = "")
cat("Session: ", session_id, "\n", sep = "")
cat("Output: ", output_dir, "\n", sep = "")
cat("mode: ", run_mode, " n_perm: ", n_perm, " n_perm_robust: ", n_perm_robust,
    " seed: ", seed, "\n", sep = "")

if (!file.exists(m03_path) || !file.exists(m04_path)) {
  stop("Missing frozen ranked tables")
}

cds_sha <- file_sha256(m03_path)
igs_sha <- file_sha256(m04_path)
expected_cds_sha <- if (!is.null(cfg$expected_cds_sha256)) cfg$expected_cds_sha256 else FROZEN_SHA256$CDS
expected_igs_sha <- if (!is.null(cfg$expected_igs_sha256)) cfg$expected_igs_sha256 else FROZEN_SHA256$IGS
if (!identical(cds_sha, expected_cds_sha)) {
  stop("CDS SHA-256 mismatch: got ", cds_sha, " expected ", expected_cds_sha)
}
if (!identical(igs_sha, expected_igs_sha)) {
  stop("IGS SHA-256 mismatch: got ", igs_sha, " expected ", expected_igs_sha)
}

fmt_p <- function(x) {
  ifelse(is.finite(x), format(x, scientific = TRUE, digits = 6), NA)
}

empty_perm <- function(z_mat, u_mat, k_min, min_overlap) {
  obs_T <- locus_T(z_mat, k_min = k_min)
  n_l <- length(obs_T$T)
  pair_obs <- pairwise_rank_concordance(u_mat, min_overlap = min_overlap)$pairs
  pair_keys <- paste(pair_obs$species_a, pair_obs$species_b, sep = "\t")
  list(
    n_perm = 0L,
    seed = NA_integer_,
    k_min = k_min,
    T_obs = obs_T$T,
    K = obs_T$K,
    Q_obs = Q_high(obs_T$T),
    C_obs = C_global(u_mat, min_overlap = min_overlap),
    M_obs = M_max(obs_T$T),
    n_ge_Q = NA_integer_,
    n_ge_C = NA_integer_,
    n_ge_M = NA_integer_,
    pair_keys = pair_keys,
    pair_obs = setNames(pair_obs$r_midrank_pearson, pair_keys),
    n_ge_pair = setNames(rep(NA_integer_, nrow(pair_obs)), pair_keys),
    p_pair = setNames(rep(NA_real_, nrow(pair_obs)), pair_keys),
    se_pair = setNames(rep(NA_real_, nrow(pair_obs)), pair_keys),
    n_ge_T = rep(NA_integer_, n_l),
    n_ge_maxT = rep(NA_integer_, n_l),
    p_Q = NA_real_,
    p_C = NA_real_,
    p_M = NA_real_,
    p_T = rep(NA_real_, n_l),
    p_T_bh = rep(NA_real_, n_l),
    p_maxT = rep(NA_real_, n_l),
    se_Q = NA_real_,
    se_C = NA_real_,
    se_M = NA_real_,
    se_T = rep(NA_real_, n_l),
    se_maxT = rep(NA_real_, n_l),
    Q_null = numeric(),
    C_null = numeric(),
    M_null = numeric()
  )
}

descriptive_threshold_overlay <- function(std_df, locus_tbl) {
  tmp <- std_df
  tmp$locus_id <- tmp$locus_canonical
  out <- list()
  k <- 1L
  for (i in seq_len(nrow(HOTSPOT_CUTOFF_MAP))) {
    tp <- HOTSPOT_CUTOFF_MAP$top_percent[i]
    q <- HOTSPOT_CUTOFF_MAP$quantile_prob[i]
    ann <- identify_hotspots_nonzero_quantile(tmp, q, 7L)
    sharing <- sharing_by_locus(ann)
    classes <- summarize_sharing_classes(sharing$n_species)
    merged <- merge(
      data.frame(locus_canonical = locus_tbl$locus_canonical,
                 T_l = locus_tbl$T_l,
                 stringsAsFactors = FALSE),
      sharing,
      by.x = "locus_canonical",
      by.y = "locus_id",
      all.x = TRUE
    )
    merged$n_species[is.na(merged$n_species)] <- 0L
    out[[k]] <- data.frame(
      region_type = unique(std_df$region_type)[[1]],
      top_percent = tp,
      quantile_prob = q,
      locus_canonical = merged$locus_canonical,
      T_l = merged$T_l,
      n_hotspot_species_descriptive = as.integer(merged$n_species),
      hotspot_species_list_descriptive = ifelse(is.na(merged$species_list), "", merged$species_list),
      n_hotspot_calls_cutoff = sum(ann$is_hotspot),
      n_unique_hotspot_loci_cutoff = classes$n_unique,
      n_private_cutoff = classes$private,
      n_shared2_cutoff = classes$shared2,
      n_shared3plus_cutoff = classes$shared3plus,
      used_in_primary_inference = FALSE,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  do.call(rbind, out)
}

analyze_region <- function(region_type, raw_df, path, sha) {
  std <- standardize_rank_table(raw_df, region_type)
  expected_rows <- if (region_type == "CDS") as.integer(cfg$expected_cds_rows) else as.integer(cfg$expected_igs_rows)
  expected_nz <- if (region_type == "CDS") as.integer(cfg$expected_cds_nonzero) else as.integer(cfg$expected_igs_nonzero)
  if (!identical(nrow(std), expected_rows)) {
    stop(region_type, " row count ", nrow(std), " != freeze gate ", expected_rows)
  }
  n_pos <- as.integer(sum(std$frequency_per_kb > 0))
  if (!identical(n_pos, expected_nz)) {
    stop(region_type, " nonzero count ", n_pos, " != freeze gate ", expected_nz)
  }
  mats <- build_species_locus_matrices(std)
  scores <- score_matrices_by_species(mats$frequency)
  audit <- audit_rows(region_type, path, sha, std, mats, scores, k_min = k_min)
  long_rank <- rank_matrix_long(std, mats, scores, region_type)
  pair_obs <- pairwise_rank_concordance(scores$u, min_overlap = min_overlap)
  pair_obs$pairs$region_type <- region_type
  t0 <- proc.time()[["elapsed"]]
  if (skip_perm) {
    perm <- empty_perm(scores$z, scores$u, k_min, min_overlap)
  } else {
    cat("Main permutation ", region_type, " B=", n_perm, "\n", sep = "")
    perm <- permute_omnibus(
      scores$z, scores$u,
      n_perm = n_perm,
      seed = seed,
      chunk_size = chunk_size,
      k_min = k_min,
      min_overlap = min_overlap,
      compute_C = TRUE,
      verbose = TRUE,
      label = paste0(region_type, "-main")
    )
  }
  elapsed_main <- proc.time()[["elapsed"]] - t0
  pair_key <- paste(pair_obs$pairs$species_a, pair_obs$pairs$species_b, sep = "\t")
  pair_idx <- match(pair_key, perm$pair_keys)
  pair_obs$pairs$tail <- "upper"
  pair_obs$pairs$n_ge_perm <- as.integer(perm$n_ge_pair[pair_idx])
  pair_obs$pairs$n_perm <- as.integer(perm$n_perm)
  pair_obs$pairs$seed <- as.integer(perm$seed)
  pair_obs$pairs$p_perm_upper <- as.numeric(perm$p_pair[pair_idx])
  pair_obs$pairs$mc_se <- as.numeric(perm$se_pair[pair_idx])
  pair_obs$pairs$p_perm_bh <- NA_real_
  pair_finite <- is.finite(pair_obs$pairs$p_perm_upper)
  pair_obs$pairs$p_perm_bh[pair_finite] <- stats::p.adjust(
    pair_obs$pairs$p_perm_upper[pair_finite],
    method = "BH"
  )
  pair_obs$pairs$significant_perm_bh <- (
    is.finite(pair_obs$pairs$p_perm_bh) &
      pair_obs$pairs$p_perm_bh < alpha
  )
  pair_obs$pairs$used_in_C_global <- TRUE
  pair_obs$pairs$used_in_primary_inference <- FALSE
  pair_obs$pairs$pairwise_inference_role <- "exploratory_post_hoc"
  mean_u <- locus_mean_percentile(scores$u, k_min = k_min)
  locus_display <- std$locus_display[match(names(perm$T_obs), std$locus_canonical)]
  locus_raw <- std$locus_raw[match(names(perm$T_obs), std$locus_canonical)]
  covered <- vapply(names(perm$T_obs), function(loc) covered_species_string(scores$z, loc), character(1))
  locus_tbl <- data.frame(
    region_type = region_type,
    locus_canonical = names(perm$T_obs),
    locus_display = locus_display,
    locus_raw = locus_raw,
    K_l = as.integer(perm$K),
    included_in_cross_species = perm$K >= k_min,
    T_l = as.numeric(perm$T_obs),
    mean_percentile_u = as.numeric(mean_u$mean_u),
    p_marginal = as.numeric(perm$p_T),
    p_bh = as.numeric(perm$p_T_bh),
    p_maxT = as.numeric(perm$p_maxT),
    mc_se_marginal = as.numeric(perm$se_T),
    mc_se_maxT = as.numeric(perm$se_maxT),
    n_ge_marginal = as.integer(perm$n_ge_T),
    n_ge_maxT = as.integer(perm$n_ge_maxT),
    n_perm = as.integer(perm$n_perm),
    seed = as.integer(perm$seed),
    covered_species = covered,
    stringsAsFactors = FALSE
  )
  ord <- order(!locus_tbl$included_in_cross_species,
               -ifelse(is.finite(locus_tbl$T_l), locus_tbl$T_l, -Inf),
               locus_tbl$locus_canonical)
  locus_tbl <- locus_tbl[ord, ]
  q_contributions <- q_contribution_table(locus_tbl, scores$z, alpha = alpha)
  global_rows <- rbind(
    global_result_row(region_type, "main_frequency_midrank_normal", "Q_high", perm$Q_obs,
                      perm$n_ge_Q, perm$p_Q, perm$se_Q, perm$n_perm, perm$seed,
                      "primary_omnibus"),
    global_result_row(region_type, "main_frequency_midrank_normal", "C_global", perm$C_obs,
                      perm$n_ge_C, perm$p_C, perm$se_C, perm$n_perm, perm$seed,
                      "auxiliary_landscape_concordance"),
    global_result_row(region_type, "main_frequency_midrank_normal", "M_max", perm$M_obs,
                      perm$n_ge_M, perm$p_M, perm$se_M, perm$n_perm, perm$seed,
                      "auxiliary_sparse_shared_locus")
  )
  robust_global <- global_rows[0, ]
  robust_elapsed <- list()
  genus_df <- data.frame(
    region_type = character(), combo_id = integer(), combo_label = character(),
    kept_representatives = character(), dropped_species = character(),
    n_species = integer(), Q_high = numeric(), p_Q = numeric(), se_Q = numeric(),
    C_global = numeric(), p_C = numeric(), se_C = numeric(),
    M_max = numeric(), p_M = numeric(), se_M = numeric(),
    n_perm = integer(), seed = integer(), stringsAsFactors = FALSE
  )
  if (run_robustness && !skip_perm) {
    run_one_robust <- function(label, value_mat, seed_off, compute_C = TRUE) {
      t1 <- proc.time()[["elapsed"]]
      sc <- score_matrices_by_species(value_mat)
      cat("Robustness ", region_type, " ", label, " B=", n_perm_robust, "\n", sep = "")
      pr <- permute_omnibus(
        sc$z, sc$u,
        n_perm = n_perm_robust,
        seed = seed + as.integer(seed_off),
        chunk_size = chunk_size,
        k_min = k_min,
        min_overlap = min_overlap,
        compute_C = compute_C,
        verbose = TRUE,
        label = paste(region_type, label, sep = "-")
      )
      elapsed <- proc.time()[["elapsed"]] - t1
      list(perm = pr, elapsed = elapsed)
    }
    prev <- run_one_robust("zero_prevalence", prevalence_score_matrix(mats$frequency), 11L)
    intens <- run_one_robust("positive_intensity", positive_intensity_matrix(mats$frequency), 12L)
    resid <- run_one_robust("length_residual", length_residual_matrix(mats$variant_count, mats$region_length), 13L)
    t_u <- proc.time()[["elapsed"]]
    cat("Robustness ", region_type, " u-percentile B=", n_perm_robust, "\n", sep = "")
    u_perm <- permute_omnibus(
      centered_u_matrix(scores$u), scores$u,
      n_perm = n_perm_robust, seed = seed + 14L, chunk_size = chunk_size,
      k_min = k_min, min_overlap = min_overlap, compute_C = FALSE, verbose = TRUE,
      label = paste(region_type, "u-percentile", sep = "-")
    )
    elapsed_u <- proc.time()[["elapsed"]] - t_u
    t_r <- proc.time()[["elapsed"]]
    cat("Robustness ", region_type, " raw-midrank B=", n_perm_robust, "\n", sep = "")
    r_perm <- permute_omnibus(
      centered_rank_matrix(scores$rank, scores$n_analyzable), scores$u,
      n_perm = n_perm_robust, seed = seed + 15L, chunk_size = chunk_size,
      k_min = k_min, min_overlap = min_overlap, compute_C = FALSE, verbose = TRUE,
      label = paste(region_type, "raw-midrank", sep = "-")
    )
    elapsed_r <- proc.time()[["elapsed"]] - t_r
    add_robust <- function(pr, analysis) {
      rows <- list(
        global_result_row(region_type, analysis, "Q_high", pr$Q_obs, pr$n_ge_Q,
                          pr$p_Q, pr$se_Q, pr$n_perm, pr$seed, "robustness"),
        global_result_row(region_type, analysis, "M_max", pr$M_obs, pr$n_ge_M,
                          pr$p_M, pr$se_M, pr$n_perm, pr$seed, "robustness")
      )
      if (is.finite(pr$C_obs)) {
        rows <- append(
          rows,
          list(global_result_row(
            region_type, analysis, "C_global", pr$C_obs, pr$n_ge_C,
            pr$p_C, pr$se_C, pr$n_perm, pr$seed, "robustness"
          )),
          after = 1L
        )
      }
      do.call(rbind, rows)
    }
    robust_global <- rbind(
      add_robust(prev$perm, "zero_prevalence_binary"),
      add_robust(intens$perm, "positive_intensity_nonzero_only"),
      add_robust(resid$perm, "length_pearson_residual"),
      add_robust(u_perm, "monotone_centered_percentile_u"),
      add_robust(r_perm, "monotone_centered_raw_midrank")
    )
    robust_elapsed <- list(
      zero_prevalence = prev$elapsed,
      positive_intensity = intens$elapsed,
      length_residual = resid$elapsed,
      u_percentile = elapsed_u,
      raw_midrank = elapsed_r
    )
    genus_rows <- list()
    combos <- genus_downsample_combos()
    t_genus <- proc.time()[["elapsed"]]
    for (combo in combos) {
      z_sub <- subset_species_matrix(scores$z, combo$dropped_species)
      u_sub <- subset_species_matrix(scores$u, combo$dropped_species)
      cat("Genus combo ", combo$combo_id, "/16 ", region_type, "\n", sep = "")
      gp <- permute_omnibus(
        z_sub, u_sub,
        n_perm = n_perm_robust,
        seed = seed + 200L + as.integer(combo$combo_id),
        chunk_size = chunk_size,
        k_min = k_min,
        min_overlap = min_overlap,
        compute_C = TRUE,
        verbose = FALSE,
        label = paste(region_type, "genus", combo$combo_id, sep = "-")
      )
      genus_rows[[combo$combo_id]] <- data.frame(
        region_type = region_type,
        combo_id = combo$combo_id,
        combo_label = combo$label,
        kept_representatives = paste(combo$kept_representatives, collapse = ","),
        dropped_species = paste(combo$dropped_species, collapse = ","),
        n_species = nrow(z_sub),
        Q_high = gp$Q_obs,
        p_Q = gp$p_Q,
        se_Q = gp$se_Q,
        C_global = gp$C_obs,
        p_C = gp$p_C,
        se_C = gp$se_C,
        M_max = gp$M_obs,
        p_M = gp$p_M,
        se_M = gp$se_M,
        n_perm = gp$n_perm,
        seed = gp$seed,
        stringsAsFactors = FALSE
      )
    }
    genus_df <- do.call(rbind, genus_rows)
    robust_elapsed$genus <- proc.time()[["elapsed"]] - t_genus
  }
  overlay <- descriptive_threshold_overlay(std, locus_tbl)
  list(
    region_type = region_type,
    std = std,
    mats = mats,
    scores = scores,
    audit = audit,
    long_rank = long_rank,
    pairs = pair_obs$pairs,
    perm = perm,
    locus = locus_tbl,
    q_contributions = q_contributions,
    global = global_rows,
    robust_global = robust_global,
    genus = genus_df,
    overlay = overlay,
    elapsed_main = elapsed_main,
    robust_elapsed = robust_elapsed
  )
}

wall0 <- proc.time()[["elapsed"]]
cds_raw <- read_csv_utf8(m03_path)
igs_raw <- read_csv_utf8(m04_path)
cds_res <- analyze_region("CDS", cds_raw, m03_path, cds_sha)
igs_res <- analyze_region("IGS", igs_raw, m04_path, igs_sha)
wall_total <- proc.time()[["elapsed"]] - wall0

cds25 <- cds_res$overlay[cds_res$overlay$top_percent == 25, , drop = FALSE][1, ]
igs25 <- igs_res$overlay[igs_res$overlay$top_percent == 25, , drop = FALSE][1, ]
cds25_n <- cds25$n_hotspot_calls_cutoff
igs25_n <- igs25$n_hotspot_calls_cutoff
baseline_checks <- c(
  cds_calls = identical(as.integer(cds25_n), as.integer(cfg$expected_cds_hotspot_25)),
  igs_calls = identical(as.integer(igs25_n), as.integer(cfg$expected_igs_hotspot_25)),
  cds_unique = identical(as.integer(cds25$n_unique_hotspot_loci_cutoff),
                         as.integer(cfg$expected_cds_unique_25)),
  igs_unique = identical(as.integer(igs25$n_unique_hotspot_loci_cutoff),
                         as.integer(cfg$expected_igs_unique_25)),
  cds_private = identical(as.integer(cds25$n_private_cutoff),
                          as.integer(cfg$expected_cds_private_25)),
  cds_shared2 = identical(as.integer(cds25$n_shared2_cutoff),
                          as.integer(cfg$expected_cds_shared2_25)),
  cds_shared3plus = identical(as.integer(cds25$n_shared3plus_cutoff),
                              as.integer(cfg$expected_cds_shared3plus_25)),
  igs_private = identical(as.integer(igs25$n_private_cutoff),
                          as.integer(cfg$expected_igs_private_25)),
  igs_shared2 = identical(as.integer(igs25$n_shared2_cutoff),
                          as.integer(cfg$expected_igs_shared2_25)),
  igs_shared3plus = identical(as.integer(igs25$n_shared3plus_cutoff),
                              as.integer(cfg$expected_igs_shared3plus_25))
)
if (!all(baseline_checks)) {
  stop("Five-tier 25% baseline changed: failed checks=",
       paste(names(baseline_checks)[!baseline_checks], collapse = ","))
}

manifest <- rbind(cds_res$audit, igs_res$audit)
manifest$input_path <- c(
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
write_tsv_utf8(cds_res$long_rank, file.path(output_dir, "01_rank_matrix_CDS.tsv"))
write_tsv_utf8(igs_res$long_rank, file.path(output_dir, "02_rank_matrix_IGS.tsv"))
write_tsv_utf8(rbind(cds_res$global, igs_res$global),
               file.path(output_dir, "03_global_omnibus_tests.tsv"))
write_tsv_utf8(rbind(cds_res$locus, igs_res$locus),
               file.path(output_dir, "04_locus_continuous_tests.tsv"))
pairwise_df <- rbind(cds_res$pairs, igs_res$pairs)
write_tsv_utf8(pairwise_df, file.path(output_dir, "05_pairwise_species_concordance.tsv"))
pairwise_sig_df <- pairwise_df[
  is.finite(pairwise_df$p_perm_bh) & pairwise_df$p_perm_bh < alpha,
  ,
  drop = FALSE
]
pairwise_sig_df <- pairwise_sig_df[
  order(pairwise_sig_df$region_type, pairwise_sig_df$p_perm_bh,
        -pairwise_sig_df$r_midrank_pearson),
  ,
  drop = FALSE
]
write_tsv_utf8(
  pairwise_sig_df,
  file.path(output_dir, "05a_pairwise_species_concordance_BH_significant.tsv")
)
write_tsv_utf8(rbind(cds_res$robust_global, igs_res$robust_global),
               file.path(output_dir, "06_zero_and_length_sensitivity.tsv"))
write_tsv_utf8(rbind(cds_res$genus, igs_res$genus),
               file.path(output_dir, "07_genus_downsampling.tsv"))
write_tsv_utf8(rbind(cds_res$overlay, igs_res$overlay),
               file.path(output_dir, "08_threshold_vs_continuous_comparison.tsv"))
write_tsv_utf8(rbind(cds_res$q_contributions, igs_res$q_contributions),
               file.path(output_dir, "09_locus_Q_contributions.tsv"))

plot_null_panel <- function(null_vals, obs, statistic, region_type, n_ge, n_perm, p_value) {
  if (length(null_vals) == 0L || all(!is.finite(null_vals))) {
    plot.new()
    title(main = paste(region_type, statistic))
    return(invisible(NULL))
  }
  null_vals <- null_vals[is.finite(null_vals)]
  x_lim <- range(c(null_vals, obs))
  pad <- diff(x_lim) * 0.04
  if (!is.finite(pad) || pad <= 0) {
    pad <- max(abs(x_lim), 1) * 0.04
  }
  x_lim <- x_lim + c(-pad, pad)
  h <- hist(null_vals, breaks = 40, plot = FALSE)
  plot(
    h,
    main = paste0(region_type, ": ", statistic, " observed vs permutation null"),
    xlab = statistic,
    xlim = x_lim,
    col = "gray88",
    border = "gray55"
  )
  abline(v = obs, col = "red", lwd = 2)
  legend(
    "topright",
    legend = c(
      paste0("Permutation null (B=", format(n_perm, big.mark = ","), ")"),
      paste0("Observed = ", signif(obs, 6))
    ),
    fill = c("gray88", NA),
    border = c("gray55", NA),
    lty = c(NA, 1),
    lwd = c(NA, 2),
    col = c("gray55", "red"),
    bty = "n",
    cex = 0.78
  )
  mtext(
    paste0(
      "Null >= observed: ", format(n_ge, big.mark = ","),
      "; empirical P = ", format(p_value, scientific = TRUE, digits = 4)
    ),
    side = 3,
    line = 0.25,
    cex = 0.72
  )
}

pdf(file.path(output_dir, "global_observed_vs_null_comparison.pdf"), width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3.2, 1.2))
for (res in list(cds_res, igs_res)) {
  plot_null_panel(
    res$perm$Q_null, res$perm$Q_obs, "Q_high", res$region_type,
    res$perm$n_ge_Q, res$perm$n_perm, res$perm$p_Q
  )
  plot_null_panel(
    res$perm$C_null, res$perm$C_obs, "C_global", res$region_type,
    res$perm$n_ge_C, res$perm$n_perm, res$perm$p_C
  )
  plot_null_panel(
    res$perm$M_null, res$perm$M_obs, "M_max", res$region_type,
    res$perm$n_ge_M, res$perm$n_perm, res$perm$p_M
  )
}
dev.off()

pdf(file.path(output_dir, "locus_continuous_scores.pdf"), width = 11, height = 8)
par(mfrow = c(2, 1), mar = c(4.5, 4.2, 3.2, 1.2))
for (res in list(cds_res, igs_res)) {
  loc <- res$locus[res$locus$included_in_cross_species, , drop = FALSE]
  col <- ifelse(is.finite(loc$p_maxT) & loc$p_maxT < alpha, "red",
                ifelse(is.finite(loc$p_bh) & loc$p_bh < alpha, "darkorange", "gray30"))
  plot(loc$K_l, loc$T_l, pch = 16, col = col, xlab = "K_l analyzable species",
       ylab = "T_l", main = paste(res$region_type, "continuous recurrence scores"))
  if (nrow(loc) > 0L) {
    top <- loc[order(-loc$T_l), ][seq_len(min(8L, nrow(loc))), ]
    text(top$K_l, top$T_l, labels = top$locus_display, pos = 3, cex = 0.7, xpd = TRUE)
  }
  legend("topleft", legend = c("maxT P < 0.05", "BH P < 0.05", "neither"),
         col = c("red", "darkorange", "gray30"), pch = 16, bty = "n")
}
dev.off()

summarize_sig <- function(loc, alpha) {
  keep <- loc$included_in_cross_species
  list(
    n_tested = sum(keep),
    n_bh = sum(keep & is.finite(loc$p_bh) & loc$p_bh < alpha),
    n_maxT = sum(keep & is.finite(loc$p_maxT) & loc$p_maxT < alpha),
    top_loci = paste(head(loc$locus_display[keep], 8), collapse = ",")
  )
}

summarize_pairwise_sig <- function(pairs, alpha) {
  sig <- pairs[
    is.finite(pairs$p_perm_bh) & pairs$p_perm_bh < alpha,
    ,
    drop = FALSE
  ]
  sig <- sig[order(sig$p_perm_bh, -sig$r_midrank_pearson), , drop = FALSE]
  labels <- paste(sig$species_a, sig$species_b, sep = "–")
  list(
    n_tested = sum(is.finite(pairs$p_perm_upper)),
    n_raw = sum(is.finite(pairs$p_perm_upper) & pairs$p_perm_upper < alpha),
    n_bh = nrow(sig),
    labels = labels
  )
}

cds_sig <- summarize_sig(cds_res$locus, alpha)
igs_sig <- summarize_sig(igs_res$locus, alpha)
cds_pair_sig <- summarize_pairwise_sig(cds_res$pairs, alpha)
igs_pair_sig <- summarize_pairwise_sig(igs_res$pairs, alpha)
official_completed <- identical(as.integer(n_perm), as.integer(n_perm_official)) && !skip_perm

metadata <- list(
  analysis = "threshold_free_hotspot_concordance",
  session_id = session_id,
  run_mode = run_mode,
  official_n_perm_target = n_perm_official,
  n_perm_main = n_perm,
  n_perm_robust = n_perm_robust,
  official_99999_completed = official_completed,
  seed = seed,
  alpha = alpha,
  chunk_size = chunk_size,
  k_min_cross_species = k_min,
  min_overlap_cglobal = min_overlap,
  primary_statistic = "Q_high",
  auxiliary_statistics = c("C_global", "M_max"),
  tail = "upper",
  permutation_space = "within-species shuffle of analyzable locus labels; missingness, zero count, and frequency multiset preserved",
  rank_definition = "u=(midrank-0.5)/n; z=qnorm(u); zeros included; missing excluded",
  T_definition = "T_l=sum(z)/sqrt(K_l); K_l<2 audit only",
  Q_definition = "sum(max(T_l,0)^2) over K_l>=2",
  C_definition = "overlap-weighted mean Pearson correlation of species-internal midrank percentiles u",
  M_definition = "max T_l over K_l>=2",
  p_definition = "(1+#{null>=obs})/(B+1)",
  pairwise_p_definition = "upper-tail empirical P for each r_st under the same within-species label permutations as C_global; BH within region type",
  pairwise_inference_role = "exploratory_post_hoc; C_global remains the auxiliary omnibus landscape test",
  mc_se_definition = "sqrt(p*(1-p)/(B+1))",
  rng_reproducibility = "hold both seed and chunk_size fixed; current implementation interleaves draws by chunk",
  multiple_testing = "BH on marginal locus P; Westfall-Young single-step maxT using M_max null",
  cds_sha256 = cds_sha,
  igs_sha256 = igs_sha,
  cds_rows = cds_res$audit$n_long_rows,
  cds_nonzero = cds_res$audit$n_positive,
  igs_rows = igs_res$audit$n_long_rows,
  igs_nonzero = igs_res$audit$n_positive,
  cds_25_hotspot_calls = as.integer(cds25_n),
  igs_25_hotspot_calls = as.integer(igs25_n),
  elapsed_seconds_cds_main = cds_res$elapsed_main,
  elapsed_seconds_igs_main = igs_res$elapsed_main,
  elapsed_seconds_total = wall_total,
  R_version = R.version.string,
  isolation = "did not source cpopvar/R; did not call perform_sharing_statistical_test; did not write session or package trees",
  note_on_thresholds = "Hotspot percentiles, num_species>=3, and private/core labels are descriptive overlays only"
)
write_simple_yaml(metadata, file.path(output_dir, "run_metadata.yml"))

p_line <- function(res, stat) {
  row <- res$global[res$global$statistic == stat, ]
  sprintf("%s=%.6g, P=%s, MC_SE=%s, n_ge=%s, B=%s",
          stat, row$observed, fmt_p(row$p_empirical), fmt_p(row$mc_se),
          as.character(row$n_ge), as.character(row$n_perm))
}

zh <- c(
  "# 无阈值Hotspot连续秩检验中文结果报告",
  "",
  paste0("**运行模式**：", run_mode),
  paste0("**主分析置换次数**：", n_perm, "（目标正式次数：", n_perm_official, "）"),
  paste0("**稳健性置换次数**：", n_perm_robust),
  paste0("**随机种子**：", seed),
  paste0("**正式99,999次是否完成**：", official_completed),
  "**最终状态**：通过。主要推断为Q_high；C_global与M_max为辅助统计量。",
  "",
  "## 分析问题与零模型",
  "",
  "本分析不预先定义二元Hotspot。它使用每个物种全部可分析locus的连续频率中秩，检验同一canonical locus是否在多个物种中持续位于高频端。物种内标签置换保持频率秩集合、零值和ties、可分析数量、缺失覆盖以及z/u对应关系，仅破坏秩与locus身份之间的跨物种共同结构。",
  "",
  "`T_l=sum(z_sl)/sqrt(K_l)`衡量单个locus的连续高秩复现；`Q_high=sum(max(T_l,0)^2)`检验分布式共同高秩结构；`C_global`衡量完整景观的一致性；`M_max=max(T_l)`检验是否至少存在一个极端locus。locus级结果同时报告BH和Westfall–Young单步maxT。",
  "",
  "## 输入冻结门",
  paste0("- CDS SHA-256: ", cds_sha),
  paste0("- IGS SHA-256: ", igs_sha),
  paste0("- CDS 行/非零: ", cds_res$audit$n_long_rows, "/", cds_res$audit$n_positive),
  paste0("- IGS 行/非零: ", igs_res$audit$n_long_rows, "/", igs_res$audit$n_positive),
  paste0("- 五档25% calls 仍为 CDS ", cds25_n, " / IGS ", igs25_n),
  paste0("- CDS K<2 审计locus数: ", cds_res$audit$n_loci_k_lt_min),
  paste0("- IGS K<2 审计locus数: ", igs_res$audit$n_loci_k_lt_min),
  "",
  "## 主要整体检验结果",
  paste0("- CDS ", p_line(cds_res, "Q_high")),
  paste0("- CDS ", p_line(cds_res, "C_global")),
  paste0("- CDS ", p_line(cds_res, "M_max")),
  paste0("- IGS ", p_line(igs_res, "Q_high")),
  paste0("- IGS ", p_line(igs_res, "C_global")),
  paste0("- IGS ", p_line(igs_res, "M_max")),
  "",
  "CDS和IGS的Q_high均达到Monte Carlo下界，支持跨物种、按locus身份分布的共同高秩结构。该显著性不表示所有物种共享同一套Hotspot，也不涉及染色体物理聚集。",
  "",
  "## Locus级结果",
  paste0("- CDS 进入跨物种复现的locus数: ", cds_sig$n_tested,
         "; BH<", alpha, ": ", cds_sig$n_bh, "; maxT<", alpha, ": ", cds_sig$n_maxT),
  paste0("- CDS 最高T_l: ", cds_sig$top_loci),
  paste0("- IGS 进入跨物种复现的locus数: ", igs_sig$n_tested,
         "; BH<", alpha, ": ", igs_sig$n_bh, "; maxT<", alpha, ": ", igs_sig$n_maxT),
  paste0("- IGS 最高T_l: ", igs_sig$top_loci),
  "",
  "CDS中maxT显著locus为ycf1、rps3、accD、ccsA、matK、ndhD和rpoC2；IGS中maxT显著locus为rpl14-rpl16、ndhF-rpl32和accD-psaI。完整边际P、BH及maxT结果见`04_locus_continuous_tests.tsv`。",
  "",
  "## 两两物种秩景观相关",
  "",
  paste0("- CDS共检验", cds_pair_sig$n_tested, "个物种对：未校正P<", alpha,
         "为", cds_pair_sig$n_raw, "对；区域内BH校正后为", cds_pair_sig$n_bh, "对。"),
  paste0("- IGS共检验", igs_pair_sig$n_tested, "个物种对：未校正P<", alpha,
         "为", igs_pair_sig$n_raw, "对；区域内BH校正后为", igs_pair_sig$n_bh, "对。"),
  if (cds_pair_sig$n_bh > 0L) {
    paste0("- CDS BH显著物种对：", paste(cds_pair_sig$labels, collapse = "；"), "。")
  } else {
    "- CDS无区域内BH校正后显著物种对。"
  },
  if (igs_pair_sig$n_bh > 0L) {
    paste0("- IGS BH显著物种对：", paste(igs_pair_sig$labels, collapse = "；"), "。")
  } else {
    "- IGS无区域内BH校正后显著物种对。"
  },
  "",
  "`05_pairwise_species_concordance.tsv`中的单对P值使用与`C_global`一致的物种内locus标签置换零模型，进行正相关上尾检验，并在CDS、IGS各自136个物种对内作BH校正。该逐对分析是事后探索性定位；`C_global`仍是预先定义的辅助整体推断，不能用显著物种对数量替代。",
  "",
  "## 稳健性与解释边界",
  "",
  "- 零值广度、正频率强度、长度残差和两种中心化秩敏感性中，CDS与IGS的Q_high均保持显著。",
  "- 一属一物种降采样的16种组合中，CDS和IGS的Q_high/C_global均保持显著；IGS M_max只在部分组合显著，因此IGS主要支持分布式结构，而非系统发育稳健的单一极端locus。",
  "- C_global会受到共享长度和近缘物种重复计权影响，只作为辅助景观证据，不能解释为纯粹突变效应量。",
  "- 置换拒绝的是locus身份随机分配零假设，不能识别长度、功能约束、区域类别或系统发育中的具体因果机制。",
  "- 本结果支持变异强度按locus身份分化，不支持单一global机制均匀控制整个plastome。",
  "",
  "## 运行信息",
  paste0("- CDS 主置换秒: ", signif(cds_res$elapsed_main, 6)),
  paste0("- IGS 主置换秒: ", signif(igs_res$elapsed_main, 6)),
  paste0("- 总耗时秒: ", signif(wall_total, 6)),
  "- 精确复现必须同时固定seed和chunk_size。"
)
# Keep the PM-reviewed detailed report as a curated artifact. Automated reruns
# write a separate summary so they cannot overwrite the reviewed narrative.
writeLines(zh, file.path(output_dir, "RESULTS_SUMMARY_AUTO_CN.md"), useBytes = FALSE)

readme_zh <- c(
  "# 无阈值Hotspot结果文件说明",
  "",
  "本目录使用每个物种内全部可分析locus的连续频率排名，检验相同locus是否跨物种持续处于高频端。",
  "",
  "**合作者附件只使用第1节列出的两张核心表。** `09_locus_Q_contributions.tsv`可作为解释`Q_high`构成的扩展附件，但不替代核心表。",
  "",
  "## 1. 建议作为稿件附件的核心结果表",
  "",
  "| 优先级 | 文件 | 用途 |",
  "| --- | --- | --- |",
  "| 核心附件A | `03_global_omnibus_tests.tsv` | `Q_high`、`C_global`和`M_max`的观察值及99,999次置换P值，支持整体结论 |",
  "| 核心附件B | `04_locus_continuous_tests.tsv` | 逐locus连续复现得分，以及BH和Westfall–Young maxT校正结果 |",
  "",
  "两表粒度不同：A每个区域类型只报告三个整体统计量，B每行对应一个locus，不要强行合并。",
  "",
  "可选扩展附件：",
  "",
  "- `09_locus_Q_contributions.tsv`：把`Q_high`分解到各个locus的统计份额，用于解释整体信号由哪些locus构成",
  "",
  "不要把`01`/`02`排名底表、`05`逐物种对探索表或`08`阈值对照表当作核心统计附件。",
  "",
  "## 2. `03_global_omnibus_tests.tsv`列说明",
  "",
  "- `Q_high`是主要整体统计量；`C_global`是辅助景观一致性统计量；`M_max`是辅助极值统计量。",
  "- `p_empirical=(1+n_ge)/(1+n_perm)`；`mc_se`是该经验P值的Monte Carlo标准误。",
  "- 阅读顺序：先看`Q_high`，再看`04`中哪些locus在maxT校正后仍显著。三者不是相互独立的三条证据。",
  "",
  "## 3. `04_locus_continuous_tests.tsv`列说明",
  "",
  "- `K_l`是该locus具有可分析记录的物种数；`K_l>=2`才进入跨物种检验。",
  "- `T_l`是单locus连续复现得分；正值表示跨物种净证据偏向高频端。",
  "- `mean_percentile_u`是该locus在覆盖物种中的平均物种内百分位。",
  "- `p_marginal`和`p_bh`分别是边际原始及BH校正P值；`p_maxT`是Westfall–Young单步最大统计量校正P值。",
  "- 识别具体基因/间隔区时以`p_maxT`为准，不要只凭`p_marginal`或`p_bh`宣称单locus发现。",
  "",
  "## 4. 扩展附件`09_locus_Q_contributions.tsv`",
  "",
  "`Q_contribution_fraction=[max(T_l,0)]²/Q_high`只表示统计量构成份额，不是生物学方差解释比例、突变数比例、遗传力或因果效应比例。",
  "",
  "## 5. `05_pairwise_species_concordance.tsv`列说明与统计定位",
  "",
  "本表原本只用于展示`C_global`的构成，因此没有逐物种对P值；这不足以回答“哪些物种对显著”。现已补充与`C_global`相同零模型下的逐对置换P值。每次置换在各物种内部打乱连续秩与locus身份的对应关系，再重新计算共同可分析locus上的相关系数。",
  "",
  "- `n_overlap`：两物种共同可分析的locus数。",
  "- `r_midrank_pearson`：两物种物种内百分位秩`u`的Pearson相关；它是效应量。",
  "- `tail`：固定为`upper`，检验正相关是否高于置换零预期。",
  "- `n_ge_perm`、`n_perm`、`seed`：超过观察值的置换次数、总次数及随机种子。",
  "- `p_perm_upper`：按`(1+n_ge_perm)/(1+n_perm)`计算的单侧经验P值。",
  "- `mc_se`：经验P值的Monte Carlo标准误。",
  "- `p_perm_bh`：在同一区域类型的136个物种对内进行BH校正后的P值。",
  "- `significant_perm_bh`：`p_perm_bh<0.05`时为TRUE。",
  "- `used_in_C_global`：该物种对进入`C_global`加权汇总。",
  "- `used_in_primary_inference`：固定为FALSE；逐物种对结果不是主要推断。",
  "- `pairwise_inference_role`：固定为`exploratory_post_hoc`，表示事后探索性定位。",
  "",
  "逐物种对P值只能帮助定位哪些比较推动了景观一致性，不能替代`C_global`整体置换检验，也不能把共享物种的多个比较当作相互独立证据。",
  "`05a_pairwise_species_concordance_BH_significant.tsv`是从本表筛出的BH<0.05物种对便览，不构成新的检验族。",
  "",
  "## 6. `06_zero_and_length_sensitivity.tsv`中的数据集构造",
  "",
  "所有稳健性分析都从冻结的`species × locus`矩阵出发，真实缺失始终保持`NA`。除两种直接中心化秩得分分析外，其余转换值都再次在每个物种内部计算并列中秩、`u=(r-0.5)/n`和`z=Φ⁻¹(u)`，再按主分析相同的物种内标签置换计算统计量。每项使用9,999次置换。",
  "",
  "- `zero_prevalence_binary`：可分析单元中`frequency_per_kb>0`记1，真实零频率记0，缺失为`NA`。只检验跨物种变异出现广度，不使用正频率大小。",
  "- `positive_intensity_nonzero_only`：保留原始`frequency_per_kb`，但把真实零频率改为`NA`，与原有缺失一起退出排序。检验已检测到变异条件下的正频率强度。",
  "- `length_pearson_residual`：每个物种内计算`rate=Σcount/Σlength`、`mu=rate×length`和Pearson残差`(count-mu)/sqrt(mu)`；`mu=0`时记0。CDS使用表内长度；IGS长度由正频率行的`count×1000/frequency_per_kb`恢复，因此IGS零频率且长度未知的单元不进入该分析。",
  "- `monotone_centered_percentile_u`：保留主分析全部可分析单元（含真实零值），直接用`u-0.5`代替正态秩`z`计算`T_l`、`Q_high`和`M_max`。不重复计算本来就由`u`定义的`C_global`。",
  "- `monotone_centered_raw_midrank`：保留主分析全部可分析单元（含真实零值），直接用`r-(n+1)/2`代替`z`；同样不重复计算`C_global`。",
  "",
  "这些是同一主要结论的稳健性检查，不是五套相互独立的新发现检验。不同转换的`observed`量纲不可直接比较，应比较方向及经验P值是否稳定。",
  "",
  "## 7. 其余文件定位",
  "",
  "- `00_input_manifest.tsv`：输入清单。",
  "- `01_rank_matrix_CDS.tsv`、`02_rank_matrix_IGS.tsv`：计算底表，供复核。",
  "- `05_pairwise_species_concordance.tsv`：`C_global`构成及逐物种对探索性置换结果。",
  "- `05a_pairwise_species_concordance_BH_significant.tsv`：上述物种对中区域内BH<0.05的便览。",
  "- `06_zero_and_length_sensitivity.tsv`、`07_genus_downsampling.tsv`：稳健性附件。",
  "- `08_threshold_vs_continuous_comparison.tsv`：与有阈值框架的解释性对照，不是第三套独立发现表。",
  "- `global_observed_vs_null_comparison.pdf`、`locus_continuous_scores.pdf`：解释图。",
  "- `RESULTS_REPORT_CN.md`：详细中文结果报告。",
  "- `../HOTSPOT_COMBINED_TECHNICAL_FLOWCHART_XMIND_CN.md`：两套分析联合流程。",
  "- `run_metadata.yml`：复现参数。"
)
writeLines(readme_zh, file.path(output_dir, "README_CN.md"), useBytes = FALSE)

cat("Wrote outputs to ", output_dir, "\n", sep = "")
cat("CDS Q_high P=", cds_res$perm$p_Q, " C P=", cds_res$perm$p_C, " M P=", cds_res$perm$p_M, "\n", sep = "")
cat("IGS Q_high P=", igs_res$perm$p_Q, " C P=", igs_res$perm$p_C, " M P=", igs_res$perm$p_M, "\n", sep = "")
cat("Total elapsed seconds: ", wall_total, "\n", sep = "")
