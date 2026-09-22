# Manuscript-facing hotspot threshold and sharing-test helpers.
# Pure statistics/utilities: no cpopvar task_params, no legacy driver phyper.

HOTSPOT_CUTOFF_MAP <- data.frame(
  top_percent = c(10, 15, 20, 25, 30),
  quantile_prob = c(0.90, 0.85, 0.80, 0.75, 0.70),
  stringsAsFactors = FALSE
)

MANUSCRIPT_BASELINE_25 <- list(
  CDS = list(n_calls = 109L, n_unique = 50L, private = 24L, shared2 = 6L, shared3plus = 20L),
  IGS = list(n_calls = 183L, n_unique = 63L, private = 21L, shared2 = 10L, shared3plus = 32L)
)

SESSION_COMPUTED_25 <- list(
  CDS = list(n_calls = 109L, n_unique = 50L, private = 24L, shared2 = 6L, shared3plus = 20L),
  IGS = list(n_calls = 183L, n_unique = 63L, private = 21L, shared2 = 10L, shared3plus = 32L)
)

as_logical_flag <- function(x) {
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES")
}

quantile_prob_for_top_percent <- function(top_percent) {
  key <- as.numeric(top_percent)
  hit <- HOTSPOT_CUTOFF_MAP$quantile_prob[HOTSPOT_CUTOFF_MAP$top_percent == key]
  if (length(hit) != 1L) {
    stop("Unsupported top_percent: ", paste(top_percent, collapse = ","),
         ". Allowed: ", paste(HOTSPOT_CUTOFF_MAP$top_percent, collapse = ", "))
  }
  as.numeric(hit)
}

file_sha256 <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  hashed <- tryCatch({
    out <- system2("sha256sum", shQuote(normalizePath(path, winslash = "/", mustWork = TRUE)),
                   stdout = TRUE, stderr = FALSE)
    sub("\\s.*$", "", out[[1]])
  }, error = function(e) NA_character_)
  if (length(hashed) != 1L || is.na(hashed) || !nzchar(hashed)) {
    paste0("md5:", unname(tools::md5sum(path)))
  } else {
    hashed
  }
}

file_manifest_row <- function(path, role, notes = "") {
  info <- file.info(path)
  data.frame(
    role = role,
    path = path,
    exists = file.exists(path),
    bytes = if (is.na(info$size)) NA_real_ else as.numeric(info$size),
    mtime_unix = if (is.na(info$mtime)) NA_real_ else as.numeric(info$mtime),
    sha256 = file_sha256(path),
    notes = notes,
    stringsAsFactors = FALSE
  )
}

strip_poigs_prefix <- function(x) {
  sub("^poiGS_", "", as.character(x))
}

species_key <- function(x) {
  gsub(" ", "_", trimws(as.character(x)))
}

read_csv_utf8 <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  fileEncoding = "UTF-8")
}

read_tsv_utf8 <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE,
                    fileEncoding = "UTF-8")
}

write_tsv_utf8 <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE,
                     na = "", fileEncoding = "UTF-8")
}

write_csv_utf8 <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

standardize_frequency_table <- function(df, region_type, locus_col,
                                        display_col = NULL,
                                        stored_threshold_col = "threshold",
                                        stored_flag_col = "is_potential_hotspot") {
  required <- c("species", locus_col, "frequency_per_kb")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing columns in ", region_type, " table: ", paste(missing, collapse = ", "))
  }
  display <- if (!is.null(display_col) && display_col %in% names(df)) {
    as.character(df[[display_col]])
  } else if (identical(locus_col, "poiGS_ID")) {
    strip_poigs_prefix(df[[locus_col]])
  } else {
    as.character(df[[locus_col]])
  }
  out <- data.frame(
    species = as.character(df$species),
    locus_id = as.character(df[[locus_col]]),
    locus_display = display,
    region_type = region_type,
    frequency_per_kb = as.numeric(df$frequency_per_kb),
    variant_count = if ("variant_count" %in% names(df)) as.numeric(df$variant_count) else NA_real_,
    stored_threshold = if (stored_threshold_col %in% names(df)) as.numeric(df[[stored_threshold_col]]) else NA_real_,
    stored_is_hotspot = if (stored_flag_col %in% names(df)) as_logical_flag(df[[stored_flag_col]]) else NA,
    stringsAsFactors = FALSE
  )
  if (anyNA(out$frequency_per_kb)) {
    stop(region_type, " frequency_per_kb contains NA")
  }
  out
}

identify_hotspots_nonzero_quantile <- function(df,
                                               quantile_prob,
                                               quantile_type = 7L,
                                               species_col = "species",
                                               freq_col = "frequency_per_kb") {
  if (!quantile_prob %in% HOTSPOT_CUTOFF_MAP$quantile_prob) {
    stop("quantile_prob must be one of: ",
         paste(HOTSPOT_CUTOFF_MAP$quantile_prob, collapse = ", "))
  }
  species <- as.character(df[[species_col]])
  freq <- as.numeric(df[[freq_col]])
  threshold <- rep(Inf, nrow(df))
  n_nonzero <- integer(nrow(df))
  n_total <- integer(nrow(df))
  for (sp in unique(species)) {
    idx <- which(species == sp)
    n_total[idx] <- length(idx)
    nz <- freq[idx] > 0
    n_nonzero[idx] <- sum(nz)
    if (!any(nz)) {
      next
    }
    thr <- as.numeric(stats::quantile(freq[idx][nz], probs = quantile_prob,
                                      type = quantile_type, names = FALSE, na.rm = TRUE))
    threshold[idx] <- thr
  }
  is_hotspot <- freq >= threshold
  hotspot_count_by_species <- tapply(is_hotspot, species, sum)
  actual_selected_fraction <- as.numeric(hotspot_count_by_species[species]) / n_nonzero
  actual_selected_fraction[!is.finite(actual_selected_fraction)] <- 0
  out <- df
  out$quantile_prob <- quantile_prob
  out$quantile_type <- as.integer(quantile_type)
  out$threshold <- threshold
  out$is_hotspot <- is_hotspot
  out$n_total_loci <- n_total
  out$n_nonzero <- n_nonzero
  out$n_hotspot_species <- as.integer(hotspot_count_by_species[species])
  out$actual_selected_fraction <- actual_selected_fraction
  out$nominal_top_fraction <- 1 - quantile_prob
  out
}

sharing_by_locus <- function(hotspot_df, locus_col = "locus_id") {
  hs <- hotspot_df[hotspot_df$is_hotspot, , drop = FALSE]
  if (nrow(hs) == 0) {
    return(data.frame(
      locus_id = character(), n_species = integer(), species_list = character(),
      stringsAsFactors = FALSE
    ))
  }
  split_sp <- split(as.character(hs$species), hs[[locus_col]])
  data.frame(
    locus_id = names(split_sp),
    n_species = vapply(split_sp, function(z) length(unique(z)), integer(1)),
    species_list = vapply(split_sp, function(z) paste(sort(unique(z)), collapse = ","), character(1)),
    stringsAsFactors = FALSE
  )
}

summarize_sharing_classes <- function(n_species) {
  n_species <- as.integer(n_species)
  list(
    n_unique = length(n_species),
    private = sum(n_species == 1L),
    shared2 = sum(n_species == 2L),
    shared3plus = sum(n_species >= 3L),
    max_shared_species = if (length(n_species)) max(n_species) else 0L
  )
}

hypergeometric_upper_tail <- function(x, K, N, n) {
  x <- as.numeric(x); K <- as.numeric(K); N <- as.numeric(N); n <- as.numeric(n)
  n_out <- max(length(x), length(K), length(N), length(n))
  x <- rep(x, length.out = n_out)
  K <- rep(K, length.out = n_out)
  N <- rep(N, length.out = n_out)
  n <- rep(n, length.out = n_out)
  p <- rep(NA_real_, n_out)
  for (i in seq_len(n_out)) {
    if (any(!is.finite(c(x[i], K[i], N[i], n[i])))) {
      next
    }
    if (N[i] <= 0 || n[i] < 0 || K[i] < 0 || n[i] > N[i] || K[i] > N[i]) {
      next
    }
    if (x[i] <= 0) {
      p[i] <- 1
      next
    }
    # P(X >= x) for X ~ Hypergeometric(N, K, n)
    p[i] <- stats::phyper(x[i] - 1, K[i], N[i] - K[i], n[i], lower.tail = FALSE)
  }
  p
}

bh_adjust_by_family <- function(p, family) {
  p <- as.numeric(p)
  family <- as.character(family)
  adj <- rep(NA_real_, length(p))
  for (fam in unique(family)) {
    idx <- which(family == fam)
    adj[idx] <- stats::p.adjust(p[idx], method = "BH")
  }
  adj
}

locus_test_table <- function(annotated_df, alpha = 0.05) {
  eligible_df <- annotated_df[
    is.finite(annotated_df$frequency_per_kb) & annotated_df$frequency_per_kb > 0,
    ,
    drop = FALSE
  ]
  N_nonzero <- nrow(eligible_df)
  n_hotspot_calls <- sum(eligible_df$is_hotspot)
  loci <- sort(unique(annotated_df$locus_id))
  rows <- lapply(loci, function(loc) {
    all_sub <- annotated_df[annotated_df$locus_id == loc, , drop = FALSE]
    sub <- eligible_df[eligible_df$locus_id == loc, , drop = FALSE]
    K_nonzero <- nrow(sub)
    x <- sum(sub$is_hotspot)
    species_called <- sort(unique(sub$species[sub$is_hotspot]))
    expected <- if (N_nonzero > 0) {
      n_hotspot_calls * K_nonzero / N_nonzero
    } else {
      NA_real_
    }
    data.frame(
      locus_id = loc,
      locus_display = unique(all_sub$locus_display)[[1]],
      region_type = unique(all_sub$region_type)[[1]],
      N_nonzero = N_nonzero,
      K_nonzero = K_nonzero,
      n_hotspot_calls = n_hotspot_calls,
      x = x,
      expected_x_nonzero = expected,
      enrichment_ratio_nonzero = if (is.finite(expected) && expected > 0) {
        x / expected
      } else {
        NA_real_
      },
      p_hyper_nonzero_raw = hypergeometric_upper_tail(
        x, K_nonzero, N_nonzero, n_hotspot_calls
      ),
      n_species = length(species_called),
      species_list = paste(species_called, collapse = ","),
      shared_ge3 = length(species_called) >= 3L,
      hypergeom_role = "coarse_sensitivity_only",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out
}

stratified_permutation_pvalues <- function(annotated_df,
                                           n_perm = 99999L,
                                           seed = 20251104L,
                                           eligible_nonzero = TRUE,
                                           return_null_distribution = FALSE) {
  species <- sort(unique(as.character(annotated_df$species)))
  loci <- sort(unique(as.character(annotated_df$locus_id)))
  obs_x <- setNames(integer(length(loci)), loci)
  tab <- tapply(annotated_df$is_hotspot, annotated_df$locus_id, sum)
  obs_x[names(tab)] <- as.integer(tab)

  eligible <- vector("list", length(species))
  names(eligible) <- species
  n_hot <- setNames(integer(length(species)), species)
  for (sp in species) {
    sub <- annotated_df[annotated_df$species == sp, , drop = FALSE]
    keep <- if (eligible_nonzero) sub$frequency_per_kb > 0 else rep(TRUE, nrow(sub))
    eligible[[sp]] <- as.character(sub$locus_id[keep])
    n_hot[[sp]] <- as.integer(sum(sub$is_hotspot))
    if (n_hot[[sp]] > length(eligible[[sp]])) {
      stop("Species ", sp, " has more hotspots than eligible loci")
    }
  }

  ge_count <- setNames(integer(length(loci)), loci)
  null_counts <- if (isTRUE(return_null_distribution)) {
    matrix(
      0L,
      nrow = length(species) + 1L,
      ncol = length(loci),
      dimnames = list(as.character(0:length(species)), loci)
    )
  } else {
    NULL
  }
  set.seed(as.integer(seed))
  for (b in seq_len(as.integer(n_perm))) {
    perm_x <- setNames(integer(length(loci)), loci)
    for (sp in species) {
      cand <- eligible[[sp]]
      nh <- n_hot[[sp]]
      if (nh <= 0L || length(cand) == 0L) {
        next
      }
      picked <- cand[sample.int(length(cand), nh, replace = FALSE)]
      perm_x[picked] <- perm_x[picked] + 1L
    }
    ge_count <- ge_count + as.integer(perm_x >= obs_x)
    if (isTRUE(return_null_distribution)) {
      idx <- cbind(as.integer(perm_x) + 1L, seq_along(loci))
      null_counts[idx] <- null_counts[idx] + 1L
    }
  }
  p_perm <- (ge_count + 1) / (as.integer(n_perm) + 1)
  pvalues <- data.frame(
    locus_id = loci,
    x_obs = as.integer(obs_x[loci]),
    p_perm_raw = as.numeric(p_perm[loci]),
    n_perm = as.integer(n_perm),
    perm_seed = as.integer(seed),
    perm_eligible = if (eligible_nonzero) "nonzero_frequency" else "all_analyzable_pairs",
    stringsAsFactors = FALSE
  )
  if (!isTRUE(return_null_distribution)) {
    return(pvalues)
  }
  null_rows <- lapply(seq_along(loci), function(j) {
    data.frame(
      locus_id = loci[j],
      x_null = 0:length(species),
      n_permutations = as.integer(null_counts[, j]),
      probability = as.numeric(null_counts[, j]) / as.integer(n_perm),
      stringsAsFactors = FALSE
    )
  })
  list(
    pvalues = pvalues,
    null_distribution = do.call(rbind, null_rows)
  )
}

overall_summary_row <- function(annotated_df, top_percent, quantile_prob) {
  sharing <- sharing_by_locus(annotated_df)
  classes <- summarize_sharing_classes(sharing$n_species)
  data.frame(
    region_type = unique(annotated_df$region_type)[[1]],
    top_percent = top_percent,
    quantile_prob = quantile_prob,
    n_analyzable_pairs = nrow(annotated_df),
    n_nonzero_pairs = sum(annotated_df$frequency_per_kb > 0),
    n_hotspot_calls = sum(annotated_df$is_hotspot),
    n_unique_loci = classes$n_unique,
    n_private = classes$private,
    n_shared2 = classes$shared2,
    n_shared3plus = classes$shared3plus,
    max_shared_species = classes$max_shared_species,
    stringsAsFactors = FALSE
  )
}

per_species_summary <- function(annotated_df, top_percent, quantile_prob) {
  species <- sort(unique(annotated_df$species))
  rows <- lapply(species, function(sp) {
    sub <- annotated_df[annotated_df$species == sp, , drop = FALSE]
    data.frame(
      region_type = unique(sub$region_type)[[1]],
      top_percent = top_percent,
      quantile_prob = quantile_prob,
      species = sp,
      n_total_loci = nrow(sub),
      n_nonzero = sum(sub$frequency_per_kb > 0),
      threshold = unique(sub$threshold)[[1]],
      n_hotspot_calls = sum(sub$is_hotspot),
      actual_selected_fraction = if (sum(sub$frequency_per_kb > 0) > 0)
        sum(sub$is_hotspot) / sum(sub$frequency_per_kb > 0) else 0,
      nominal_top_fraction = 1 - quantile_prob,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

count_match <- function(observed, expected) {
  identical(as.integer(observed), as.integer(expected))
}

reconciliation_counts <- function(annotated_25, stored_flag, label) {
  hs <- annotated_25[annotated_25$is_hotspot, , drop = FALSE]
  sharing <- sharing_by_locus(annotated_25)
  classes <- summarize_sharing_classes(sharing$n_species)
  stored_n <- if (all(is.na(stored_flag))) NA_integer_ else as.integer(sum(stored_flag, na.rm = TRUE))
  data.frame(
    source = label,
    region_type = unique(annotated_25$region_type)[[1]],
    n_calls = nrow(hs),
    n_unique = classes$n_unique,
    n_private = classes$private,
    n_shared2 = classes$shared2,
    n_shared3plus = classes$shared3plus,
    stored_flag_calls = stored_n,
    stringsAsFactors = FALSE
  )
}

try_write_xlsx <- function(sheets, path) {
  ok <- requireNamespace("openxlsx", quietly = TRUE)
  if (!ok) {
    return(list(ok = FALSE, path = NA_character_,
                note = "openxlsx is not installed; wrote TSV sheets instead"))
  }
  wb <- openxlsx::createWorkbook()
  for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, sheets[[nm]])
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  list(ok = TRUE, path = path, note = "Wrote xlsx with openxlsx")
}

read_simple_yaml_flat <- function(path) {
  if (!file.exists(path)) {
    return(list())
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  out <- list()
  for (raw in lines) {
    line <- sub("\\s+#.*$", "", raw)
    if (!nzchar(trimws(line)) || startsWith(trimws(line), "#")) {
      next
    }
    if (startsWith(line, " ") || startsWith(line, "\t") || startsWith(trimws(line), "-")) {
      next
    }
    if (!grepl(":", line, fixed = TRUE)) {
      next
    }
    parts <- strsplit(line, ":", fixed = TRUE)[[1]]
    key <- trimws(parts[[1]])
    val <- trimws(paste(parts[-1], collapse = ":"))
    val <- gsub('^["\']|["\']$', "", val)
    if (val %in% c("", "null", "NULL", "~")) {
      out[[key]] <- NULL
    } else if (val %in% c("true", "TRUE")) {
      out[[key]] <- TRUE
    } else if (val %in% c("false", "FALSE")) {
      out[[key]] <- FALSE
    } else if (grepl("^-?[0-9]+$", val)) {
      out[[key]] <- as.integer(val)
    } else if (grepl("^-?[0-9]+\\.[0-9]+$", val)) {
      out[[key]] <- as.numeric(val)
    } else {
      out[[key]] <- val
    }
  }
  out
}
