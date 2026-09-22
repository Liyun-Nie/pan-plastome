# Threshold-free continuous-rank hotspot concordance helpers.
# Pure statistics: no cpopvar task_params, no M03/M04 controllers,
# and no legacy sharing-test driver.
#
# Pre-registered roles:
#   Q_high  = primary omnibus (distributed shared high-rank structure)
#   C_global = auxiliary landscape concordance
#   M_max    = auxiliary sparse shared-locus signal
#   T_l      = locus continuous recurrence score; K_l < 2 is audit-only

FROZEN_SHA256 <- list(
  CDS = "6b97c74d50514e391c44c4a67a8aa4a843b3696775458f598f961e84135a6a65",
  IGS = "dd11bcfa56fa34792f36fef79e526f8e5af9de184c028050480fcfe1493254cf"
)

GENUS_PAIRS <- list(
  Glycine = c("Glycine_max", "Glycine_soja"),
  Oryza = c("Oryza_rufipogon", "Oryza_sativa"),
  Solanum = c("Solanum_brevicaule", "Solanum_candolleanum"),
  Gossypium = c("Gossypium_barbadense", "Gossypium_hirsutum")
)

assert_no_cpopvar_attached <- function() {
  if ("package:cpopvar" %in% search()) {
    stop("cpopvar must not be attached for threshold-free analysis")
  }
  loaded <- loadedNamespaces()
  if ("cpopvar" %in% loaded) {
    stop("cpopvar namespace is loaded; unload before threshold-free analysis")
  }
  invisible(TRUE)
}

canonical_igs_id <- function(raw_id) {
  x <- strip_poigs_prefix(raw_id)
  parts <- strsplit(x, "-", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) == 2L) {
      paste(sort(p), collapse = "-")
    } else {
      paste(p, collapse = "-")
    }
  }, character(1))
}

canonical_cds_id <- function(raw_id) {
  as.character(raw_id)
}

standardize_rank_table <- function(df, region_type) {
  region_type <- toupper(as.character(region_type)[[1]])
  if (region_type == "CDS") {
    required <- c("species", "gene", "frequency_per_kb", "variant_count", "region_length")
    missing <- setdiff(required, names(df))
    if (length(missing) > 0L) {
      stop("CDS ranked table missing columns: ", paste(missing, collapse = ", "))
    }
    if (!all(tolower(df$region_type) == "cds")) {
      stop("CDS ranked table is not CDS-only")
    }
    if (!all(tolower(df$var_type) %in% c("snp", "snv"))) {
      stop("CDS ranked table is not SNV/snp-only")
    }
    locus_raw <- as.character(df$gene)
    locus_canonical <- canonical_cds_id(locus_raw)
    locus_display <- locus_canonical
    region_length <- as.numeric(df$region_length)
    length_source <- ifelse(is.finite(region_length) & region_length > 0,
                            "ranked_table", "missing")
  } else if (region_type == "IGS") {
    required <- c("species", "poiGS_ID", "frequency_per_kb", "variant_count")
    missing <- setdiff(required, names(df))
    if (length(missing) > 0L) {
      stop("IGS ranked table missing columns: ", paste(missing, collapse = ", "))
    }
    if (!all(toupper(df$region_type) == "IGS")) {
      stop("IGS ranked table is not IGS-only")
    }
    locus_raw <- as.character(df$poiGS_ID)
    locus_canonical <- canonical_igs_id(locus_raw)
    locus_display <- strip_poigs_prefix(locus_raw)
    freq <- as.numeric(df$frequency_per_kb)
    count <- as.numeric(df$variant_count)
    region_length <- ifelse(is.finite(freq) & freq > 0 & is.finite(count) & count > 0,
                            count * 1000 / freq, NA_real_)
    length_source <- ifelse(is.finite(region_length),
                            "recovered_count_over_frequency",
                            "unknown_zero_or_missing")
  } else {
    stop("region_type must be CDS or IGS")
  }
  freq <- as.numeric(df$frequency_per_kb)
  count <- as.numeric(df$variant_count)
  if (any(!is.finite(freq))) {
    stop(region_type, " frequency_per_kb contains non-finite values; missing must not appear as NA in the long table")
  }
  if (any(!is.finite(count))) {
    stop(region_type, " variant_count contains non-finite values")
  }
  if (any(freq < 0) || any(count < 0)) {
    stop(region_type, " has negative frequency or variant_count")
  }
  zero_freq <- freq == 0
  zero_count <- count == 0
  if (any(zero_freq != zero_count)) {
    stop(region_type, " zero frequency and zero variant_count disagree")
  }
  data.frame(
    species = as.character(df$species),
    locus_raw = locus_raw,
    locus_canonical = locus_canonical,
    locus_display = locus_display,
    region_type = region_type,
    frequency_per_kb = freq,
    variant_count = count,
    region_length = as.numeric(region_length),
    length_source = length_source,
    is_zero = freq == 0,
    analyzable = TRUE,
    stringsAsFactors = FALSE
  )
}

audit_key_conflicts <- function(std_df) {
  keys <- paste(std_df$species, std_df$locus_canonical, sep = "\t")
  dup <- duplicated(keys) | duplicated(keys, fromLast = TRUE)
  if (!any(dup)) {
    return(std_df[0, , drop = FALSE])
  }
  std_df[dup, , drop = FALSE]
}

stop_on_key_conflicts <- function(std_df, region_type) {
  conflicts <- audit_key_conflicts(std_df)
  if (nrow(conflicts) > 0L) {
    stop(region_type, " species x canonical locus is not unique; n_conflict_rows=",
         nrow(conflicts), ". Formal analysis halted.")
  }
  invisible(TRUE)
}

build_species_locus_matrices <- function(std_df) {
  stop_on_key_conflicts(std_df, unique(std_df$region_type)[[1]])
  species <- sort(unique(as.character(std_df$species)))
  loci <- sort(unique(as.character(std_df$locus_canonical)))
  n_s <- length(species)
  n_l <- length(loci)
  empty <- matrix(NA_real_, n_s, n_l, dimnames = list(species, loci))
  freq <- empty
  count <- empty
  length_mat <- empty
  idx <- cbind(match(std_df$species, species), match(std_df$locus_canonical, loci))
  if (anyNA(idx)) {
    stop("Failed to map species or locus onto matrix indices")
  }
  freq[idx] <- std_df$frequency_per_kb
  count[idx] <- std_df$variant_count
  length_mat[idx] <- std_df$region_length
  analyzable <- is.finite(freq)
  list(
    species = species,
    loci = loci,
    frequency = freq,
    variant_count = count,
    region_length = length_mat,
    analyzable = analyzable,
    n_species = n_s,
    n_loci = n_l,
    n_analyzable = sum(analyzable),
    n_zero = sum(analyzable & freq == 0),
    n_positive = sum(analyzable & freq > 0),
    n_missing = sum(!analyzable)
  )
}

midrank_normal_scores <- function(x) {
  ok <- is.finite(x)
  n <- as.integer(sum(ok))
  out_r <- rep(NA_real_, length(x))
  out_u <- out_r
  out_z <- out_r
  if (n < 1L) {
    return(list(rank = out_r, u = out_u, z = out_z, n = n))
  }
  r <- rank(x[ok], ties.method = "average")
  u <- (r - 0.5) / n
  z <- stats::qnorm(u)
  if (any(!is.finite(z))) {
    stop("Non-finite normal scores; midrank percentiles must stay inside (0,1)")
  }
  out_r[ok] <- r
  out_u[ok] <- u
  out_z[ok] <- z
  list(rank = out_r, u = out_u, z = out_z, n = n)
}

score_matrices_by_species <- function(value_mat) {
  n_s <- nrow(value_mat)
  n_l <- ncol(value_mat)
  r_mat <- matrix(NA_real_, n_s, n_l, dimnames = dimnames(value_mat))
  u_mat <- r_mat
  z_mat <- r_mat
  n_analyzable <- integer(n_s)
  names(n_analyzable) <- rownames(value_mat)
  for (i in seq_len(n_s)) {
    sc <- midrank_normal_scores(value_mat[i, ])
    r_mat[i, ] <- sc$rank
    u_mat[i, ] <- sc$u
    z_mat[i, ] <- sc$z
    n_analyzable[i] <- sc$n
  }
  list(rank = r_mat, u = u_mat, z = z_mat, n_analyzable = n_analyzable)
}

locus_T <- function(score_mat, k_min = 2L) {
  K <- as.integer(colSums(is.finite(score_mat)))
  names(K) <- colnames(score_mat)
  sum_score <- colSums(score_mat, na.rm = TRUE)
  T <- sum_score / sqrt(K)
  T[!is.finite(K) | K < as.integer(k_min)] <- NA_real_
  list(T = T, K = K, sum_score = sum_score)
}

Q_high <- function(T) {
  T <- as.numeric(T)
  sum(pmax(T, 0)^2, na.rm = TRUE)
}

M_max <- function(T) {
  T <- as.numeric(T)
  ok <- is.finite(T)
  if (!any(ok)) {
    return(NA_real_)
  }
  max(T[ok])
}

q_contribution_table <- function(locus_tbl, z_mat, alpha = 0.05) {
  required <- c(
    "region_type", "locus_canonical", "locus_display", "K_l",
    "included_in_cross_species", "T_l", "mean_percentile_u",
    "p_marginal", "p_bh", "p_maxT"
  )
  missing <- setdiff(required, names(locus_tbl))
  if (length(missing) > 0L) {
    stop("locus_tbl missing columns: ", paste(missing, collapse = ", "))
  }
  if (is.null(colnames(z_mat)) || is.null(rownames(z_mat))) {
    stop("z_mat must have species row names and locus column names")
  }
  keep <- locus_tbl$included_in_cross_species & is.finite(locus_tbl$T_l)
  out <- locus_tbl[keep, required, drop = FALSE]
  idx <- match(out$locus_canonical, colnames(z_mat))
  if (anyNA(idx)) {
    stop("Failed to match locus_tbl loci to z_mat")
  }
  z_sub <- z_mat[, idx, drop = FALSE]
  positive_species <- vapply(seq_len(ncol(z_sub)), function(j) {
    z <- z_sub[, j]
    paste(rownames(z_sub)[is.finite(z) & z > 0], collapse = ",")
  }, character(1))
  n_positive <- colSums(is.finite(z_sub) & z_sub > 0)
  n_zero <- colSums(is.finite(z_sub) & z_sub == 0)
  q_component <- pmax(out$T_l, 0)^2
  q_total <- sum(q_component)
  q_fraction <- if (q_total > 0) q_component / q_total else rep(NA_real_, nrow(out))
  ord <- order(-q_component, -out$T_l, out$locus_canonical)
  out <- out[ord, , drop = FALSE]
  q_component <- q_component[ord]
  q_fraction <- q_fraction[ord]
  n_positive <- as.integer(n_positive[ord])
  n_zero <- as.integer(n_zero[ord])
  positive_species <- positive_species[ord]
  contribution_rank <- rep(NA_integer_, nrow(out))
  positive_component <- q_component > 0
  contribution_rank[positive_component] <- rank(
    -q_component[positive_component],
    ties.method = "min"
  )
  max_t <- max(out$T_l)
  data.frame(
    region_type = out$region_type,
    Q_contribution_rank = contribution_rank,
    locus_canonical = out$locus_canonical,
    locus_display = out$locus_display,
    K_l = as.integer(out$K_l),
    n_species_above_median = n_positive,
    n_species_at_median = n_zero,
    species_above_median = positive_species,
    T_l = out$T_l,
    mean_percentile_u = out$mean_percentile_u,
    Q_component = q_component,
    Q_total = rep(q_total, nrow(out)),
    Q_contribution_fraction = q_fraction,
    Q_contribution_percent = 100 * q_fraction,
    Q_cumulative_fraction = cumsum(q_fraction),
    Q_cumulative_percent = 100 * cumsum(q_fraction),
    is_M_max_locus = abs(out$T_l - max_t) < sqrt(.Machine$double.eps),
    p_marginal = out$p_marginal,
    p_bh = out$p_bh,
    p_maxT = out$p_maxT,
    maxT_significant = is.finite(out$p_maxT) & out$p_maxT < alpha,
    stringsAsFactors = FALSE
  )
}

pairwise_rank_concordance <- function(u_mat, min_overlap = 3L) {
  species <- rownames(u_mat)
  n_s <- length(species)
  rows <- vector("list", n_s * (n_s - 1L) / 2L)
  k <- 1L
  csum <- 0
  wsum <- 0
  for (i in seq_len(n_s - 1L)) {
    for (j in (i + 1L):n_s) {
      ok <- is.finite(u_mat[i, ]) & is.finite(u_mat[j, ])
      n_ij <- as.integer(sum(ok))
      r <- NA_real_
      if (n_ij >= as.integer(min_overlap)) {
        r <- suppressWarnings(stats::cor(u_mat[i, ok], u_mat[j, ok], method = "pearson"))
      }
      rows[[k]] <- data.frame(
        species_a = species[i],
        species_b = species[j],
        n_overlap = n_ij,
        r_midrank_pearson = r,
        stringsAsFactors = FALSE
      )
      if (is.finite(r)) {
        csum <- csum + r * n_ij
        wsum <- wsum + n_ij
      }
      k <- k + 1L
    }
  }
  pairs <- do.call(rbind, rows)
  list(
    C = if (wsum > 0) csum / wsum else NA_real_,
    pairs = pairs,
    weight_sum = wsum,
    min_overlap = as.integer(min_overlap)
  )
}

C_global <- function(u_mat, min_overlap = 3L) {
  pairwise_rank_concordance(u_mat, min_overlap = min_overlap)$C
}

empirical_p_plus_one <- function(n_ge, n_perm) {
  n_ge <- as.numeric(n_ge)
  B <- as.numeric(n_perm)
  (1 + n_ge) / (B + 1)
}

monte_carlo_se <- function(p, n_perm) {
  p <- as.numeric(p)
  B <- as.numeric(n_perm)
  sqrt(p * (1 - p) / (B + 1))
}

row_pearson <- function(A, B) {
  if (ncol(A) < 3L) {
    return(rep(NA_real_, nrow(A)))
  }
  Am <- A - rowMeans(A)
  Bm <- B - rowMeans(B)
  num <- rowSums(Am * Bm)
  den <- sqrt(rowSums(Am * Am) * rowSums(Bm * Bm))
  ifelse(den > 0, num / den, NA_real_)
}

draw_perm_rows <- function(B, n_loci, col_idx, values) {
  n_s <- length(values)
  mat <- matrix(NA_real_, nrow = B, ncol = n_loci)
  if (n_s <= 0L) {
    return(mat)
  }
  if (n_s == 1L) {
    mat[, col_idx] <- values
    return(mat)
  }
  for (b in seq_len(B)) {
    mat[b, col_idx] <- values[sample.int(n_s)]
  }
  mat
}

precompute_species_slots <- function(score_mat) {
  n_s <- nrow(score_mat)
  out <- vector("list", n_s)
  for (i in seq_len(n_s)) {
    ok <- which(is.finite(score_mat[i, ]))
    out[[i]] <- list(
      col_idx = ok,
      values = score_mat[i, ok]
    )
  }
  out
}

precompute_overlaps <- function(u_mat, min_overlap = 3L) {
  n_s <- nrow(u_mat)
  pairs <- list()
  k <- 1L
  for (i in seq_len(n_s - 1L)) {
    for (j in (i + 1L):n_s) {
      ok <- which(is.finite(u_mat[i, ]) & is.finite(u_mat[j, ]))
      if (length(ok) >= as.integer(min_overlap)) {
        pairs[[k]] <- list(i = i, j = j, ok = ok, w = length(ok))
        k <- k + 1L
      }
    }
  }
  pairs
}

permute_omnibus <- function(z_mat,
                            u_mat,
                            n_perm,
                            seed,
                            chunk_size = 250L,
                            k_min = 2L,
                            min_overlap = 3L,
                            compute_C = TRUE,
                            verbose = FALSE,
                            label = "") {
  n_perm <- as.integer(n_perm)
  chunk_size <- as.integer(chunk_size)
  k_min <- as.integer(k_min)
  if (n_perm < 1L) {
    stop("n_perm must be positive")
  }
  n_l <- ncol(z_mat)
  obs_T <- locus_T(z_mat, k_min = k_min)
  T_obs <- obs_T$T
  K <- obs_T$K
  Q_obs <- Q_high(T_obs)
  M_obs <- M_max(T_obs)
  C_obs <- if (isTRUE(compute_C)) C_global(u_mat, min_overlap = min_overlap) else NA_real_
  sqrt_K <- sqrt(K)
  sqrt_K[K < k_min] <- NA_real_
  finite_T <- is.finite(T_obs)
  slots_z <- precompute_species_slots(z_mat)
  slots_u <- precompute_species_slots(u_mat)
  overlaps <- if (isTRUE(compute_C)) {
    precompute_overlaps(u_mat, min_overlap = min_overlap)
  } else {
    list()
  }
  pair_keys <- vapply(overlaps, function(pr) {
    paste(rownames(u_mat)[pr$i], rownames(u_mat)[pr$j], sep = "\t")
  }, character(1))
  pair_obs <- vapply(overlaps, function(pr) {
    suppressWarnings(stats::cor(
      u_mat[pr$i, pr$ok],
      u_mat[pr$j, pr$ok],
      method = "pearson"
    ))
  }, numeric(1))
  names(pair_obs) <- pair_keys
  n_ge_Q <- 0
  n_ge_C <- 0
  n_ge_M <- 0
  n_ge_pair <- integer(length(overlaps))
  names(n_ge_pair) <- pair_keys
  n_ge_T <- rep(0, n_l)
  n_ge_maxT <- rep(0, n_l)
  names(n_ge_T) <- colnames(z_mat)
  names(n_ge_maxT) <- colnames(z_mat)
  Q_null <- rep(NA_real_, n_perm)
  C_null <- rep(NA_real_, n_perm)
  M_null <- rep(NA_real_, n_perm)
  set.seed(as.integer(seed))
  pos <- 1L
  while (pos <= n_perm) {
    B <- min(chunk_size, n_perm - pos + 1L)
    Zsum <- matrix(0, nrow = B, ncol = n_l)
    Ulist <- vector("list", length(slots_u))
    for (s in seq_along(slots_z)) {
      col_idx <- slots_z[[s]]$col_idx
      z_vals <- slots_z[[s]]$values
      u_vals <- slots_u[[s]]$values
      if (length(col_idx) != length(u_vals)) {
        stop("z and u analyzable slots disagree for species index ", s)
      }
      Zs <- matrix(NA_real_, nrow = B, ncol = n_l)
      Us <- matrix(NA_real_, nrow = B, ncol = n_l)
      n_ok <- length(col_idx)
      if (n_ok == 1L) {
        Zs[, col_idx] <- z_vals
        if (isTRUE(compute_C)) {
          Us[, col_idx] <- u_vals
        }
      } else if (n_ok > 1L && isTRUE(compute_C)) {
        for (b in seq_len(B)) {
          ord <- sample.int(n_ok)
          Zs[b, col_idx] <- z_vals[ord]
          Us[b, col_idx] <- u_vals[ord]
        }
      } else if (n_ok > 1L) {
        for (b in seq_len(B)) {
          Zs[b, col_idx] <- z_vals[sample.int(n_ok)]
        }
      }
      if (isTRUE(compute_C)) {
        Ulist[[s]] <- Us
      }
      Zs[!is.finite(Zs)] <- 0
      Zsum <- Zsum + Zs
    }
    T_chunk <- sweep(Zsum, 2, sqrt_K, "/")
    T_chunk[, !finite_T] <- NA_real_
    Q_chunk <- rowSums(pmax(T_chunk, 0)^2, na.rm = TRUE)
    M_chunk <- apply(T_chunk, 1, function(v) {
      ok <- is.finite(v)
      if (!any(ok)) NA_real_ else max(v[ok])
    })
    if (isTRUE(compute_C)) {
      c_acc <- rep(0, B)
      w_acc <- rep(0, B)
      for (pair_idx in seq_along(overlaps)) {
        pr <- overlaps[[pair_idx]]
        r <- row_pearson(Ulist[[pr$i]][, pr$ok, drop = FALSE],
                         Ulist[[pr$j]][, pr$ok, drop = FALSE])
        fin <- is.finite(r)
        c_acc[fin] <- c_acc[fin] + r[fin] * pr$w
        w_acc[fin] <- w_acc[fin] + pr$w
        if (is.finite(pair_obs[pair_idx])) {
          n_ge_pair[pair_idx] <- n_ge_pair[pair_idx] +
            sum(fin & r >= pair_obs[pair_idx])
        }
      }
      C_chunk <- ifelse(w_acc > 0, c_acc / w_acc, NA_real_)
    } else {
      C_chunk <- rep(NA_real_, B)
    }
    idx <- pos:(pos + B - 1L)
    Q_null[idx] <- Q_chunk
    C_null[idx] <- C_chunk
    M_null[idx] <- M_chunk
    n_ge_Q <- n_ge_Q + sum(Q_chunk >= Q_obs, na.rm = TRUE)
    if (isTRUE(compute_C)) {
      n_ge_C <- n_ge_C + sum(C_chunk >= C_obs, na.rm = TRUE)
    }
    n_ge_M <- n_ge_M + sum(M_chunk >= M_obs, na.rm = TRUE)
    ge_T <- sweep(T_chunk, 2, T_obs, ">=")
    ge_T[, !finite_T] <- FALSE
    n_ge_T <- n_ge_T + colSums(ge_T, na.rm = TRUE)
    ge_max <- matrix(M_chunk, nrow = B, ncol = n_l) >= matrix(T_obs, nrow = B, ncol = n_l, byrow = TRUE)
    ge_max[, !finite_T] <- FALSE
    n_ge_maxT <- n_ge_maxT + colSums(ge_max, na.rm = TRUE)
    if (isTRUE(verbose)) {
      cat(sprintf("  %s permutation %d / %d\n", label, pos + B - 1L, n_perm))
    }
    pos <- pos + B
  }
  p_Q <- empirical_p_plus_one(n_ge_Q, n_perm)
  p_C <- if (isTRUE(compute_C)) empirical_p_plus_one(n_ge_C, n_perm) else NA_real_
  p_M <- empirical_p_plus_one(n_ge_M, n_perm)
  p_pair <- if (isTRUE(compute_C)) {
    empirical_p_plus_one(n_ge_pair, n_perm)
  } else {
    numeric()
  }
  p_pair[!is.finite(pair_obs)] <- NA_real_
  names(p_pair) <- pair_keys
  p_T <- empirical_p_plus_one(n_ge_T, n_perm)
  p_maxT <- empirical_p_plus_one(n_ge_maxT, n_perm)
  p_T[!finite_T] <- NA_real_
  p_maxT[!finite_T] <- NA_real_
  names(p_T) <- colnames(z_mat)
  names(p_maxT) <- colnames(z_mat)
  names(n_ge_T) <- colnames(z_mat)
  names(n_ge_maxT) <- colnames(z_mat)
  p_T_bh <- rep(NA_real_, n_l)
  names(p_T_bh) <- colnames(z_mat)
  if (any(finite_T)) {
    p_T_bh[finite_T] <- stats::p.adjust(p_T[finite_T], method = "BH")
  }
  list(
    n_perm = n_perm,
    seed = as.integer(seed),
    k_min = k_min,
    T_obs = T_obs,
    K = K,
    Q_obs = Q_obs,
    C_obs = C_obs,
    M_obs = M_obs,
    n_ge_Q = n_ge_Q,
    n_ge_C = n_ge_C,
    n_ge_M = n_ge_M,
    pair_keys = pair_keys,
    pair_obs = pair_obs,
    n_ge_pair = n_ge_pair,
    p_pair = p_pair,
    se_pair = monte_carlo_se(p_pair, n_perm),
    n_ge_T = n_ge_T,
    n_ge_maxT = n_ge_maxT,
    p_Q = p_Q,
    p_C = p_C,
    p_M = p_M,
    p_T = p_T,
    p_T_bh = p_T_bh,
    p_maxT = p_maxT,
    se_Q = monte_carlo_se(p_Q, n_perm),
    se_C = if (isTRUE(compute_C)) monte_carlo_se(p_C, n_perm) else NA_real_,
    se_M = monte_carlo_se(p_M, n_perm),
    se_T = monte_carlo_se(p_T, n_perm),
    se_maxT = monte_carlo_se(p_maxT, n_perm),
    Q_null = Q_null,
    C_null = C_null,
    M_null = M_null
  )
}

centered_u_matrix <- function(u_mat) {
  u_mat - 0.5
}

centered_rank_matrix <- function(rank_mat, n_analyzable) {
  out <- rank_mat
  for (i in seq_len(nrow(rank_mat))) {
    n <- n_analyzable[i]
    if (n > 0L) {
      ok <- is.finite(rank_mat[i, ])
      out[i, ok] <- rank_mat[i, ok] - (n + 1) / 2
    }
  }
  out
}

prevalence_score_matrix <- function(freq_mat) {
  out <- matrix(NA_real_, nrow = nrow(freq_mat), ncol = ncol(freq_mat),
                dimnames = dimnames(freq_mat))
  ok <- is.finite(freq_mat)
  out[ok] <- as.numeric(freq_mat[ok] > 0)
  out
}

positive_intensity_matrix <- function(freq_mat) {
  out <- freq_mat
  out[is.finite(freq_mat) & freq_mat == 0] <- NA_real_
  out
}

length_residual_matrix <- function(count_mat, length_mat) {
  out <- matrix(NA_real_, nrow = nrow(count_mat), ncol = ncol(count_mat),
                dimnames = dimnames(count_mat))
  for (i in seq_len(nrow(count_mat))) {
    ok <- is.finite(count_mat[i, ]) & is.finite(length_mat[i, ]) & length_mat[i, ] > 0
    if (!any(ok)) {
      next
    }
    count <- count_mat[i, ok]
    len <- length_mat[i, ok]
    rate <- sum(count) / sum(len)
    mu <- rate * len
    denom <- sqrt(pmax(mu, 0))
    resid <- ifelse(denom > 0, (count - mu) / denom, 0)
    out[i, ok] <- resid
  }
  out
}

genus_downsample_combos <- function(genus_pairs = GENUS_PAIRS) {
  genera <- names(genus_pairs)
  grid <- expand.grid(lapply(genus_pairs, function(x) x), stringsAsFactors = FALSE)
  names(grid) <- genera
  combos <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    kept_pair <- as.character(unlist(grid[i, , drop = TRUE]))
    dropped <- character()
    for (g in genera) {
      dropped <- c(dropped, setdiff(genus_pairs[[g]], kept_pair))
    }
    combos[[i]] <- list(
      combo_id = i,
      kept_representatives = kept_pair,
      dropped_species = dropped,
      label = paste(sprintf("%s=%s", genera, kept_pair), collapse = ";")
    )
  }
  combos
}

subset_species_matrix <- function(mat, drop_species) {
  keep <- setdiff(rownames(mat), drop_species)
  mat[keep, , drop = FALSE]
}

locus_mean_percentile <- function(u_mat, k_min = 2L) {
  K <- as.integer(colSums(is.finite(u_mat)))
  mu <- colMeans(u_mat, na.rm = TRUE)
  mu[K < as.integer(k_min)] <- NA_real_
  list(mean_u = mu, K = K)
}

covered_species_string <- function(score_mat, locus) {
  z <- score_mat[, locus]
  paste(rownames(score_mat)[is.finite(z)], collapse = ",")
}

global_result_row <- function(region_type, analysis, statistic, observed, n_ge, p, se,
                              n_perm, seed, role, tail = "upper") {
  data.frame(
    region_type = region_type,
    analysis = analysis,
    statistic = statistic,
    role = role,
    tail = tail,
    observed = observed,
    n_ge = as.integer(n_ge),
    n_perm = as.integer(n_perm),
    seed = as.integer(seed),
    p_empirical = p,
    mc_se = se,
    stringsAsFactors = FALSE
  )
}

write_simple_yaml <- function(x, path, indent = 0L) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lines <- yaml_lines(x, indent)
  writeLines(lines, path, useBytes = FALSE)
}

yaml_scalar <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return("~")
  }
  if (is.logical(x)) {
    return(if (isTRUE(x)) "true" else "false")
  }
  if (is.numeric(x)) {
    if (!is.finite(x)) return(as.character(x))
    return(format(x, digits = 15, scientific = FALSE, trim = TRUE))
  }
  val <- as.character(x)
  val <- gsub("\\", "\\\\", val, fixed = TRUE)
  val <- gsub("\"", "\\\"", val, fixed = TRUE)
  if (grepl("[:#\\n]", val) || grepl(" ", val) || val == "") {
    return(paste0('"', val, '"'))
  }
  val
}

yaml_lines <- function(x, indent = 0L) {
  pad <- paste(rep("  ", indent), collapse = "")
  if (!is.list(x)) {
    if (length(x) > 1L) {
      return(paste0(pad, "- ", vapply(x, yaml_scalar, character(1))))
    }
    return(paste0(pad, yaml_scalar(x)))
  }
  nms <- names(x)
  if (is.null(nms)) {
    out <- character()
    for (item in x) {
      if (is.list(item)) {
        out <- c(out, paste0(pad, "-"), yaml_lines(item, indent + 1L))
      } else {
        out <- c(out, paste0(pad, "- ", yaml_scalar(item)))
      }
    }
    return(out)
  }
  out <- character()
  for (nm in nms) {
    val <- x[[nm]]
    if (is.list(val)) {
      out <- c(out, paste0(pad, nm, ":"))
      out <- c(out, yaml_lines(val, indent + 1L))
    } else if (length(val) > 1L) {
      out <- c(out, paste0(pad, nm, ":"))
      out <- c(out, paste0(pad, "  - ", vapply(val, yaml_scalar, character(1))))
    } else {
      out <- c(out, paste0(pad, nm, ": ", yaml_scalar(val)))
    }
  }
  out
}

rank_matrix_long <- function(std_df, mats, scores, region_type) {
  species <- mats$species
  loci <- mats$loci
  grid <- expand.grid(species = species, locus_canonical = loci, stringsAsFactors = FALSE)
  sp <- as.character(grid$species)
  loc <- as.character(grid$locus_canonical)
  idx <- cbind(sp, loc)
  std_key <- paste(std_df$species, std_df$locus_canonical, sep = "\t")
  j <- match(paste(sp, loc, sep = "\t"), std_key)
  analyzable <- as.logical(mats$analyzable[idx])
  freq <- as.numeric(mats$frequency[idx])
  is_zero <- freq == 0
  is_zero[!analyzable] <- NA
  data.frame(
    region_type = region_type,
    species = sp,
    locus_canonical = loc,
    locus_raw = std_df$locus_raw[j],
    locus_display = std_df$locus_display[j],
    analyzable = analyzable,
    frequency_per_kb = freq,
    variant_count = as.numeric(mats$variant_count[idx]),
    region_length = as.numeric(mats$region_length[idx]),
    is_zero = is_zero,
    is_missing = !analyzable,
    midrank = scores$rank[idx],
    u_percentile = scores$u[idx],
    z_normal = scores$z[idx],
    length_source = std_df$length_source[j],
    stringsAsFactors = FALSE
  )
}

audit_rows <- function(region_type, path, sha, std_df, mats, scores, k_min = 2L) {
  K <- as.integer(colSums(mats$analyzable))
  names(K) <- mats$loci
  k1 <- names(K)[K < as.integer(k_min)]
  n_prefix <- if (region_type == "IGS") sum(grepl("^poiGS_", std_df$locus_raw)) else NA_integer_
  n_alias <- NA_integer_
  if (region_type == "IGS") {
    n_alias <- nrow(unique(std_df[, c("locus_raw", "locus_canonical")])) - length(unique(std_df$locus_canonical))
  }
  length_ok <- is.finite(std_df$region_length) & std_df$region_length > 0
  info <- file.info(path)
  data.frame(
    role = paste0(tolower(region_type), "_full_ranked"),
    region_type = region_type,
    input_path = path,
    exists = file.exists(path),
    bytes = if (is.na(info$size)) NA_real_ else as.numeric(info$size),
    mtime_unix = if (is.na(info$mtime)) NA_real_ else as.numeric(info$mtime),
    sha256 = sha,
    n_long_rows = nrow(std_df),
    n_species = mats$n_species,
    n_loci = mats$n_loci,
    n_analyzable = mats$n_analyzable,
    n_zero = mats$n_zero,
    n_positive = mats$n_positive,
    n_missing_cells = mats$n_missing,
    n_duplicate_keys = nrow(audit_key_conflicts(std_df)),
    n_loci_k_lt_min = length(k1),
    k_min_cross_species = as.integer(k_min),
    k1_loci = paste(k1, collapse = ","),
    n_poigs_prefix = n_prefix,
    n_direction_alias_collapses = n_alias,
    n_rows_with_length = sum(length_ok),
    n_rows_length_unknown = sum(!length_ok),
    stringsAsFactors = FALSE
  )
}
