library(tidyverse)
library(MMUPHin)
library(metafor)

code_dir <- file.path(
  base_symbol,
  "downstream_analysis/src/github_code_availability/submission_version_261102/FigS2_sensitivity_analysis"
)
sensitivity_root <- file.path(
  base_symbol,
  "downstream_analysis/output/submission_version_261102/FigS2_sensitivity_analysis"
)
source(file.path(code_dir, "00.sensitivity_functions.R"))

analysis_root <- file.path(sensitivity_root, "downsampling")
dir.create(
  file.path(analysis_root, "results"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(analysis_root, "generated_data"),
  recursive = TRUE,
  showWarnings = FALSE
)
seed_base <- 2026102900L
n_rep <- 1000L
f_levels <- c(1, 0.75, 0.50, 0.38, 0.25)

# Fixed candidate universe and historical masking; only LO is subsampled.
regression_geometry <- function(x, group_column = 2L) {
  stopifnot(all(is.finite(x)), nrow(x) > ncol(x))
  scale <- sqrt(colSums(x^2))
  scale[scale == 0] <- 1
  xx <- sweep(x, 2, scale, "/")
  q <- qr(xx, tol = 1e-09)
  qn <- qr(xx[, -group_column, drop = FALSE], tol = 1e-09)
  stopifnot(q$rank == qn$rank + 1L)
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  gr <- as.numeric(qr.resid(qn, x[, group_column]))
  info <- sum(gr^2)
  stopifnot(info > 1e-12)
  a <- gr / info
  h <- rowSums(Q^2)

  list(
    x = x,
    Q = Q,
    a = a,
    aQ = as.numeric(crossprod(Q, a)),
    information = info,
    h = h,
    rank = q$rank,
    df = nrow(x) -
      q$rank
  )
}

# Pool candidate effects with HS, then apply the original QC and BH rules.
hs_qc <- function(beta, variance) {
  good <- is.finite(beta) & is.finite(variance) & variance > 1e-16
  b <- beta
  b[!good] <- 0
  w <- 1 / variance
  w[!good] <- 0
  k <- colSums(good)
  sw <- colSums(w)
  fe <- colSums(w * b) / sw
  Q <- colSums(w * sweep(b, 2, fe, "-")^2)
  tau2 <- pmax(0, (Q - k) / sw)
  wr <- 1 / sweep(variance, 2, tau2, "+")
  wr[!good] <- 0
  swr <- colSums(wr)
  pooled <- colSums(wr * b) / swr
  se <- sqrt(1 / swr)
  weights <- sweep(wr, 2, swr, "/")
  max_weight <- weights[1, ]
  if (nrow(weights) > 1L)

    for (i in 2:nrow(weights)) max_weight <- pmax(max_weight, weights[i, ], na.rm = TRUE)
  vt <- (k - 1) * sw / (sw^2 - colSums(w^2))
  I2 <- 100 * tau2 / (tau2 + vt)
  I2[k < 2L | !is.finite(I2)] <- NA_real_
  pooled[k < 2L] <- se[k < 2L] <- NA_real_
  nonzero <- good & b != 0
  nz <- colSums(nonzero)
  same <- sweep(sign(b), 2, sign(pooled), "==") & nonzero
  consistency <- colSums(same, na.rm = TRUE) / nz
  pass <- nz >= 3L &
    max_weight < 0.4 &
    is.finite(pooled) &
    ((I2 < 50 &
      consistency >= 0.6) |
      (I2 >= 50 &
        I2 <
          80 & consistency >= 0.75))
  pass[is.na(pass)] <- FALSE
  p <- 2 * pnorm(-abs(pooled / se))
  q <- rep(NA_real_, length(p))
  q[pass] <- p.adjust(p[pass], "BH")
  selected <- pass & is.finite(q) & q < 0.1

  list(
    beta = pooled,
    se = se,
    k = k,
    tau2 = tau2,
    I2 = I2,
    weights = weights,
    max_weight = max_weight,
    n_nonzero = nz,
    consistency = consistency,
    pass = pass,
    p = p,
    q = q,
    selected = selected
  )
}

# Refit one sampled block with the original masking and zero replacement.
fit_block_historical <- function(b, ix) {
  p <- length(b$features)
  beta <- rep(NA_real_, p)
  variance <- rep(NA_real_, p)
  raw <- b$raw[ix, , drop = FALSE]
  count <- colSums(raw > 0)
  masked <- raw
  masked[, count > 0 & count < 5] <- 0
  floor <- apply(
    masked,
    2L,
    function(z) if (any(z > 0))
      min(z[z > 0]) / 2
    else NA_real_
  )
  y <- log2(pmax(masked, matrix(floor, nrow(masked), p, byrow = TRUE)))
  y[, !is.finite(floor)] <- 0
  x <- b$x[ix, , drop = FALSE]
  group <- b$metadata$Group[ix]
  cell <- table(factor(group, levels = c("Control", "CRC")))
  active <- count >= 5L & is.finite(floor)
  reason <- rep("inactive", p)
  if (any(cell < 3L) || nrow(x) <= ncol(x)) {
    reason[active] <- "design_or_cell_failure"
    return(list(
      beta = beta,
      variance = variance,
      active = active,
      reason = reason,
      n = nrow(x),
      df = NA_integer_
    ))
  }

  g <- tryCatch(regression_geometry(x), error = function(e) NULL)
  if (is.null(g)) {
    reason[active] <- "design_rank_failure"
    return(list(
      beta = beta,
      variance = variance,
      active = active,
      reason = reason,
      n = nrow(x),
      df = NA_integer_
    ))
  }

  if (g$df < 1L) {
    reason[active] <- "no_residual_df"
    return(list(
      beta = beta,
      variance = variance,
      active = active,
      reason = reason,
      n = nrow(x),
      df = g$df
    ))
  }

  ids <- which(active)
  if (!length(ids)) {
    return(list(
      beta = beta,
      variance = variance,
      active = active,
      reason = reason,
      n = nrow(x),
      df = g$df
    ))
  }

  yy <- y[, ids, drop = FALSE]
  projection <- crossprod(g$Q, yy)
  coef <- as.numeric(crossprod(g$aQ, projection))
  residual <- yy - g$Q %*% projection
  v <- colSums(residual^2) / g$df / g$information
  good <- is.finite(coef) & is.finite(v) & v > 1e-16
  beta[ids[good]] <- coef[good]
  variance[ids[good]] <- v[good]
  reason[ids[good]] <- "success"
  reason[ids[!good]] <- "invalid_variance"

  list(
    beta = beta,
    variance = variance,
    active = active,
    reason = reason,
    n = nrow(x),
    df = g$df
  )
}

# Reselect concordant shared taxa and summarise their effect differences.
algorithm_dir <- function(beta, variance, ages) {
  eo <- hs_qc(beta[ages == "EO", , drop = FALSE], variance[ages == "EO", , drop = FALSE])
  lo <- hs_qc(beta[ages == "LO", , drop = FALSE], variance[ages == "LO", , drop = FALSE])
  raw_shared <- which(eo$selected & lo$selected)
  opposite <- raw_shared[sign(eo$beta[raw_shared]) != sign(lo$beta[raw_shared])]
  shared <- setdiff(raw_shared, opposite)
  d <- abs(eo$beta[shared]) - abs(lo$beta[shared])
  T <- if (length(d))
    median(d)
  else NA_real_
  out <- data.frame(
    T = T,
    n_shared = length(shared),
    n_shared_raw = length(raw_shared),
    n_opposite = length(opposite),
    n_EO = sum(eo$selected),
    n_LO = sum(lo$selected),
    n_QC_EO = sum(eo$pass),
    n_QC_LO = sum(lo$pass),
    n_estimable_EO = sum(is.finite(eo$beta)),
    n_estimable_LO = sum(is.finite(lo$beta)),
    median_abs_EO_shared = if (length(shared))
      median(abs(eo$beta[shared]))
    else NA_real_,
    median_abs_LO_shared = if (length(shared))
      median(abs(lo$beta[shared]))
    else NA_real_,
    positive_n = if (length(d))
      sum(d > 0)
    else NA_integer_,
    stringsAsFactors = FALSE
  )
  attr(out, "shared") <- shared
  attr(out, "selected_EO") <- which(eo$selected)
  attr(out, "selected_LO") <- which(lo$selected)
  out
}

# Permute controls and CRC cases separately within each LO cohort.
make_lo_perms <- function(lo_blocks) {
  lapply(
    lo_blocks,
    function(b) {
      g <- b$metadata$Group
      control <- which(g == "Control")
      crc <- which(g == "CRC")

      list(
        control = sample(control),
        crc = sample(crc),
        n_control = length(control),
        n_crc = length(crc)
      )
    }
  )
}

# Measure baseline overlap while retaining failed-draw states.
shared_overlap <- function(n_shared, error, members, baseline) {
  if (isTRUE(error) || is.na(n_shared))
    return(c(n_overlap = NA_real_, jaccard = NA_real_))
  if (n_shared == 0L)
    return(c(n_overlap = 0, jaccard = 0))
  stopifnot(length(unique(members)) == n_shared)
  c(
    n_overlap = length(intersect(members, baseline)),
    jaccard = length(intersect(members, baseline)) / length(union(members, baseline))
  )
}

# Take the requested fraction from each permuted disease group.
lo_indices_for_level <- function(perms, f) {
  lapply(
    perms,
    function(z) {
      nc <- max(1L, round(f * z$n_control))
      nr <- max(1L, round(f * z$n_crc))
      nc <- min(nc, z$n_control)
      nr <- min(nr, z$n_crc)
      c(z$control[seq_len(nc)], z$crc[seq_len(nr)])
    }
  )
}

features <- fit_meta.crc_EO$meta_fits$feature
main_fits <- list(EO = fit_meta.crc_EO, LO = fit_meta.crc_LO)
cohorts <- intersect(names(main_fits$EO$maaslin_fits), names(main_fits$LO$maaslin_fits))
metadata <- metadata_allcohorts %>%
  filter(Group %in% c("Control", "CRC"))
input_ids <- metadata$Sample_ID[metadata$Cohort %in% union(
  names(main_fits$EO$maaslin_fits),
  names(main_fits$LO$maaslin_fits)
)]
raw_matrix <- matrix(0, length(features), length(input_ids), dimnames = list(features, input_ids))

for (tab in meta_list_metadata_intersect_samples) {
  f <- match(features, tab$clade_name)
  ids <- intersect(names(tab), input_ids)
  if (!length(ids)) next
  raw_matrix[!is.na(f), ids] <- as.matrix(as.data.frame(tab)[f[!is.na(f)], ids, drop = FALSE])
}

raw_matrix <- sweep(raw_matrix, 2, colSums(raw_matrix), "/") * 100
blocks <- list()
cohort_sequence <- c(
  "Feng_2015",
  "FUSCC-SHSD",
  "Gupta_2019",
  "Hannigan_2018",
  "Liu_2022",
  "Piccinno_2025_IIGM_CZ",
  "Piccinno_2025_IIGM_IT",
  "Piccinno_2025_IIGM_TU",
  "ThomasAM_2018b",
  "Vogtmann_2016",
  "Wirbel_2019",
  "YachidaS_2019",
  "Yang_2020",
  "Yu_2017",
  "Zeller_2014"
)

for (age in c("EO", "LO")) {
  # Preserve the original cohort order because LO permutations consume RNG sequentially.
  cohort_order <- cohort_sequence[cohort_sequence %in% names(main_fits[[age]]$maaslin_fits)]
  if (age == "EO") cohort_order <- c(
    setdiff(cohort_order, c("Feng_2015", "ThomasAM_2018b")),
    "Feng_2015",
    "ThomasAM_2018b"
  )

  for (cohort in cohort_order[cohort_order %in% cohorts]) {
    d <- metadata %>%
      filter(Cohort == cohort, Age_class == age)
    d$Group <- factor(d$Group, levels = c("Control", "CRC"))
    covariates <- c("Age", "Sex", "BMI")
    covariates <- covariates[vapply(
      d[covariates],
      function(z) length(unique(z[!is.na(z)])) > 1L,
      logical(1)
    )]
    if ("Sex" %in% covariates) d$Sex <- factor(d$Sex)
    formula <- reformulate(c("Group", covariates))
    x <- model.matrix(formula, model.frame(formula, d, na.action = na.pass))
    keep <- complete.cases(x)
    d <- d[keep, ]
    blocks[[paste(age, cohort, sep = "__")]] <- list(
      x = x[keep, , drop = FALSE],
      raw = t(raw_matrix[, d$Sample_ID, drop = FALSE]),
      metadata = d,
      age = age,
      cohort = cohort,
      features = features
    )
  }
}

ages <- vapply(blocks, `[[`, "", "age")
eo_rows <- which(ages == "EO")
lo_rows <- which(ages == "LO")
n_blocks <- length(blocks)
beta <- variance <- matrix(
  NA_real_,
  n_blocks,
  length(features),
  dimnames = list(names(blocks), features)
)

for (i in seq_along(blocks)) {
  fitted <- fit_block_historical(blocks[[i]], seq_len(nrow(blocks[[i]]$x)))
  beta[i, ] <- fitted$beta
  variance[i, ] <- fitted$variance
}

baseline <- algorithm_dir(beta, variance, ages)
baseline_shared <- features[attr(baseline, "shared")]
eo_beta <- beta[eo_rows, , drop = FALSE]
eo_variance <- variance[eo_rows, , drop = FALSE]
lo_blocks <- blocks[lo_rows]
p <- list(
  blocks = blocks,
  features = features,
  cohorts = cohorts,
  ages = ages,
  eo_rows = eo_rows,
  lo_rows = lo_rows,
  eo_beta = eo_beta,
  eo_variance = eo_variance,
  baseline = baseline,
  baseline_shared = baseline_shared
)
counts <- data.frame(
  cohort = vapply(blocks, `[[`, "", "cohort"),
  age = ages,
  n = vapply(blocks, function(b) nrow(b$x), integer(1)),
  Control = vapply(blocks, function(b) sum(b$metadata$Group == "Control"), integer(1)),
  CRC = vapply(blocks, function(b) sum(b$metadata$Group == "CRC"), integer(1))
)
write_csv(counts, file.path(analysis_root, "results/baseline_block_counts.csv"))
stopifnot(length(features) == 1063L, length(cohorts) == 12L, baseline$n_shared == 24L)

base_control <- sum(vapply(lo_blocks, function(b) sum(b$metadata$Group == "Control"), integer(1)))

base_crc <- sum(vapply(lo_blocks, function(b) sum(b$metadata$Group == "CRC"), integer(1)))

rows <- list(data.frame(
  draw = 0L,
  f = 1,
  T = baseline$T,
  n_shared = baseline$n_shared,
  n_shared_raw = baseline$n_shared_raw,
  n_opposite = baseline$n_opposite,
  n_EO = baseline$n_EO,
  n_LO = baseline$n_LO,
  n_QC_EO = baseline$n_QC_EO,
  n_QC_LO = baseline$n_QC_LO,
  n_estimable_EO = baseline$n_estimable_EO,
  n_estimable_LO = baseline$n_estimable_LO,
  median_abs_EO_shared = baseline$median_abs_EO_shared,
  median_abs_LO_shared = baseline$median_abs_LO_shared,
  positive_n = baseline$positive_n,
  n_LO_total = base_control +
    base_crc,
  n_LO_control = base_control,
  n_LO_crc = base_crc,
  error = FALSE,
  error_message = NA_character_,
  stringsAsFactors = FALSE
))

baseline_shared <- p$baseline_shared

write.csv(
  data.frame(feature = baseline_shared),
  file.path(analysis_root, "results/baseline_shared_species.csv"),
  row.names = FALSE
)

membership <- list()
failures <- list()
cohort_counts <- list()
start <- proc.time()[3]

for (draw in seq_len(n_rep)) {
  set.seed(seed_base + draw)
  perms <- make_lo_perms(lo_blocks)

  for (f in f_levels[-1]) {
    res <- tryCatch(
      {
        lo_idx <- lo_indices_for_level(perms, f)
        beta <- matrix(
          NA_real_,
          n_blocks,
          length(features),
          dimnames = list(names(blocks), features)
        )
        variance <- beta
        beta[eo_rows, ] <- eo_beta
        variance[eo_rows, ] <- eo_variance
        failure <- list()

        for (i in seq_along(lo_rows)) {
          fit <- fit_block_historical(lo_blocks[[i]], lo_idx[[i]])
          beta[lo_rows[i], ] <- fit$beta
          variance[lo_rows[i], ] <- fit$variance
          rr <- table(fit$reason)
          failure[[i]] <- data.frame(
            f = f,
            draw = draw,
            cohort = lo_blocks[[i]]$cohort,
            n = fit$n,
            n_active = sum(fit$active),
            design_or_cell_failure = sum(rr["design_or_cell_failure"], na.rm = TRUE),
            design_rank_failure = sum(rr["design_rank_failure"], na.rm = TRUE),
            no_residual_df = sum(rr["no_residual_df"], na.rm = TRUE),
            invalid_variance = sum(rr["invalid_variance"], na.rm = TRUE),
            inactive = sum(rr["inactive"], na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        }

        a <- algorithm_dir(beta, variance, ages)

        list(
          ok = TRUE,
          lo_idx = lo_idx,
          a = a,
          failure = do.call(rbind, failure),
          error = NA_character_
        )
      },
      error = function(e) list(ok = FALSE, error = conditionMessage(e))
    )
    if (!res$ok) {
      rows[[length(rows) + 1L]] <- data.frame(
        draw = draw,
        f = f,
        T = NA_real_,
        n_shared = NA_integer_,
        n_shared_raw = NA_integer_,
        n_opposite = NA_integer_,
        n_EO = NA_integer_,
        n_LO = NA_integer_,
        n_QC_EO = NA_integer_,
        n_QC_LO = NA_integer_,
        n_estimable_EO = NA_integer_,
        n_estimable_LO = NA_integer_,
        median_abs_EO_shared = NA_real_,
        median_abs_LO_shared = NA_real_,
        positive_n = NA_integer_,
        n_LO_total = NA_integer_,
        n_LO_control = NA_integer_,
        n_LO_crc = NA_integer_,
        error = TRUE,
        error_message = res$error,
        stringsAsFactors = FALSE
      )
      next
    }

    a <- res$a
    lo_idx <- res$lo_idx
    n_control <- sum(vapply(
      seq_along(lo_blocks),
      function(i) sum(lo_blocks[[i]]$metadata$Group[lo_idx[[i]]] ==
        "Control"),
      integer(1)
    ))
    n_crc <- sum(vapply(
      seq_along(lo_blocks),
      function(i) sum(lo_blocks[[i]]$metadata$Group[lo_idx[[i]]] ==
        "CRC"),
      integer(1)
    ))
    rows[[length(rows) + 1L]] <- data.frame(
      draw = draw,
      f = f,
      T = a$T,
      n_shared = a$n_shared,
      n_shared_raw = a$n_shared_raw,
      n_opposite = a$n_opposite,
      n_EO = a$n_EO,
      n_LO = a$n_LO,
      n_QC_EO = a$n_QC_EO,
      n_QC_LO = a$n_QC_LO,
      n_estimable_EO = a$n_estimable_EO,
      n_estimable_LO = a$n_estimable_LO,
      median_abs_EO_shared = a$median_abs_EO_shared,
      median_abs_LO_shared = a$median_abs_LO_shared,
      positive_n = a$positive_n,
      n_LO_total = n_control +
        n_crc,
      n_LO_control = n_control,
      n_LO_crc = n_crc,
      error = FALSE,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    failures[[length(failures) + 1L]] <- res$failure
    cohort_counts[[length(cohort_counts) + 1L]] <- data.frame(
      f = f,
      draw = draw,
      cohort = vapply(
        lo_blocks,
        `[[`,
        "",
        "cohort"
      ),
      n_LO = vapply(lo_idx, length, integer(1)),
      n_control = vapply(
        seq_along(lo_blocks),
        function(i) sum(lo_blocks[[i]]$metadata$Group[lo_idx[[i]]] == "Control"),
        integer(1)
      ),
      n_crc = vapply(
        seq_along(lo_blocks),
        function(i) sum(lo_blocks[[i]]$metadata$Group[lo_idx[[i]]] == "CRC"),
        integer(1)
      ),
      stringsAsFactors = FALSE
    )
    sh <- attr(a, "shared")
    if (length(sh)) {
      membership[[length(membership) + 1L]] <- data.frame(
        f = f,
        draw = draw,
        feature = features[sh],
        stringsAsFactors = FALSE
      )
    }
  }

  if (draw %% 100L == 0L)
    cat(
      "DOWNSAMPLE_PROGRESS",
      draw,
      "/",
      n_rep,
      "elapsed_sec",
      round(proc.time()[3] - start, 1),
      "\n"
    )
}

draws <- do.call(rbind, rows)

write.csv(draws, file.path(analysis_root, "results/downsample_draws.csv"), row.names = FALSE)

membership_df <- if (length(membership)) do.call(rbind, membership) else data.frame(
  f = numeric(),
  draw = integer(),
  feature = character()
)

split_key <- paste(draws$f, draws$draw)

overlap <- do.call(
  rbind,
  lapply(
    which(draws$f != 1),
    function(i) {
      z <- draws[i, ]
      a <- unique(membership_df$feature[membership_df$f == z$f & membership_df$draw == z$draw])
      ov <- shared_overlap(z$n_shared, z$error, a, baseline_shared)

      data.frame(
        f = z$f,
        draw = z$draw,
        n_overlap = ov[1],
        jaccard = ov[2],
        row.names = NULL
      )
    }
  )
)

write.csv(
  overlap,
  file.path(analysis_root, "results/downsample_overlap.csv"),
  row.names = FALSE
)

summarise_level <- function(d, f, baseline_T) {
  ok <- is.finite(d$T)
  q <- if (sum(ok) >= 20L)
    quantile(d$T[ok], c(0.025, 0.5, 0.975), names = FALSE)
  else rep(NA_real_, 3)
  med <- if (any(ok))
    median(d$T[ok])
  else NA_real_
  pos <- sum(d$T > 0, na.rm = TRUE)
  neg <- sum(d$T < 0, na.rm = TRUE)
  zero <- sum(d$T == 0, na.rm = TRUE)
  empty <- sum(d$n_shared == 0, na.rm = TRUE)
  ci <- if (nrow(d) > 1L)
    binom.test(pos, nrow(d))$conf.int
  else c(NA_real_, NA_real_)

  data.frame(
    f = f,
    n_rep = nrow(d),
    n_valid = sum(ok),
    n_empty = empty,
    n_error = sum(d$error, na.rm = TRUE),
    T_q025 = q[1],
    T_median = med,
    T_q975 = q[3],
    T_mean = mean(d$T, na.rm = TRUE),
    T_abs_dev_median = median(
      abs(d$T -
        baseline_T),
      na.rm = TRUE
    ),
    positive_n = pos,
    negative_n = neg,
    zero_n = zero,
    positive_fraction_all = pos / nrow(d),
    positive_fraction_nonempty = if (sum(ok) > 0)
      pos / sum(ok)
    else NA_real_,
    positive_MC_low = ci[1],
    positive_MC_high = ci[2],
    shared_median = median(d$n_shared, na.rm = TRUE),
    shared_q025 = if (nrow(d) > 1L)
      unname(quantile(d$n_shared, 0.025, na.rm = TRUE))
    else NA_real_,
    shared_q975 = if (nrow(d) > 1L)
      unname(quantile(d$n_shared, 0.975, na.rm = TRUE))
    else NA_real_,
    shared_min = min(d$n_shared, na.rm = TRUE),
    shared_max = max(d$n_shared, na.rm = TRUE),
    opposite_median = median(d$n_opposite, na.rm = TRUE),
    n_LO_median = median(d$n_LO_total, na.rm = TRUE),
    n_LO_min = min(d$n_LO_total, na.rm = TRUE),
    n_LO_max = max(d$n_LO_total, na.rm = TRUE),
    n_LO_control_median = median(d$n_LO_control, na.rm = TRUE),
    n_LO_crc_median = median(d$n_LO_crc, na.rm = TRUE),
    n_EO_selected_median = median(d$n_EO, na.rm = TRUE),
    n_LO_selected_median = median(d$n_LO, na.rm = TRUE),
    median_abs_LO_shared_median = median(d$median_abs_LO_shared, na.rm = TRUE),
    median_abs_EO_shared_median = median(d$median_abs_EO_shared, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

summary <- do.call(
  rbind,
  lapply(
    f_levels,
    function(f) {
      summarise_level(draws[draws$f == f, , drop = FALSE], f, baseline$T)
    }
  )
)

ov_summary <- aggregate(cbind(n_overlap, jaccard) ~ f, overlap, function(x) median(x, na.rm = TRUE))

names(ov_summary)[-1] <- c("baseline_overlap_median", "baseline_jaccard_median")

summary <- merge(
  summary,
  ov_summary,
  by = "f",
  all.x = TRUE,
  sort = FALSE
)

summary$baseline_overlap_median[summary$f == 1] <- length(baseline_shared)

summary$baseline_jaccard_median[summary$f == 1] <- 1

summary$baseline_T <- baseline$T
summary$interpretation <- "LO downsampling sensitivity; EO fixed; empirical resampling range, not CI or P"

write.csv(
  summary,
  file.path(analysis_root, "results/downsample_summary.csv"),
  row.names = FALSE
)

freq <- as.data.frame(table(membership_df$f, membership_df$feature), stringsAsFactors = FALSE)

names(freq) <- c("f", "feature", "times_shared")

freq <- freq[freq$times_shared > 0, ]
freq$f <- as.numeric(as.character(freq$f))

write.csv(
  freq,
  file.path(analysis_root, "results/downsample_shared_frequencies.csv"),
  row.names = FALSE
)

if (length(failures)) write.csv(
  do.call(rbind, failures),
  file.path(analysis_root, "results/downsample_block_failures.csv"),
  row.names = FALSE
)

if (length(cohort_counts)) write.csv(
  do.call(rbind, cohort_counts),
  file.path(analysis_root, "results/downsample_cohort_counts.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    draws = draws,
    summary = summary,
    baseline = baseline,
    f_levels = f_levels,
    seed_base = seed_base,
    n_rep = n_rep,
    membership = membership_df,
    baseline_shared = baseline_shared,
    overlap = overlap
  ),
  file.path(analysis_root, "generated_data/downsample_results.rds")
)

print(summary)

figure_data <- read_csv(
  file.path(analysis_root, "results/downsample_summary.csv"),
  show_col_types = FALSE
) %>%
  filter(f %in% c(1, 0.75, 0.50, 0.38)) %>%
  arrange(desc(f))
write_csv(figure_data, file.path(analysis_root, "results/figure_data.csv"))
