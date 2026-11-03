library(dplyr)
library(metafor)

# Build eligible model metadata and its fixed-design projection terms.
prepare_design <- function(dat, target) {
  keep <- dat$metadata$complete
  if (target != "Common_age") keep <- keep & dat$metadata$Group == sub("_age$", "", target)
  mm <- droplevels(dat$metadata[keep, , drop = FALSE])

  # Set the covariates and fixed-design projection terms.
  rhs <- if (target == "Common_age") ~ D + Age10 + Sex + BMI_c else ~ Age10 + Sex + BMI_c
  if (length(unique(mm$Sex)) < 2) rhs <- update(rhs, ~ . - Sex)
  if (sd(mm$BMI_c) == 0) rhs <- update(rhs, ~ . - BMI_c)
  X <- model.matrix(rhs, mm)
  rank_ok <- qr(X)$rank == ncol(X)
  influence <- if (rank_ok) solve(crossprod(X)) %*% t(X) else NULL
  h <- if (rank_ok) colSums(t(X) * influence) else rep(NA_real_, nrow(X))
  age_groups <- if (target == "Common_age") split(mm$Age, mm$Group) else list(mm$Age)

  # Check cohort eligibility before fitting individual species.
  reasons <- c(
    group_n_lt10 = min(lengths(age_groups)) < 10,
    unique_ages_lt5 = min(vapply(age_groups, function(a) length(unique(a)), integer(1))) < 5,
    age_span_lt15 = min(vapply(age_groups, function(a) diff(range(a)), numeric(1))) < 15,
    rank_failure = !rank_ok,
    residual_df_lt10 = nrow(X) - ncol(X) < 10,
    leverage_one = !rank_ok || max(h) >= 1 - 1e-8
  )
  audit <- data.frame(
    Cohort = unique(mm$Cohort),
    Target = target,
    n = nrow(mm),
    n_control = sum(mm$Group == "Control"),
    n_crc = sum(mm$Group == "CRC"),
    parameters = ncol(X),
    residual_df = nrow(X) - ncol(X),
    max_leverage = if (rank_ok) max(h) else NA_real_,
    eligible = !any(reasons),
    reasons = paste(names(reasons)[reasons], collapse = ";"),
    formula = paste(deparse(rhs), collapse = "")
  )

  list(
    metadata = mm,
    keep = keep,
    X = X,
    influence = influence,
    leverage = h,
    audit = audit
  )
}

# Fit cohort age slopes and calculate HC3 variances.
fit_target <- function(dat, obj, target, features, pseudocount) {
  mm <- obj$metadata
  cohort <- unique(mm$Cohort)
  audit <- data.frame(
    feature = features,
    Cohort = cohort,
    Target = target,
    eligible = FALSE,
    reason = "absent_in_cohort",
    residual_df = obj$audit$residual_df,
    detected_control = NA_integer_,
    detected_crc = NA_integer_,
    variance_control = NA_real_,
    variance_crc = NA_real_
  )
  if (!obj$audit$eligible) {
    audit$reason <- paste0("cohort:", obj$audit$reasons)
    return(list(local = data.frame(), audit = audit))
  }

  present <- intersect(features, rownames(dat$raw))
  xx <- dat$raw[present, obj$keep, drop = FALSE]
  Y <- t(xx)

  for (j in seq_along(present)) Y[Y[, j] == 0, j] <- pseudocount[[present[j]]]
  Y <- log2(Y)

  # Apply the original detection and response-variance thresholds.
  detected <- rowSums(xx > .001)
  dc <- if (any(mm$Group == "Control"))
    rowSums(xx[, mm$Group == "Control", drop = FALSE] > .001) else
    rep(0L, length(present))
  dr <- if (any(mm$Group == "CRC"))
    rowSums(xx[, mm$Group == "CRC", drop = FALSE] > .001) else
    rep(0L, length(present))
  vc <- if (any(mm$Group == "Control"))
    apply(Y[mm$Group == "Control", , drop = FALSE], 2, var) else
    rep(NA_real_, length(present))
  vr <- if (any(mm$Group == "CRC"))
    apply(Y[mm$Group == "CRC", , drop = FALSE], 2, var) else
    rep(NA_real_, length(present))
  eligible <- if (target == "Common_age") dc >= 5 & dr >= 5 & vc > 1e-10 & vr > 1e-10 else
    detected >= 5 & detected / nrow(mm) >= .05 & apply(Y, 2, var) > 1e-12
  index <- match(present, audit$feature)
  audit$detected_control[index] <- dc
  audit$detected_crc[index] <- dr
  audit$variance_control[index] <- vc
  audit$variance_crc[index] <- vr
  audit$reason[index] <- ifelse(
    eligible,
    "pending_variance_check",
    if (target == "Common_age")
      "either_group_detection_lt5_or_variance_le1e-10" else
      "group_detection_lt5_or_prevalence_lt005_or_constant"
  )

  # Fit all retained responses against the same design.
  selected <- which(eligible)
  YY <- Y[, selected, drop = FALSE]
  B <- obj$influence %*% YY
  E <- YY - obj$X %*% B
  local <- list()

  for (j in seq_along(selected)) {
    feature <- present[selected[j]]
    a <- match(feature, audit$feature)

    # HC3 uses leverage-adjusted residuals for the coefficient covariance.
    Z <- sweep(obj$influence, 2, E[, j] / (1 - obj$leverage), "*")
    V <- tcrossprod(Z)
    variance <- V["Age10", "Age10"]
    if (!is.finite(variance) || variance <= 1e-12) {
      audit$reason[a] <- "age_coefficient_variance_le1e-12_or_invalid"
      next
    }

    audit$eligible[a] <- TRUE
    audit$reason[a] <- "ok"
    pos <- xx[feature, ] > .001
    local[[length(local) + 1L]] <- data.frame(
      feature = feature,
      Cohort = cohort,
      Target = target,
      beta = unname(B["Age10", j]),
      variance = variance,
      n = nrow(mm),
      n_control = sum(mm$Group == "Control"),
      n_crc = sum(mm$Group == "CRC"),
      residual_df = obj$audit$residual_df,
      residual_df_control = if (target == "Control_age") obj$audit$residual_df else NA_real_,
      residual_df_crc = if (target == "CRC_age") obj$audit$residual_df else NA_real_,
      detected_control = dc[selected[j]],
      detected_crc = dr[selected[j]],
      detected = sum(pos),
      prevalence = mean(pos),
      age_min = min(mm$Age),
      age_max = max(mm$Age),
      unique_ages = length(unique(mm$Age)),
      detected_age_span = diff(range(mm$Age[pos])),
      detected_unique_ages = length(unique(mm$Age[pos])),
      max_leverage = max(obj$leverage),
      formula = obj$audit$formula
    )
  }

  list(local = bind_rows(local), audit = audit)
}

# Pool cohort slopes with REML, standard KH inference and meta-QC.
pool_target <- function(dd, feature, target) {
  dd <- dd[is.finite(dd$beta) & is.finite(dd$variance) & dd$variance > 1e-12, , drop = FALSE]
  row <- data.frame(
    feature = feature,
    Species = gsub("_", " ", sub(".*\\|s__", "", feature)),
    Target = target,
    analysis_version = "v7_unified_df10_20261023",
    method = "standard_KH",
    status = "fewer_than_3_cohorts",
    k = nrow(dd),
    n = sum(dd$n),
    n_control = sum(dd$n_control),
    n_crc = sum(dd$n_crc),
    beta = NA_real_,
    se = NA_real_,
    ci_lo = NA_real_,
    ci_hi = NA_real_,
    p = NA_real_,
    meta_df = if (nrow(dd)) nrow(dd) - 1 else NA_integer_,
    regression_df_threshold = 10L,
    min_local_residual_df = if (nrow(dd)) min(dd$residual_df) else NA_real_,
    I2 = NA_real_,
    tau2 = NA_real_,
    kh_scale = NA_real_,
    max_weight_pct = NA_real_,
    nonzero_k = NA_integer_,
    direction_consistency = NA_real_,
    qc_pass = FALSE,
    qc_reasons = "not_estimable",
    cohorts = paste(dd$Cohort, collapse = ";"),
    warning = ""
  )
  if (nrow(dd) < 3) return(list(row = row, weights = data.frame()))
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      rma.uni(
        yi = dd$beta,
        vi = dd$variance,
        method = "REML",
        test = "knha",
        control = list(
          threshold = 1e-12,
          maxiter = 2000,
          stepadj = .5
        )
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    row$status <- "meta_error"
    row$warning <- conditionMessage(fit)
    return(list(row = row, weights = data.frame()))
  }


  # Summarise pooling weights, direction agreement and meta-QC.
  w <- 1 / (dd$variance + fit$tau2)
  beta <- as.numeric(fit$b)
  nz <- dd$beta != 0
  consistency <- mean(sign(dd$beta[nz]) == sign(beta))
  reasons <- c(
    nonzero_cohorts_lt3 = sum(nz) < 3,
    max_weight_ge40 = max(w / sum(w)) >= .4,
    I2_ge80 = fit$I2 >= 80,
    inconsistent_direction = !is.finite(consistency) || consistency < if (fit$I2 < 50) .60 else .75
  )
  row$status <- "ok"
  row$beta <- beta
  row$se <- fit$se
  row$ci_lo <- fit$ci.lb
  row$ci_hi <- fit$ci.ub
  row$p <- fit$pval
  row$I2 <- fit$I2
  row$tau2 <- fit$tau2
  row$kh_scale <- sum(w * (dd$beta - beta)^2) / (nrow(dd) - 1)
  row$max_weight_pct <- 100 * max(w / sum(w))
  row$nonzero_k <- sum(nz)
  row$direction_consistency <- consistency
  row$qc_pass <- !any(reasons)
  row$qc_reasons <- paste(names(reasons)[reasons], collapse = ";")
  row$warning <- paste(unique(warnings), collapse = ";")

  list(
    row = row,
    weights = data.frame(
      feature = feature,
      Target = target,
      Cohort = dd$Cohort,
      weight_pct = 100 * w / sum(w)
    )
  )
}
