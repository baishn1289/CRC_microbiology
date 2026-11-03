library(dplyr)
library(tidyr)
library(tibble)
library(metafor)
library(Maaslin2)
library(ggplot2)
library(grid)

# Define the site, age and interaction comparisons.
build_site_sets <- function(metadata, levels) {
  sets <- list()

  for (cl in names(levels)) {
    pairs <- combn(levels[[cl]], 2, simplify = FALSE)
    targets <- c(
      lapply(pairs, function(s) list(
        mode = "site",
        sites = s,
        age = "EO"
      )),
      lapply(pairs, function(s) list(
        mode = "site",
        sites = s,
        age = "LO"
      )),
      lapply(pairs, function(s) list(
        mode = "interaction",
        sites = s,
        age = "joint"
      )),
      lapply(levels[[cl]], function(s) list(
        mode = "age",
        sites = s,
        age = "joint"
      ))
    )

    for (tg in targets) {
      field <- paste0("Tumor_Site_", cl)
      ss <- tg$sites
      mode <- tg$mode
      md <- metadata %>%
        filter(Group == "CRC", complete_primary, .data[[field]] %in% ss) %>%
        mutate(Site = factor(.data[[field]], levels = ss))
      if (mode == "site") md <- md %>%
        filter(Age_class == tg$age)
      formula <- switch(mode,
        site = ~ Site + Age + Sex + BMI,
        age = ~ Age_class + Sex + BMI,
        interaction = ~ Age_class * Site + Sex + BMI
      )
      target <- if (length(ss) == 2) paste0(ss[2], "_vs_", ss[1]) else ss
      scope <- paste(cl, mode, tg$age, target, sep = "__")
      eligible <- character()

      for (co in unique(metadata$Cohort)) {
        z <- droplevels(md[md$Cohort == co, , drop = FALSE])
        cells <- switch(mode,
          site = table(factor(z$Site, levels = ss)),
          age = table(factor(z$Age_class, levels = c("LO", "EO"))),
          interaction = table(
            factor(z$Age_class, levels = c("LO", "EO")),
            factor(z$Site, levels = ss)
          )
        )
        x <- tryCatch(model.matrix(formula, z), error = function(e) NULL)
        rank <- if (is.null(x)) NA_integer_ else qr(x)$rank
        df <- if (is.null(x)) NA_integer_ else nrow(x) - rank
        ok <- all(cells >= 3) && !is.null(x) && rank == ncol(x) && df >= 10
        if (ok) eligible <- c(eligible, co)
      }

      if (length(eligible) >= 2) sets[[scope]] <- list(
        scope = scope,
        classification = cl,
        mode = mode,
        sites = ss,
        age = tg$age,
        formula = formula,
        metadata = md %>%
          filter(Cohort %in% eligible),
        evidence_tier = if (length(eligible) >= 3) "formal_k3" else "exploratory_k2",
        variant = "primary"
      )
    }
  }

  list(sets = sets)
}

# Apply comparison-specific filtering and the original log transform.
prepare_site_set <- function(s, abundance, metadata) {
  md <- s$metadata
  technical <- bind_rows(lapply(
    unique(md$Cohort),
    function(co) {
      x <- abundance[, md$Sample_ID[md$Cohort == co], drop = FALSE] > .001

      data.frame(
        Cohort = co,
        feature = rownames(x),
        n_samples = ncol(x),
        n_detected = rowSums(x),
        prevalence = rowMeans(x),
        pass = rowSums(x) >= 3 & rowMeans(x) >= .05
      )
    }
  ))
  features <- technical %>%
    group_by(feature) %>%
    summarise(n_pass = sum(pass), .groups = "drop") %>%
    filter(n_pass >= min(3L, n_distinct(md$Cohort))) %>%
    pull(feature) %>%
    sort()
  if (!length(features)) stop("Empty technical family: ", s$scope)
  s$features <- features
  s$technical <- technical
  s$log <- list()
  pseudo_rows <- list()

  for (co in unique(md$Cohort)) {
    ids <- metadata$Sample_ID[metadata$Cohort == co & metadata$Group %in% c("CRC", "Control")]
    x <- abundance[features, ids, drop = FALSE]
    floor <- apply(x, 1, function(v) if (any(v > 0)) min(v[v > 0]) / 2 else .0005)

    for (j in seq_len(nrow(x))) x[j, x[j, ] == 0] <- floor[j]
    s$log[[co]] <- t(log2(x))
    pseudo_rows[[co]] <- data.frame(
      Cohort = co,
      feature = features,
      pseudocount = floor
    )
  }

  s$pseudocounts <- bind_rows(pseudo_rows)

  s
}

# Pool cohort contrasts using the existing HS model.
pool_site_models <- function(local, summary) {
  groups <- split(
    seq_len(nrow(local)),
    paste(local$scope, local$contrast, local$feature, sep = "\t")
  )
  rows <- vector("list", length(groups))

  for (g in seq_along(groups)) {
    z <- local[groups[[g]], , drop = FALSE]
    b <- z[z$status == "ok" & is.finite(z$estimate) & is.finite(z$se) & z$se > 0, , drop = FALSE]
    k <- nrow(b)
    stat <- matrix(NA_real_, 1, 10)
    status <- "fewer_than_2_cohorts"
    if (k >= 2) {
      hs <- tryCatch(
        suppressWarnings(rma.uni(
          yi = b$estimate,
          sei = b$se,
          method = "HS",
          test = "z"
        )),
        error = function(e) e
      )
      fits <- list(hs)

      for (j in 1L) if (!inherits(fits[[j]], "error")) {
        f <- fits[[j]]
        stat[j, ] <- c(
          as.numeric(f$b),
          f$se,
          f$ci.lb,
          f$ci.ub,
          f$pval,
          f$tau2,
          f$I2,
          max(weights(f)),
          mean(sign(b$estimate) == sign(as.numeric(f$b))),
          f$QEp
        )
        status[j] <- if (k >= 3) "ok" else "descriptive_k2"
      } else status[j] <- paste0("meta_error:", conditionMessage(fits[[j]]))
    }

    colnames(stat) <- c(
      "estimate",
      "se",
      "ci_low",
      "ci_high",
      "p_value",
      "tau2",
      "I2",
      "max_weight_percent",
      "consistency",
      "Q_p"
    )
    rows[[g]] <- data.frame(
      scope = z$scope[1],
      configuration = "LOG_primary",
      contrast = z$contrast[1],
      feature = z$feature[1],
      method = "HS",
      inference = "z",
      k = k,
      n_samples = sum(b$n_model),
      stat,
      cohorts = paste(b$Cohort, collapse = ";"),
      status = status
    )
    if (g %% 2000L == 0L) message("Pooled models: ", g, "/", length(groups))
  }

  bind_rows(rows) %>%
    left_join(summary %>%
      select(scope, classification, variant, evidence_tier), by = "scope")
}

# Fit local contrasts using the original MaAsLin2 model.
fit_maaslin_contrasts <- function(
  metadata,
  response_matrix,
  formula,
  contrasts,
  detected_matrix,
  label,
  min_detected = 5L
) {
  mf <- model.frame(formula, metadata, na.action = na.pass)
  keep <- complete.cases(mf)
  md <- metadata[keep, , drop = FALSE]
  mf <- droplevels(mf[keep, , drop = FALSE])
  x <- model.matrix(formula, mf)
  y <- response_matrix[keep, , drop = FALSE]
  detected <- detected_matrix[keep, , drop = FALSE]
  original_features <- colnames(y)
  df <- nrow(x) - qr(x)$rank
  if (qr(x)$rank < ncol(x) || df < 10L) {
    return(bind_rows(lapply(
      names(contrasts),
      function(comparison) data.frame(
        feature = original_features,
        contrast = comparison,
        estimate = NA_real_,
        se = NA_real_,
        p_value = 1,
        ci_low = NA_real_,
        ci_high = NA_real_,
        df_residual = df,
        n_model = nrow(x),
        n_detected = colSums(detected),
        status = "rank_or_df_failure"
      )
    )))
  }

  if (any(!is.finite(y))) stop("Non-finite response matrix.")
  feature_map <- data.frame(
    package_id = sprintf("F%04d", seq_len(ncol(y))),
    feature = original_features
  )
  term_map <- data.frame(
    package_term = c("(Intercept)", sprintf("D%02d", seq_len(ncol(x) - 1L))),
    original_term = colnames(x)
  )
  package_data <- as.data.frame(y)
  names(package_data) <- feature_map$package_id
  package_metadata <- as.data.frame(x[, -1, drop = FALSE])
  names(package_metadata) <- term_map$package_term[-1]
  rownames(package_data) <- rownames(package_metadata) <- md$Sample_ID

  abundance_cutoff <- min(y) - 1
  model_folder <- file.path(models_dir, gsub("[^A-Za-z0-9_-]", "_", label))
  model_file <- file.path(model_folder, "fits", "models.rds")
  {
    dir.create(model_folder, recursive = TRUE, showWarnings = FALSE)
    write.csv(
      feature_map,
      file.path(model_folder, "feature_mapping.csv"),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
    write.csv(
      term_map,
      file.path(model_folder, "design_terms.csv"),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
    write.csv(
      data.frame(Sample_ID = md$Sample_ID, package_metadata),
      file.path(model_folder, "input_design.csv"),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
    capture.output(
      package_fit <- Maaslin2::Maaslin2(
        input_data = package_data,
        input_metadata = package_metadata,
        output = model_folder,
        analysis_method = "LM",
        normalization = "NONE",
        transform = "NONE",
        standardize = FALSE,
        min_abundance = abundance_cutoff,
        min_prevalence = 0,
        min_variance = 0,
        fixed_effects = names(package_metadata),
        random_effects = NULL,
        correction = "BH",
        max_significance = 1,
        cores = 1,
        save_models = TRUE,
        plot_heatmap = FALSE,
        plot_scatter = FALSE,
        max_pngs = 0
      ),
      file = file.path(model_folder, "console.txt")
    )
    rm(package_fit)
    message("MaAsLin2 fitted: ", label, " / n=", nrow(y), " / features=", ncol(y))
  }

  fits <- readRDS(model_file)

  transformed <- read.delim(
    file.path(model_folder, "features", "filtered_data_norm_transformed.tsv"),
    check.names = FALSE,
    colClasses = c(feature = "character")
  )
  retained <- intersect(feature_map$package_id, names(transformed))
  transformed_matrix <- as.matrix(transformed[match(
    md$Sample_ID,
    transformed[[1]]
  ), retained, drop = FALSE])
  expected_matrix <- as.matrix(package_data[, retained, drop = FALSE])
  input_error <- max(abs(transformed_matrix - expected_matrix))
  if (!is.finite(input_error) ||
    input_error > 1e-10)
    stop("MaAsLin2 changed the response or sample order.")
  rows <- list()
  index <- 1L

  for (f in seq_along(original_features)) {
    fit <- fits[[feature_map$package_id[f]]]
    n_detected <- sum(detected[, f])

    for (comparison in names(contrasts)) {
      row <- data.frame(
        feature = original_features[f],
        contrast = comparison,
        estimate = NA_real_,
        se = NA_real_,
        p_value = 1,
        ci_low = NA_real_,
        ci_high = NA_real_,
        df_residual = df,
        n_model = nrow(y),
        n_detected = n_detected,
        status = ifelse(n_detected < min_detected, "fewer_than_5_detected", "invalid_variance")
      )
      if (inherits(fit, "glm")) {
        actual_x <- model.matrix(fit)
        if (ncol(actual_x) != ncol(x) || max(abs(unname(actual_x) - unname(x))) > 1e-10)
          stop("MaAsLin2 changed the design matrix.")
        weights <- setNames(rep(0, ncol(x)), colnames(x))
        if (!all(names(contrasts[[comparison]]) %in% names(weights))) stop("Missing contrast term.")
        weights[names(contrasts[[comparison]])] <- contrasts[[comparison]]
        covariance <- vcov(fit)
        estimate <- sum(weights * coef(fit))
        variance <- as.numeric(t(weights) %*% covariance %*% weights)
        if (n_detected >= min_detected && is.finite(variance) && variance > 1e-16) {
          row$estimate <- estimate
          row$se <- sqrt(variance)
          row$df_residual <- df.residual(fit)
          row$p_value <- 2 * pt(-abs(estimate / row$se), row$df_residual)
          row$ci_low <- estimate - qt(0.975, row$df_residual) * row$se
          row$ci_high <- estimate + qt(0.975, row$df_residual) * row$se
          row$status <- "ok"
        }
      }

      rows[[index]] <- row
      index <- index + 1L
    }
  }

  bind_rows(rows)
}

# Apply the final cohort-support, heterogeneity and weight rules.
apply_meta_qc <- function(meta, local, main_qc, weight_limit = 50) {
  stopifnot(!anyDuplicated(meta[, c("scope", "feature")]))
  decisions <- list()

  for (sc in unique(meta$scope)) {
    m <- meta[meta$scope == sc, , drop = FALSE]
    l <- local[local$scope == sc & local$status == "ok" & is.finite(local$estimate) &
      is.finite(local$se) & local$se > 0, , drop = FALSE]
    object <- list(
      meta_fits = data.frame(
        coef = m$estimate,
        I2 = m$I2,
        weight_max = m$max_weight_percent,
        row.names = m$feature
      ),
      maaslin_fits = list(data.frame(
        feature = l$feature,
        coef = l$estimate,
        batch = l$Cohort
      ))
    )

    z <- suppressMessages(main_qc(
      object,
      i2_tier1 = 50,
      i2_tier2 = 80,
      weight_limit = weight_limit,
      base_consistency = .60,
      rescue_consistency = .75
    ))
    decisions[[sc]] <- z %>%
      transmute(
        scope = sc,
        feature = Pathway,
        qc_nonzero_cohorts = n_total_cohorts,
        qc_same_direction = n_same_dir,
        qc_consistency_score = consistency_score,
        qc_source_pass = pass_qc,
        qc_source_status = qc_status
      )
  }

  support <- local %>%
    filter(status == "ok", is.finite(estimate), estimate != 0, is.finite(se), se > 0) %>%
    group_by(scope, feature) %>%
    summarise(qc_nonzero_actual = n_distinct(Cohort), .groups = "drop")
  out <- meta %>%
    left_join(bind_rows(decisions), by = c("scope", "feature")) %>%
    left_join(support, by = c("scope", "feature")) %>%
    mutate(
      qc_nonzero_cohorts = coalesce(qc_nonzero_cohorts, qc_nonzero_actual, 0L),
      qc_source_pass = coalesce(qc_source_pass, FALSE),
      qc_estimable = status %in% c("ok", "descriptive_k2") &
        is.finite(estimate) &
        is.finite(se) &
        se > 0 &
        is.finite(p_value) & p_value >= 0 & p_value <= 1,
      qc_enough_cohorts = qc_nonzero_cohorts >= 3L,
      qc_weight_ok = is.finite(max_weight_percent) & max_weight_percent < weight_limit,
      qc_I2_ok = is.finite(I2) & I2 >= 0 & I2 < 80,
      qc_direction_ok = case_when(
        is.finite(I2) & I2 < 50 & !is.na(qc_consistency_score) ~ qc_consistency_score >= .60,
        is.finite(I2) &
          I2 >= 50 &
          I2 < 80 &
          !is.na(qc_consistency_score) ~ qc_consistency_score >= .75,
        TRUE ~ FALSE
      ),
      qc_pass = qc_estimable &
        qc_source_pass &
        qc_enough_cohorts &
        qc_weight_ok &
        qc_I2_ok &
        qc_direction_ok
    )

  flags <- c("qc_estimable", "qc_enough_cohorts", "qc_weight_ok", "qc_I2_ok", "qc_direction_ok")
  reasons <- c(
    "not_estimable",
    "fewer_than_3_nonzero_cohorts",
    paste0("max_weight_ge", weight_limit, "_or_invalid"),
    "I2_ge80_or_invalid",
    "insufficient_direction_consistency"
  )
  out$qc_failure_reasons <- apply(
    as.matrix(out[, flags]),
    1,
    function(v) paste(reasons[!v], collapse = ";")
  )
  out
}

# Adjust P values within the estimable family for each comparison.
estimable_bh <- function(d) {
  stopifnot(!anyDuplicated(d[, c("scope", "feature")]), !anyNA(d$qc_pass_50))
  strict <- d$evidence_tier == "formal_k3"
  explore <- d$evidence_tier == "exploratory_k2" & grepl("__site__EO__", d$scope, fixed = TRUE)
  d$policy_tier <- ifelse(
    strict,
    "multi_ge3_fullQC50",
    ifelse(explore, "EO_two_cohort_exploratory", "excluded_two_cohort_interaction")
  )
  valid <- d$status %in% c("ok", "descriptive_k2") & is.finite(d$estimate) &
    is.finite(d$se) & d$se > 0 & is.finite(d$p_value) & d$p_value >= 0 & d$p_value <= 1
  stopifnot(all(valid[strict & d$qc_pass_50]), all(d$k[strict & d$qc_pass_50] >= 3L))

  d$in_BH_family <- (strict & d$qc_pass_50) | (explore & valid & d$k == 2L)
  d$reportable <- d$in_BH_family
  d$P_for_BH <- ifelse(d$in_BH_family, d$p_value, NA_real_)
  d$q_comparison <- rep(NA_real_, nrow(d))
  d$n_comparison <- rep(0L, nrow(d))

  for (sc in unique(d$scope)) {
    rows <- which(d$scope == sc)
    members <- rows[d$in_BH_family[rows]]
    d$n_comparison[rows] <- length(members)
    d$q_comparison[members] <- p.adjust(d$P_for_BH[members], "BH")
  }

  d
}

plotcase_forest <- function(df_raw, x_limits, breaks, y_limits) {
  plot_df <- df_raw
  data_df <- plot_df %>%
    filter(kind != "header")
  ggplot(plot_df, aes(y = row_id)) +
    geom_vline(
      xintercept = 0,
      linetype = 2,
      linewidth = .3,
      color = "grey65"
    ) +
    geom_errorbar(
      data = data_df,
      aes(xmin = ci_low, xmax = ci_high),
      orientation = "y",
      width = .12,
      linewidth = .35,
      color = "grey30"
    ) +
    {
      if (any(data_df$kind == "cohort")) geom_point(
        data = data_df %>%
          filter(kind == "cohort"),
        aes(x = estimate, size = weight),
        shape = 22,
        fill = "grey55",
        color = "grey20",
        stroke = .2
      )
    } +
    geom_point(
      data = data_df %>%
        filter(kind == "pooled"),
      aes(x = estimate),
      shape = 18,
      size = 2.1,
      color = "black"
    ) +
    scale_size_area(max_size = 2.8, guide = "none") +
    scale_x_continuous(
      limits = x_limits,
      breaks = breaks,
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    labs(x = expression("Beta" ~ "(log"[2] ~ "RA)"), y = NULL) +
    theme_classic(base_size = 7, base_family = "sans") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(color = "black", size = 6.8),
      axis.title.x = element_text(size = 7),
      axis.line.x = element_line(linewidth = .3),
      axis.ticks.x = element_line(linewidth = .3),
      plot.margin = margin(0, 1, 0, 1, unit = "mm")
    )
}

label_column <- function(d, y_limits, text_size = 2.4) {
  ggplot(d, aes(x = .01, y = row_id)) +
    geom_text(
      aes(label = Level, fontface = ifelse(kind == "pooled", "bold", "plain")),
      hjust = 0,
      size = text_size
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

panel_title <- function(
  tag,
  title,
  x,
  y,
  italic = FALSE,
  size = 8
) {
  list(
    textGrob(
      tag,
      x = x,
      y = y,
      just = c("left", "top"),
      gp = gpar(fontsize = 9, fontface = "bold")
    ),
    textGrob(
      title,
      x = x + .021,
      y = y,
      just = c("left", "top"),
      gp = gpar(fontsize = size, fontface = if (italic) "italic" else "plain")
    )
  )
}
