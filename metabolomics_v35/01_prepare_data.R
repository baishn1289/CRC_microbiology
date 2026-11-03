library(data.table)

code_dir <- file.path(
  base_symbol,
  "downstream_analysis/src/github_code_availability/submission_version_261102/metabolomics_v35"
)
output_dir <- file.path(
  base_symbol,
  "downstream_analysis/output/submission_version_261102/metabolomics_v35"
)
analysis_dir <- file.path(base_symbol, "FUSCC_metabonomics_analysis")
v35_dir <- output_dir

for (folder in c("objects", "tables", "plots")) dir.create(
  file.path(output_dir, folder),
  recursive = TRUE,
  showWarnings = FALSE
)

# Frozen peak audit contains the original raw peaks and 557-participant metadata.
inputs <- readRDS(file.path(analysis_dir, "v30/peak_audit_inputs_v30.rds"))
meta <- as.data.table(inputs$meta)
raw_all <- as.data.table(inputs$raw_all)
qc_cols <- grep("^QC[0-9]+$", names(raw_all), value = TRUE)
stopifnot(nrow(meta) == 557L, length(qc_cols) == 62L, !anyDuplicated(meta$Sample_ID))
meta[, Disease := factor(Disease, levels = c("Control", "CRC"))]
meta[, Age_group := factor(Age_group, levels = c("LO", "EO"))]
meta[, Sex := factor(Sex, levels = c("female", "male"))]
members <- list(
  U_P = c("Arginine", "Citrulline", "Ornithine_related", "N-Acetylornithine", "Proline"),
  D0 = c("5-Aminopentanoate", "N-Acetylputrescine", "N-Carbamoylputrescine")
)
member_order <- unlist(members, use.names = FALSE)
mapping <- fread(file.path(analysis_dir, "v30/tables/frozen_peak_mapping_v30.csv"))[
  variant == "representative_core" & canonical_metabolite %in% member_order
]
mapping <- mapping[match(member_order, canonical_metabolite)]
stopifnot(nrow(mapping) == 8L, !anyDuplicated(mapping$peak_id))
decisions <- fread(file.path(
  analysis_dir,
  "extended_metabonomics_results_v29/tables/targeted_mode_dedup_decisions.csv"
))

# Match the original selected rows using annotation, mode and all sample intensities.
selected <- decisions[selected_mode == TRUE]
selected[, peak_id := NA_character_]
selected[, mapping_method := "name_mode_and_557_intensities"]

for (i in seq_len(nrow(selected))) {
  candidates <- raw_all[metabolite_identification == selected$metabolite_identification[i] &
    mode == selected$mode[i]]
  target <- as.numeric(selected[i, meta$Sample_ID, with = FALSE])
  if (!nrow(candidates)) {
    mode_rows <- raw_all[mode == selected$mode[i]]
    first_values <- mode_rows[[meta$Sample_ID[1]]]
    close_first <- if (is.na(target[1])) is.na(first_values) else
      is.finite(first_values) & abs(first_values - target[1]) <= 1e-8 * max(1, abs(target[1]))
    candidates <- mode_rows[close_first]
    selected$mapping_method[i] <- "historical_label_encoding_fallback_full_intensity_match"
  }

  same <- vapply(
    seq_len(nrow(candidates)),
    function(j)
      isTRUE(all.equal(
        as.numeric(candidates[j, meta$Sample_ID, with = FALSE]),
        target,
        tolerance = 1e-10
      )),
    logical(1)
  )
  if (sum(same) != 1L) stop(
    "Cannot map selected row ",
    i,
    ": ",
    selected$metabolite_identification[i],
    " mode=",
    selected$mode[i],
    " matching peaks=",
    sum(same),
    " candidates=",
    nrow(candidates)
  )
  selected$peak_id[i] <- candidates$peak_id[same]
}

stopifnot(nrow(selected) == 2014L, !anyDuplicated(selected$peak_id))
selected_raw <- raw_all[match(selected$peak_id, peak_id)]
intensity <- as.matrix(selected_raw[, meta$Sample_ID, with = FALSE])
rownames(intensity) <- selected$peak_id
detected <- is.finite(intensity) & intensity > 0
group_detection <- sapply(
  unique(meta$Group4),
  function(g) rowMeans(detected[, meta$Group4 == g, drop = FALSE])
)
qc_mat <- as.matrix(selected_raw[, qc_cols, with = FALSE])
qc_cv <- 100 * apply(qc_mat, 1, sd, na.rm = TRUE) / rowMeans(qc_mat, na.rm = TRUE)
feature_qc <- selected_raw[, .(
  peak_id,
  metabolite_identification,
  mode,
  retention_time,
  mass_to_charge,
  database_identifier,
  chemical_formula
)]
feature_qc[, `:=`(
  v29_label = selected$metabolite_identification,
  mapping_method = selected$mapping_method
)]
feature_qc[, `:=`(
  missing_fraction = rowMeans(is.na(intensity)),
  detection_fraction = rowMeans(detected),
  max_group_detection = apply(group_detection, 1, max),
  qc_cv_percent = qc_cv,
  qc_n_positive = rowSums(is.finite(qc_mat) & qc_mat > 0)
)]
feature_qc[, keep_primary := missing_fraction <= 0.5 & max_group_detection >= 0.7]
stopifnot(
  sum(feature_qc$keep_primary) == 2013L,
  all(mapping$peak_id %in% feature_qc[keep_primary == TRUE, peak_id])
)

fwrite(feature_qc, file.path(output_dir, "tables/feature_qc.csv"))

# Use the same eligible reference peaks; unused exploratory peaks are omitted.
all_ids <- unique(c(feature_qc[keep_primary == TRUE, peak_id], mapping$peak_id))
all_raw <- raw_all[match(all_ids, peak_id)]
positive_mat <- as.matrix(all_raw[, meta$Sample_ID, with = FALSE])
rownames(positive_mat) <- all_ids

for (i in seq_len(nrow(positive_mat))) {
  x <- positive_mat[i, ]
  positive_min <- min(x[is.finite(x) & x > 0])
  bad <- !is.finite(x) | x <= 0
  positive_mat[i, bad] <- positive_min / 2
}

stopifnot(all(is.finite(positive_mat)), all(positive_mat > 0))
log_pqn <- log2(positive_mat)
factor_rows <- list()

for (mode_name in c("positive", "negative")) {
  reference_ids <- feature_qc[keep_primary == TRUE & mode == mode_name, peak_id]
  reference <- apply(positive_mat[reference_ids, , drop = FALSE], 1, median)
  factors <- apply(sweep(positive_mat[reference_ids, , drop = FALSE], 1, reference, "/"), 2, median)
  mode_ids <- all_raw[mode == mode_name, peak_id]
  log_pqn[mode_ids, ] <- log2(sweep(positive_mat[mode_ids, , drop = FALSE], 2, factors, "/"))
  factor_rows[[mode_name]] <- data.table(
    Sample_ID = meta$Sample_ID,
    mode = mode_name,
    pqn_factor = factors
  )
}

fwrite(rbindlist(factor_rows), file.path(output_dir, "tables/PQN_sample_factors.csv"))

# Peak z scores use all 557 samples; canonical and panel means are not rescaled.
peak_z <- t(scale(t(log_pqn)))
canonical <- peak_z[mapping$peak_id, , drop = FALSE]
rownames(canonical) <- mapping$canonical_metabolite
score <- t(vapply(members, function(m) colMeans(canonical[m, , drop = FALSE]), numeric(nrow(meta))))
colnames(score) <- meta$Sample_ID
member_log2 <- log_pqn[mapping$peak_id, meta$Sample_ID, drop = FALSE]
rownames(member_log2) <- mapping$canonical_metabolite
background_z <- peak_z[feature_qc[keep_primary == TRUE, peak_id], , drop = FALSE]
set_ids <- lapply(members, function(m) mapping[match(m, canonical_metabolite), peak_id])

# Separate covariate slopes for the two age strata.
strata <- c("EO", "LO")
x_parts <- lapply(
  strata,
  function(g) model.matrix(~ Disease + Age + Sex + BMI, data = meta[Age_group == g])
)
names(x_parts) <- strata
design <- matrix(0, nrow(meta), sum(vapply(x_parts, ncol, integer(1))))
rownames(design) <- meta$Sample_ID
colnames(design) <- unlist(lapply(
  strata,
  function(g) paste(g, make.names(colnames(x_parts[[g]])), sep = "_")
))

for (g in strata) {
  columns <- grep(paste0("^", g, "_"), colnames(design))
  design[meta$Age_group == g, columns] <- x_parts[[g]]
}

stopifnot(qr(design)$rank == ncol(design), ncol(design) == 10L)
contrast_mat <- matrix(
  0,
  ncol(design),
  3L,
  dimnames = list(colnames(design), c("EO", "LO", "EO_minus_LO"))
)
contrast_mat["EO_DiseaseCRC", c("EO", "EO_minus_LO")] <- 1
contrast_mat["LO_DiseaseCRC", "LO"] <- 1
contrast_mat["LO_DiseaseCRC", "EO_minus_LO"] <- -1

# Clinical coding matches the original complete-case designs.
crc <- copy(meta[Disease == "CRC"])
crc[, Stage_label := trimws(as.character(Stage_raw))]
crc[, Stage_status := fifelse(
  is.na(Stage_label) | Stage_label == "",
  "Missing",
  fifelse(
    tolower(Stage_label) == "pcr",
    "Recorded_pCR",
    fifelse(Stage_label %in% as.character(0:4), "Known_TNM", "Invalid")
  )
)]
stopifnot(!any(crc$Stage_status == "Invalid"))
crc[, Stage_factor := factor(
  ifelse(Stage_status == "Known_TNM", paste0("Stage", Stage_label), NA_character_),
  levels = paste0("Stage", 0:4)
)]
crc[, Tumor_Site_factor := factor(
  fcase(
    as.character(Tumor_Site_raw) == "0",
    "right_colon",
    as.character(Tumor_Site_raw) == "1",
    "left_colon",
    as.character(Tumor_Site_raw) == "2",
    "rectum",
    default = NA_character_
  ),
  levels = c("right_colon", "left_colon", "rectum")
)]
crc[, Stage_binary := factor(
  fcase(
    Stage_label %in% c("0", "1", "2"),
    "Stage0_II",
    Stage_label %in% c("3", "4"),
    "StageIII_IV",
    default = NA_character_
  ),
  levels = c("Stage0_II", "StageIII_IV")
)]
crc[, Tumor_Site_binary := factor(
  fcase(
    as.character(Tumor_Site_raw) %in% c("0", "1"),
    "colon",
    as.character(Tumor_Site_raw) == "2",
    "rectum",
    default = NA_character_
  ),
  levels = c("colon", "rectum")
)]
crc[, clinical_complete := complete.cases(Stage_factor, Tumor_Site_factor, Sex, BMI)]
crc[, binary_complete := complete.cases(Age_group, Sex, BMI, Stage_binary, Tumor_Site_binary)]
stopifnot(identical(crc$clinical_complete, crc$binary_complete))
d <- copy(crc[clinical_complete == TRUE])
stopifnot(nrow(d) == 302L, sum(d$Age_group == "EO") == 139L, sum(d$Age_group == "LO") == 163L)
crc_formulas <- c(
  clinical_multicategory = "~ Age_group + Sex + BMI + Stage_factor + Tumor_Site_factor",
  clinical_binary = "~ Age_group + Sex + BMI + Stage_binary + Tumor_Site_binary"
)
crc_designs <- lapply(crc_formulas, function(formula) model.matrix(as.formula(formula), data = d))

for (model_name in names(crc_designs)) rownames(crc_designs[[model_name]]) <- d$Sample_ID

ctx <- list(
  meta = meta[, .(Sample_ID, Patient_ID, Group4, Age_group, Disease, Age, Sex, BMI)],
  members = members,
  mapping = mapping,
  canonical_z = canonical,
  member_log2 = member_log2,
  scores = score,
  background_z = background_z,
  set_ids = set_ids,
  design = design,
  contrasts = contrast_mat
)
prep <- list(feature_qc = feature_qc, log_mats = list(PQN_log2 = log_pqn))
saveRDS(
  list(
    ctx = ctx,
    prep = prep,
    crc_designs = crc_designs,
    crc_formulas = crc_formulas
  ),
  file.path(output_dir, "objects/prepared_data.rds")
)
saveRDS(ctx, file.path(output_dir, "objects/fixed_scheme_context.rds"))
fwrite(mapping, file.path(output_dir, "tables/fixed_member_peak_mapping.csv"))
