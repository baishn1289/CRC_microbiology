library(tidyverse)
library(vegan)
library(permute)

# Fit raw-TSS models with cohort fixed effects and within-cohort permutations.
rows <- list()
tests <- list()

for (age in c("EO", "LO")) {
  p <- primary[[age]]
  m <- droplevels(p$metadata)
  d <- vegan::vegdist(t(p$feature_raw), method = "bray")
  x <- model.matrix(~ Cohort + Age + Sex + BMI + Group, m)
  stopifnot(qr(x)$rank == ncol(x), nrow(m) > ncol(x))
  model <- vegan::dbrda(d ~ Cohort + Age + Sex + BMI + Group, data = m)
  set.seed(123)
  test <- anova(
    model,
    by = "margin",
    scope = "Group",
    permutations = permute::how(nperm = 999, blocks = m$Cohort)
  )
  stopifnot("Group" %in% rownames(test))
  ss_column <- intersect(c("SumOfSqs", "Variance"), colnames(test))[1]
  rows[[age]] <- data.frame(
    Age_class = age,
    analysis = "raw_TSS_cohort_fixed_all_samples",
    n = nrow(m),
    n_control = sum(m$Group == "Control"),
    n_crc = sum(m$Group == "CRC"),
    n_cohorts = nlevels(m$Cohort),
    Group_R2 = test["Group", ss_column] / model$tot.chi,
    Group_F = test["Group", "F"],
    Group_P = test["Group", "Pr(>F)"],
    residual_df = test["Residual", "Df"],
    permutations = 999,
    restriction = "within_Cohort",
    formula = "Cohort + Age + Sex + BMI + Group"
  )
  tests[[age]] <- list(model = model, test = test)
  message("FRESH RAW COHORT SENSITIVITY COMPLETE: ", age)
}

raw_cohort_sensitivity <- bind_rows(rows)
write.csv(
  raw_cohort_sensitivity,
  file.path(output_dir, "raw_cohort_sensitivity.csv"),
  row.names = FALSE
)
