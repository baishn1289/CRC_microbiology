library(dplyr)
library(Maaslin2)
library(metafor)

# Fit every current scope and cohort once; no historical model reuse.
models_dir <- file.path(output_dir, "models")
dir.create(models_dir, showWarnings = FALSE)
local_rows <- list()

for (sc in names(sets)) {
  s <- sets[[sc]]
  md <- s$metadata

  for (co in unique(md$Cohort)) {
    d <- md %>%
      filter(Cohort == co)
    coefficient <- switch(s$mode,
      site = paste0("Site", s$sites[2]),
      age = "Age_classEO",
      interaction = paste0("Age_classEO:Site", s$sites[2])
    )
    contrast <- switch(s$mode,
      site = paste(s$age, s$sites[2], "minus", s$sites[1], sep = "_"),
      age = paste0(s$sites, "_EO_minus_LO"),
      interaction = paste("Age_interaction", s$sites[2], "minus", s$sites[1], sep = "_")
    )
    detected <- t(abundance_full[s$features, d$Sample_ID, drop = FALSE] > .001)
    z <- fit_maaslin_contrasts(
      d,
      s$log[[co]][d$Sample_ID, , drop = FALSE],
      s$formula,
      setNames(list(setNames(1, coefficient)), contrast),
      detected,
      paste(sc, co, sep = "__")
    )
    z$scope <- sc
    z$configuration <- "LOG_primary"
    z$Cohort <- co
    z$classification <- s$classification
    z$variant <- "primary"
    z$evidence_tier <- s$evidence_tier
    local_rows[[paste(sc, co)]] <- z
  }
}

direct_local <- bind_rows(local_rows) %>%
  mutate(taxon = gsub("_", " ", sub(".*\\|s__", "", feature)))
model_summary <- bind_rows(lapply(
  sets,
  function(s) data.frame(
    scope = s$scope,
    classification = s$classification,
    variant = "primary",
    evidence_tier = s$evidence_tier
  )
))
direct_meta <- pool_site_models(direct_local, model_summary) %>%
  mutate(taxon = gsub("_", " ", sub(".*\\|s__", "", feature)))
write.csv(direct_local, file.path(output_dir, "local_models.csv"), row.names = FALSE)
