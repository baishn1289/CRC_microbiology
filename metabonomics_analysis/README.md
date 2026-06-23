# Metabolomics Analysis Code

This directory contains the split v29 metabolomics workflow derived from:
`FUSCC_metabonomics_analysis/src/metabonomics_analysis_template_v29.R`.

Run order is controlled by `00.run_all_metabonomics_analysis.R`. The legacy
`1.metabonomics_analysis.R` file is now only a compatibility entry point.

Modules:

1. `01.setup_and_config.R`: libraries, `base_symbol`, input/output directories, constants.
2. `02.general_helpers.R`: utility functions, LC-MS import, limma wrappers.
3. `03.targeted_panel_helpers.R`: targeted metabolite panel definitions and heatmaps.
4. `04.cross_omics_helpers.R`: species, MetaCyc, EC, and taxonomy helper functions.
5. `05.model_and_figure_helpers.R`: panel-score models and main-figure helper functions.
6. `06.metabolomics_core_analysis.R`: metadata processing, metabolite QC, limma models, metabolite plots.
7. `07.cross_omics_support_analysis.R`: cross-omics support analyses and candidate species outputs.
8. `08.composite_figures.R`: final composite figures and session info.

Edit `base_symbol` in `01.setup_and_config.R` to match the project root before running.
Outputs follow the main microbiome analysis convention and are written under
`downstream_analysis/output/metabonomics_analysis/`.
