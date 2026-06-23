
get_current_script_dir <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_hits <- command_args[startsWith(command_args, file_arg)]
  if (length(file_hits) > 0) {
    script_path <- sub(file_arg, "", file_hits[[1]], fixed = TRUE)
    return(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
  }

  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) {
      return(frame$ofile)
    }
    NA_character_
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    return(dirname(normalizePath(frame_files[[length(frame_files)]], winslash = "/", mustWork = FALSE)))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

module_dir <- get_current_script_dir()
module_files <- c(
  "01.setup_and_config.R",
  "02.general_helpers.R",
  "03.targeted_panel_helpers.R",
  "04.cross_omics_helpers.R",
  "05.model_and_figure_helpers.R",
  "06.metabolomics_core_analysis.R",
  "07.cross_omics_support_analysis.R",
  "08.composite_figures.R"
)

for (module_file in module_files) {
  message("Sourcing ", module_file)
  source(file.path(module_dir, module_file), chdir = TRUE, local = FALSE, encoding = "UTF-8")
}
