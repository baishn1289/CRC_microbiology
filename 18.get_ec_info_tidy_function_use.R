
get_ec_info_tidy <- function(ec_codes, pause = 0.5) {
  results <- list()
  for (ec in ec_codes) {
    query_id <- paste0("ec:", ec)
    Sys.sleep(pause)
    tryCatch({
      recs <- keggGet(query_id)
      rows <- lapply(recs, function(data) {
        name_list <- if(!is.null(data$NAME)) paste(data$NAME, collapse = ";; ") else NA_character_
        sys_name  <- if(!is.null(data$SYSNAME)) paste(data$SYSNAME, collapse = ";; ") else NA_character_
        reaction  <- if(!is.null(data$REACTION)) paste(data$REACTION, collapse = ";; ") else NA_character_
        class_name<- if(!is.null(data$CLASS)) paste(data$CLASS, collapse = ";; ") else NA_character_
        data.frame(EC_Number = ec,
                   NAME_raw = name_list,
                   SYSNAME_raw = sys_name,
                   REACTION_raw = reaction,
                   CLASS_raw = class_name,
                   stringsAsFactors = FALSE)
      })
      results[[ec]] <- do.call(rbind, rows)
    }, error = function(e) {
      results[[ec]] <- data.frame(EC_Number = ec, NAME_raw = NA_character_,
                                  SYSNAME_raw = NA_character_, REACTION_raw = NA_character_,
                                  CLASS_raw = NA_character_, stringsAsFactors = FALSE)
      message("Error fetching ", ec, ": ", e$message)
    })
  }
  combined <- do.call(rbind, results)

  tidy <- combined %>%
    group_by(EC_Number) %>%
    summarise(
      All_NAMES = paste(unique(na.omit(NAME_raw)), collapse = ";; "),
      All_SYSNAMES = paste(unique(na.omit(SYSNAME_raw)), collapse = ";; "),
      REACTION = paste(unique(na.omit(REACTION_raw)), collapse = ";; "),
      CLASS = paste(unique(na.omit(CLASS_raw)), collapse = ";; "),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      Name_vec = list(str_trim(unlist(str_split(All_NAMES, ";;")))),
      Sys_vec  = list(str_trim(unlist(str_split(All_SYSNAMES, ";;")))),
      Preferred_Name = if(length(Sys_vec) >= 1 && !is.na(Sys_vec[1]) && Sys_vec[1] != "") {
        Sys_vec[1]
      } else if(length(Name_vec) >= 1 && !is.na(Name_vec[1]) && Name_vec[1] != "") {
        Name_vec[1]
      } else {
        NA_character_
      }
    ) %>%
    ungroup() %>%
    select(EC_Number, Preferred_Name, All_NAMES, All_SYSNAMES, REACTION, CLASS)
  
  return(tidy)
}
