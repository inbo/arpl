is_leeg <- function(r) {
  if (is.null(r)) return(TRUE)
  if (!terra::hasMinMax(r)) suppressWarnings(terra::setMinMax(r))
  all(is.na(terra::minmax(r)))
}

get_param <- function(df_rij, col_naam, default_val = NA) {
  if (!(col_naam %in% colnames(df_rij)) || nrow(df_rij) == 0) return(default_val)
  val <- df_rij[[col_naam]][1]
  if (is.null(val) || is.na(val) || val == "" || val == "NA") return(default_val)
  return(val)
}

get_path <- function(..., paths = NULL, fallback = NULL) {
  velden <- c(...)
  for (v in velden) {
    val <- if (!is.null(paths)) paths[[v]] else NULL
    if (is.null(val) && exists("PATHS") && is.list(PATHS)) val <- PATHS[[v]]
    if (!is.null(val) && !is.na(val) && val != "" && file.exists(val)) return(val)
  }
  if (!is.null(fallback) && file.exists(fallback)) return(fallback)
  return(NULL)
}
