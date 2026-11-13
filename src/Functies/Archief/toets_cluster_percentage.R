# ==============================================================================
# toets_cluster_percentage.R
# Toetst per cluster of subcomponent/matrixlaag minstens pct_drempel bedraagt
# ==============================================================================

toets_cluster_percentage <- function(r_basis_cluster, 
                                     cl_ids, 
                                     r_component, 
                                     pct_drempel, 
                                     perimeter_m = 0, 
                                     master_grid) {
  
  if (all(is.na(suppressWarnings(terra::minmax(cl_ids))))) {
    return(master_grid * NA_real_)
  }
  
  # Eventuele perimeter focal filter (bv. 50m venster rondom het cluster)
  if (perimeter_m > 0) {
    w_mat <- terra::focalMat(cl_ids, perimeter_m, type = 'circle')
    w_mat[w_mat > 0] <- 1
    cl_zone <- terra::focal(cl_ids, w = w_mat, fun = "max", na.rm = TRUE)
    rm(w_mat)
  } else {
    cl_zone <- cl_ids
  }
  
  # Totale oppervlakte per cluster
  df_tot <- terra::zonal(r_basis_cluster, cl_zone, fun = "sum", na.rm = TRUE)
  colnames(df_tot) <- c("ID", "Totaal_Sum")
  
  # Component oppervlakte per cluster
  r_comp_clean <- terra::ifel(!is.na(r_component) & r_component > 0, r_component, 0)
  comp_in_zone <- terra::mask(r_comp_clean, cl_zone)
  df_comp <- terra::zonal(comp_in_zone, cl_zone, fun = "sum", na.rm = TRUE)
  colnames(df_comp) <- c("ID", "Comp_Sum")
  rm(comp_in_zone, r_comp_clean)
  if (perimeter_m > 0) rm(cl_zone)
  
  # Bereken percentage
  dt_calc <- merge(data.table::as.data.table(df_tot), data.table::as.data.table(df_comp), by = "ID", all.x = TRUE)
  dt_calc[is.na(Comp_Sum), Comp_Sum := 0]
  dt_calc[, Pct := (Comp_Sum / Totaal_Sum) * 100]
  
  valide_ids <- dt_calc[!is.na(ID) & Pct >= pct_drempel, ID]
  rm(df_tot, df_comp, dt_calc)
  
  if (length(valide_ids) > 0) {
    r_valide <- terra::mask(r_basis_cluster, cl_ids %in% as.integer(valide_ids))
  } else {
    r_valide <- master_grid * NA_real_
  }
  
  return(r_valide)
}