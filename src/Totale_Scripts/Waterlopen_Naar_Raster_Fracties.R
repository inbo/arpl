pad_vha_shp <- here::here("data/input/ASCI Files/Watervlakken_012_2026.shp")
pad_vha_tif <- here::here("data/input/Raster_Vlaanderen/vlaanderen_vha_cat012_fractie_10m.tif")
master_grid <- if (inherits(PATHS$master_grid, "SpatRaster")) PATHS$master_grid else terra::rast(PATHS$master_grid)

if (!file.exists(pad_vha_tif) && file.exists(pad_vha_shp)) {
  message("⏳ VHA cat. 0, 1, 2 voorbereiden (pure in-memory terra 1m->10m pipeline)...")
  t_start <- Sys.time()
  
  v_vha <- terra::vect(pad_vha_shp)
  if (!terra::same.crs(v_vha, master_grid)) {
    v_vha <- terra::project(v_vha, terra::crs(master_grid))
  }
  kol_cat <- names(v_vha)[grepl("cat", tolower(names(v_vha)))][1]
  if (!is.na(kol_cat)) {
    vals_cat <- suppressWarnings(as.integer(trimws(as.character(v_vha[[kol_cat]][[1]]))))
    v_vha <- v_vha[!is.na(vals_cat) & vals_cat %in% c(0L, 1L, 2L), ]
  }
  v_vha <- terra::disagg(v_vha[, 0])
  
  # 25 verticale stroken (~10 km breed -> past volledig in RAM, ~2-3 sec per strook)
  n_tiles <- 25
  e <- terra::ext(master_grid)
  x_breaks <- round(seq(e$xmin, e$xmax, length.out = n_tiles + 1) / 10) * 10
  x_breaks[1] <- e$xmin
  x_breaks[n_tiles + 1] <- e$xmax
  
  pb <- utils::txtProgressBar(min = 0, max = n_tiles, style = 3)
  tile_list <- vector("list", n_tiles)
  
  for (i in seq_len(n_tiles)) {
    e_strook <- terra::ext(x_breaks[i], x_breaks[i + 1], e$ymin, e$ymax)
    v_sub <- terra::crop(v_vha, e_strook)
    
    if (nrow(v_sub) > 0) {
      # Beperk ymin/ymax tot het echte bereik van de waterlopen in deze strook (afgerond op 10m)
      es <- terra::ext(v_sub)
      ymin_sub <- max(e$ymin, floor(es$ymin / 10) * 10)
      ymax_sub <- min(e$ymax, ceiling(es$ymax / 10) * 10)
      
      # 1m raster in RAM (~100-300 MB) -> snelle C++ scanline rasterize -> aggregate 10x10
      r_1m_grid  <- terra::rast(xmin = x_breaks[i], xmax = x_breaks[i + 1],
                                ymin = ymin_sub, ymax = ymax_sub,
                                resolution = 1, crs = terra::crs(master_grid))
      r_1m_bin   <- terra::rasterize(v_sub, r_1m_grid, field = 1L, background = 0L)
      r_10m_tile <- terra::aggregate(r_1m_bin, fact = 10, fun = "mean")
      r_10m_tile <- terra::ifel(r_10m_tile > 0, r_10m_tile, NA_real_)
      
      tile_list[[i]] <- r_10m_tile
      rm(r_1m_grid, r_1m_bin)
    }
    rm(v_sub)
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  rm(v_vha)
  gc()
  
  message("🧩 Stroken samenvoegen en wegschrijven naar Raster_Vlaanderen...")
  tile_list <- Filter(Negate(is.null), tile_list)
  r_10m <- terra::sprc(tile_list) |> terra::merge()
  r_10m <- terra::extend(r_10m, master_grid)
  if (!terra::compareGeom(r_10m, master_grid, stopOnError = FALSE)) {
    r_10m <- terra::resample(r_10m, master_grid, method = "near")
  }
  
  terra::writeRaster(
    r_10m,
    pad_vha_tif,
    overwrite = TRUE,
    datatype  = "FLT4S",
    NAflag    = -9999,
    gdal      = "COMPRESS=LZW"
  )
  
  rm(tile_list, r_10m)
  gc()
  message(sprintf("✅ VHA 10m-fractieraster permanent opgeslagen in %.1f sec!", 
                  as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
}

PATHS$vha_waterlopen <- pad_vha_tif