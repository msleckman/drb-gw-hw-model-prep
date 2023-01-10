source("2a_process_nhd_downscaling/src/subset_closest_nhd.R")
source("2a_process_nhd_downscaling/src/munge_split_temp_dat.R")
source("2a_process_nhd_downscaling/src/write_data.R")
source("2a_process_nhd_downscaling/src/estimate_mean_width.R")
source("2a_process_nhd_downscaling/src/combine_nhd_input_drivers.R")

p2a_targets_list <- list(

  # Subset NHDv2 reaches that overlap the NHM network to only include those 
  # that have a corresponding catchment (and meteorological data)
  tar_target(
    p2a_dendritic_nhd_reaches_along_NHM_w_cats,
    p1_dendritic_nhd_reaches_along_NHM %>%
      filter(areasqkm > 0)
  ),
  
  # Match temperature monitoring locations to "mainstem" NHDPlusv2 flowline
  # reaches, i.e. those that NHD reaches that intersect the NHM river network.
  # Note that this site-to-segment matching procedure emulates the process
  # used in delaware-model-prep for modeling water temperature in the Delaware
  # River Basin. Here, match each site to an NHDPlusv2 reach, preferring reaches 
  # for which the downstream vertex (endpoint) is close to the site point.
  tar_target(
    p2a_drb_temp_sites_w_segs,
    subset_closest_nhd(nhd_lines = p2a_dendritic_nhd_reaches_along_NHM_w_cats,
                       sites = p1_drb_temp_sites_sf)
  ),
  
  # Match temperature observational time series to "mainstem" NHDPlusv2 flowline
  # reaches by COMID. Rename columns to match the column names used in the aggregated
  # temps data file in the temperature forecasting data release. 
  tar_target(
    p2a_drb_temp_obs_w_segs,
    p1_drb_temp_obs %>%
      left_join(y = p2a_drb_temp_sites_w_segs[,c("site_id","comid")], by = "site_id")
  ),
  
  # Resolve duplicate temperature observations and summarize the temperature data
  # to return one value for each COMID-date. By setting prioritize_nwis_sites to
  # FALSE we include all observations from non-NWIS sites in the summarized data.
  tar_target(
    p2a_drb_temp_obs_by_comid,
    p2a_drb_temp_obs_w_segs %>%
      munge_split_temp_dat(., prioritize_nwis_sites = FALSE) %>%
      rename(COMID = comid)
  ),
  
  # Save temperature observations mapped to NHDPlusv2 COMIDs as a zarr data store
  tar_target(
    p2a_drb_temp_obs_w_segs_zarr,
    write_df_to_zarr(p2a_drb_temp_obs_by_comid, 
                     index_cols = c("date", "COMID"), 
                     "2a_process_nhd_downscaling/out/drb_temp_observations_nhdv2.zarr"),
    format = "file"
  ),
  
  # Estimate mean width for each "mainstem" NHDv2 reach. 
  # Note that one NHM segment, segidnat 1721 (subsegid 287_1) is not included
  # in the dendritic nhd reaches w cats data frame because the only COMID that
  # intersects the segment (COMID 4188275) does not have an NHD catchment area. 
  # So in addition to estimating width for the COMIDs represented in 
  # p2a_dendritic_nhd_reaches_along_NHM_w_cats, we also want to estimate width 
  # for COMID 4188275.
  tar_target(
    p2a_nhd_mainstem_reaches_w_width,
    {
      nhd_lines <- bind_rows(p2a_dendritic_nhd_reaches_along_NHM_w_cats,
                             filter(p1_dendritic_nhd_reaches_along_NHM,
                                    comid == "4188275"))
      estimate_mean_width(nhd_lines, 
                          estimation_method = 'nwis',
                          network_pos_variable = 'arbolate_sum',
                          ref_gages = p1_ref_gages_sf)
    }
  ),
  
  # Pull static segment attributes from PRMS SNTemp model driver data
  tar_target(
    p2a_static_inputs_prms,
    p1_sntemp_input_output %>%
      group_by(seg_id_nat) %>%
      summarize(seg_elev = unique(seg_elev),
                seg_slope = unique(seg_slope),
                seg_width = mean(seg_width, na.rm = TRUE)) %>%
      mutate(seg_id_nat = as.character(seg_id_nat)) %>%
      select(seg_id_nat, seg_elev, seg_slope, seg_width)
  ),
  
  # Pull dynamic segment attributes from PRMS SNTemp model driver data
  tar_target(
    p2a_dynamic_inputs_prms,
    p1_sntemp_input_output %>%
      mutate(seg_id_nat = as.character(seg_id_nat)) %>%
      select(seg_id_nat, date, seginc_potet)
  ),
  
  # Track the .py script as a file target so that downstream targets rebuild
  # if the .py script is edited. We do this because targets does not seem to 
  # automatically detect changes in functions within .py scripts.
  tar_target(
    p2a_py_file,
    "2a_process_nhd_downscaling/src/subset_nc_to_comid.py",
    format = "file"
  ),
  
  # Subset the DRB meteorological data to only include the NHDPlusv2 catchments 
  # (COMID) that intersect the NHM segments. `subset_nc_to_comid()` originally
  # developed by Jeff Sadler as part of the PGDL-DO project:
  # https://github.com/USGS-R/drb-do-ml/blob/main/2_process/src/subset_nc_to_comid.py
  # The resulting target is ~0.6 GB.
  tar_target(
    p2a_met_data_nhd_mainstem_reaches,
    {
      reticulate::source_python(p2a_py_file)
      subset_nc_to_comids(p1_drb_nhd_gridmet, 
                          p2a_dendritic_nhd_reaches_along_NHM_w_cats$comid) %>%
        as_tibble() %>%
        relocate(c(COMID,time), .before = "tmmx") %>%
        # format dates
        mutate(date = lubridate::as_date(time, tz = "UTC")) %>%
        # convert gridmet precip units from inches to meters, and temperature
        # units from degrees Farenheit to degrees Celsius. 
        # Note that we average the daily minimum and maximum temperatures from
        # gridmet and call that "seg_tave_air", which we assume corresponds 
        # approximately to the PRMS-SNTemp variable with the same name. 
        mutate(tmmean = rowMeans(select(., c(tmmn,tmmx))),
               seg_tave_air = ((tmmean - 32) * (5/9)),
               seg_rain = pr * 0.0254) %>%
        # rename gridmet columns to conform to PRMS-SNTemp names used in river-dl
        select(COMID, date, seg_tave_air, srad, seg_rain) %>%
        rename(seginc_swrad = srad)
    }
  ),
  
  # Compile river-dl static input drivers at NHDv2 resolution, including river 
  # width (m) slope (unitless), and min/max elevation (m). Note that these input 
  # drivers represent "mainstem" NHDv2 reaches only (i.e., those reaches that 
  # intersect the NHM fabric). Some of the columns in p2a_static_input_drivers_nhd 
  # are meant to facilitate comparison/EDA between segment attributes at the NHM and 
  # NHDPlusv2 scales (i.e., seg_slope ~ slope_len_wtd_mean; seg_elev ~ seg_elev_min; 
  # seg_width ~ seg_width_max). 
  tar_target(
    p2a_static_inputs_nhd,
    prepare_nhd_static_inputs(nhd_flowlines = p2a_nhd_mainstem_reaches_w_width,
                              prms_inputs = p2a_static_inputs_prms,
                              nhd_nhm_xwalk = p1_drb_comids_all_tribs)
  ),
  
  # Format NHD-scale static input drivers
  tar_target(
    p2a_static_inputs_nhd_formatted,
    p2a_static_inputs_nhd %>%
      select(COMID, seg_id_nat, subsegid, 
             est_width_m, min_elev_m, slope) %>%
      rename(seg_width_empirical = est_width_m,
             seg_elev = min_elev_m,
             seg_slope = slope,
             seg_id_nat = seg_id_nat)
  ),
  
  # Combine NHD-scale static input drivers with dynamic climate drivers. 
  # Note that there is currently no pot ET analog variable at the NHD-scale,
  # so we are using seginc_potet from PRMS-SNTemp and assuming that there is 
  # not much intra-segment variation in potential ET among NHD reaches that 
  # contribute to a given NHM segment. The resulting target contains 15,800 
  # days of climate data across each of 3,173 COMIDs = 50,133,400 total rows.
  tar_target(
    p2a_input_drivers_nhd,
    combine_nhd_input_drivers(nhd_static_inputs = p2a_static_inputs_nhd_formatted, 
                              climate_inputs = p2a_met_data_nhd_mainstem_reaches,
                              prms_dynamic_inputs = p2a_dynamic_inputs_prms,
                              earliest_date = "1979-01-01",
                              latest_date = "2022-04-04")
  ),
  
  # Save river-dl input drivers at NHDv2 resolution as a zarr data store
  tar_target(
    p2a_input_drivers_nhd_zarr,
    write_df_to_zarr(p2a_input_drivers_nhd, 
                     index_cols = c("date", "COMID"), 
                     "2a_process_nhd_downscaling/out/nhdv2_inputs_io.zarr"),
    format = "file"
  )
)
