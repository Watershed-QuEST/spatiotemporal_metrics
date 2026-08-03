#### READ ME ####

# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Alex Webster, 2026-07-28 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-08-03

# This script builds the spatial infrastructure needed for synthesizing synthetic stream network data and conducting Monte-Carlo sensitivity analyses for each metric.
# This script is written for the "nm" watershed in the toy data set, and should be repeated or looped for additional watersheds. The user must know the sites and stream network well enough to verify that every site snapped to the correct stream in Part A and be able to adjust lat/lon if needed. 
# This script is organized as:
#   PART A -- delineate stream networks, subwatersheds, and snapped site locations from a DEM, for every site in "All sites.xlsx" toy dataset
#   PART B -- assemble SSNbler/SSN2 .ssn objects with the tail-up covariance weighted by real flow accumulation (sampled directly off the DEM-derived raster in Part A).
#   PART C -- save ssn objects for use in downstream scripts

# Requires: "All sites.xlsx" (sheet "sites": Site, Code, Lat, Lon), "NM_BR Toy dataset.csv", and internet access for get_elev_raster(). 
# Outputs:
#   1. streams.gpkg, sites_snapped.gpkg, subwatersheds.gpkg
#   2. <network>.ssn/  (the assembled SSN2 object)
#   3. ssn_objects.rds

#### Packages ####
library(tidyverse)
library(raster)
library(sf)
library(sp)
library(elevatr)
library(mapview)
library(stars)
library(readxl)
library(whitebox)
library(SSNbler)
library(SSN2)
library(googledrive)

# raster package messes with dplyr::select functionality. This fixes the package conflicts when functions are called: 
select  <- dplyr::select
extract <- raster::extract

#### Read in data ####

# List all files in the folder
toy_files <- drive_ls(drive_get("https://drive.google.com/drive/u/1/folders/1zh0YTDM5w971iFwmw-iSyTDQQ4MyGL8-"))
# Download files
# googledrive::drive_download(file = toy_files$id[toy_files$name=="NM-BR Toy dataset.csv"], 
#                             path = "drivedata/toy.csv",
#                             overwrite = T)
googledrive::drive_download(file = toy_files$id[toy_files$name=="All sites.xlsx"], 
                            path = "drivedata/All sites.xlsx",
                            overwrite = T)
# read in csv
# toy = read.csv("drivedata/toy.csv")

#### Configure in/outputs for flow accumulation workflow ####

## This named list should make it easy to adapt this whole pipeline to a different single watershed later by adding a new entry here.
networks <- list(
  nm = list(code = "NM", label = "santafe",    epsg = 32613,   # UTM 13N
            z = 11, expand_m = 17000,
            lsn_path = "nm_lsn", ssn_path = "SantaFe.ssn")
)

## Which network(s) to actually run below. Note that following scripts are currently hard-coded to nm specifically.
active_networks <- c("nm")

flow_threshold  <- 30    # min contributing cells to be called a "stream"
snap_dist       <- 50    # meters, for snapping pour points onto the stream
min_trib_length <- 100   # meters -- dangling leaf tributary segments shorter than this are pruned (see prune_short_tributaries() below). Start small and increase just enough to clear whatever omplex-confluence nodes the Part B diagnostic flags -- check the printed prune-pass output and the before/after mapview to make sure you're not also trimming real headwater reaches.
crop_buffer_m   <- 300   # meters -- margin added around crop_extent before cropping the DEM for pass-2 hydrology.

## Manual coordinate corrections. If a site's Lat/Lon in "All sites.xlsx" puts it off its correct stream reach (or close enough that it snaps onto the wrong nearby tributary), fix it here rather than editing the spreadsheet by hand: delineate_network() prints a map for each network with a live lat/lon readout in the top-left corner as you hover -- move the mouse to the correct spot on the correct stream, read off the coordinates, add a row below, and rerun. These corrections are applied in-memory only; "All sites.xlsx" itself is never modified. 
# Example:
##   site_corrections <- tibble::tribble(
##     ~Site,   ~Lat,     ~Lon,
##     "USF07", 35.6789,  -105.9876
##   )

site_corrections <- tibble::tribble(
     ~Site,   ~Lat,     ~Lon,
     "USF01", 35.68861,  -105.82576,
     "USF02", 35.68807,  -105.82349,
     "USF12", 35.68882,  -105.82326,
     "USF25", 35.69586,  -105.81313,
     "USF26", 35.70320,  -105.80708,
     "USF03", 35.70854,  -105.80553,
     "USF08", 35.71694,  -105.80246,
     "USF04", 35.72173,  -105.79843,
     "USF05", 35.72121,  -105.79748,
     "USF22", 35.73158,  -105.79444,
     "USF11", 35.74519,  -105.77631,
     "USF09", 35.74927,  -105.77855,
     "USF21", 35.77917,  -105.77402,
     "USF14", 35.77974,  -105.77333,
     "USF13", 35.77980,  -105.77375,
     "USF28", 35.78104,  -105.77210,
     "USF29", 35.78144,  -105.77128,
     "USF18", 35.78328,  -105.77111
   )

## Sites to exclude entirely from the analysis -- e.g. old/historical sites
excluded_sites <- c("USF01", "USF02", "USF30")
all_sites_raw <- read_excel("drivedata/All sites.xlsx", sheet = "sites") %>%
  filter(!is.na(Site))
all_sites <- all_sites_raw %>% filter(!Site %in% excluded_sites)

set.seed(42)

#### PART A -- DEM Delineation ####

# make/check file structure for outputs
data_out_dir <- "data"  # working data (ssn_objects.rds) 
for (d in c("temp", "geo_output", data_out_dir)) if (!dir.exists(d)) dir.create(d)

## Define functions that fill -> breach -> D8 pointer -> D8 flow accum -> extract streams -> raster-to-vector, given a locally-stored DEM
# Called twice per network: 1. generous-extent DEM to find the watershed boundary, then 2. the cropped DEM for the final clean network
run_hydrology_steps <- function(dem_path, wd, prefix, threshold) {
  fill_path      <- file.path("temp", paste0(prefix, "_fill.tif"))
  breach_path    <- file.path("temp", paste0(prefix, "_breach.tif"))
  fdir_path      <- file.path("temp", paste0(prefix, "_flowdir.tif"))
  faccum_path    <- file.path("temp", paste0(prefix, "_flowaccum.tif"))
  streams_r_path <- file.path("temp", paste0(prefix, "_streams.tif"))
  streams_v_path <- file.path("temp", paste0(prefix, "_streams.shp"))

  wbt_fill_single_cell_pits(dem = dem_path, output = fill_path, wd = wd)
  wbt_breach_depressions(dem = fill_path, output = breach_path, wd = wd)
  wbt_d8_pointer(dem = breach_path, output = fdir_path, wd = wd)
  wbt_d8_flow_accumulation(input = breach_path, output = faccum_path, wd = wd)
  wbt_extract_streams(flow_accum = faccum_path, output = streams_r_path,
                       threshold = threshold, wd = wd)
  wbt_raster_streams_to_vector(streams = streams_r_path, d8_pntr = fdir_path,
                                output = streams_v_path, wd = wd)

  list(fdir = fdir_path, faccum = faccum_path, streams_vec = streams_v_path)
}

## This helper function is defined to trim short, dangling leaf tributary segments, used in the delineate_network function workflow below.
# A "leaf" edge here is one whose upstream endpoint isn't shared with any other edge -- i.e. nothing flows into it, so removing a short leaf either drops a genuine-but-tiny first-order stub, or, at a complex-confluence node, brings a 3+-way junction back down to the 2-in/1-out topology ssn_assemble() requires. max_iter argument allows it to iterate because removing one short stub can occasionally expose another short stub just upstream of it (a chain of two tiny rasterization-artifact segments, say).
prune_short_tributaries <- function(streams_sf, min_length_m, coord_tol = 0.01, max_iter = 10) {
  for (i in seq_len(max_iter)) {
    ends <- streams_sf %>%
      st_coordinates() %>%
      as.data.frame() %>%
      group_by(L1) %>%
      summarize(x1 = first(X), y1 = first(Y), x2 = last(X), y2 = last(Y), .groups = "drop")

    node_key <- function(x, y) paste(round(x / coord_tol) * coord_tol, round(y / coord_tol) * coord_tol)
    start_key <- node_key(ends$x1, ends$y1)
    end_key   <- node_key(ends$x2, ends$y2)
    node_degree <- table(c(start_key, end_key))

    is_leaf  <- as.integer(node_degree[start_key]) == 1 | as.integer(node_degree[end_key]) == 1
    length_m <- as.numeric(st_length(streams_sf))
    to_drop  <- which(is_leaf & length_m < min_length_m)

    if (length(to_drop) == 0) break
    cat("    prune pass", i, "- removing", length(to_drop),
        "short leaf segment(s), length(s) (m):",
        paste(round(length_m[to_drop], 1), collapse = ", "), "\n")
    streams_sf <- streams_sf[-to_drop, ]
  }
  streams_sf
}

## This helper function is defined to split a complex 3-way confluence into a 2-in/1-out topology as required by ssn_assemble(). 
# This is the R version of nudging one incoming line's endpoint a few meters downstream in QGIS as SSN2 package prompts the user to do when it encounters 3-way confluences. It detaches one edge touching the node and reattaches it a short distance along another edge at that same node, inserting a new vertex and splitting that edge in two. By default it detaches the edge with the LOWEST flowaccum touching the node (the smallest/least important tributary) and reattaches it along the edge with the HIGHEST flowaccum at that node. The user should override detach_edge/target_edge (row indices into streams_sf) if the function picks the wrong edges for a particular node -- inspect the printed flowaccum values and the mapview plot to check.
# Note that this function only detaches ONE edge per call, so a node with more than 4 edges touching it (more than 3 tributaries) needs a second call on the same node_xy afterward.
split_confluence <- function(streams_sf, node_xy, nudge_dist = 5, coord_tol = 0.01,
                              detach_edge = NULL, target_edge = NULL) {
  node_x <- node_xy[1]; node_y <- node_xy[2]
  node_key <- function(x, y) paste(round(x / coord_tol) * coord_tol, round(y / coord_tol) * coord_tol)
  target_key <- node_key(node_x, node_y)

  ends <- streams_sf %>%
    st_coordinates() %>%
    as.data.frame() %>%
    group_by(L1) %>%
    summarize(x1 = first(X), y1 = first(Y), x2 = last(X), y2 = last(Y), .groups = "drop") %>%
    mutate(touches_start = node_key(x1, y1) == target_key,
           touches_end   = node_key(x2, y2) == target_key,
           touches       = touches_start | touches_end)

  touching_idx <- ends$L1[ends$touches]
  if (length(touching_idx) < 4) {
    stop("Node (", node_x, ", ", node_y, ") has only ", length(touching_idx),
         " edge(s) touching it -- nothing to split (need >= 4; call again if",
         " a previous split left this node still with > 3).")
  }

  fa <- streams_sf$flowaccum[touching_idx]
  if (is.null(detach_edge)) detach_edge <- touching_idx[which.min(fa)]
  if (is.null(target_edge)) {
    candidates  <- touching_idx[touching_idx != detach_edge]
    target_edge <- candidates[which.max(streams_sf$flowaccum[candidates])]
  }

  target_geom   <- st_geometry(streams_sf)[[target_edge]]
  target_coords <- st_coordinates(target_geom)[, 1:2, drop = FALSE]
  node_is_first <- ends$touches_start[ends$L1 == target_edge]
  ordered <- if (node_is_first) target_coords else target_coords[nrow(target_coords):1, , drop = FALSE]

  seg_len <- sqrt(diff(ordered[, 1])^2 + diff(ordered[, 2])^2)
  cum_len <- c(0, cumsum(seg_len))
  if (max(cum_len) <= nudge_dist) {
    stop("Target edge (row ", target_edge, ") is only ", round(max(cum_len), 1),
         "m long -- shorter than nudge_dist (", nudge_dist, "m). Pick a smaller",
         " nudge_dist or specify a different target_edge.")
  }
  ins_i <- which(cum_len >= nudge_dist)[1]
  frac  <- (nudge_dist - cum_len[ins_i - 1]) / (cum_len[ins_i] - cum_len[ins_i - 1])
  new_pt <- ordered[ins_i - 1, ] + frac * (ordered[ins_i, ] - ordered[ins_i - 1, ])

  near_part <- rbind(ordered[seq_len(ins_i - 1), , drop = FALSE], new_pt)
  far_part  <- rbind(new_pt, ordered[ins_i:nrow(ordered), , drop = FALSE])

  out <- streams_sf
  st_geometry(out)[[target_edge]] <- st_linestring(near_part)

  far_row <- out[target_edge, ]
  st_geometry(far_row) <- st_sfc(st_linestring(far_part), crs = st_crs(out))
  out <- rbind(out, far_row)

  detach_coords <- st_coordinates(st_geometry(out)[[detach_edge]])[, 1:2, drop = FALSE]
  if (ends$touches_start[ends$L1 == detach_edge]) {
    detach_coords[1, ] <- new_pt
  } else {
    detach_coords[nrow(detach_coords), ] <- new_pt
  }
  st_geometry(out)[[detach_edge]] <- st_linestring(detach_coords)

  cat("  split confluence at (", round(node_x, 1), ",", round(node_y, 1),
      "): detached edge", detach_edge, "(flowaccum =",
      round(streams_sf$flowaccum[detach_edge], 1), ") onto a new node", nudge_dist,
      "m along edge", target_edge, "(flowaccum =",
      round(streams_sf$flowaccum[target_edge], 1), ")\n")
  out
}

## This defines the central function that incorporates the previous two helper functions to delineate the network:
delineate_network <- function(net_name, cfg, sites_all) {

  sites <- sites_all %>% filter(Code == cfg$code)
  if (nrow(sites) == 0) stop("No sites found for Code == ", cfg$code)

  pour <- st_as_sf(sites, coords = c("Lon", "Lat"), crs = 4326) %>%
    st_transform(crs = cfg$epsg)

  wd <- normalizePath(getwd())   # WhiteboxTools needs an absolute path; computed here so it isn't tied to one user's machine

  ## PULL A DEM of generous extent, covering all pour points at once
  dem <- get_elev_raster(pour, z = cfg$z, clip = "bbox", expand = cfg$expand_m)
  dem_path <- file.path("temp", paste0(cfg$label, "_dem.tif"))
  writeRaster(dem, dem_path, overwrite = TRUE)

  ## FIRST PASS: hydrology on the full-extent DEM, just to get watershed boundaries so the DEM can be cropped for a cleaner final network.
  step1 <- run_hydrology_steps(dem_path, wd, paste0(cfg$label, "_pass1"), flow_threshold)
  streams1 <- st_read(step1$streams_vec, quiet = TRUE); st_crs(streams1) <- cfg$epsg

  ## Apply any manual coordinate corrections for this network's sites (see site_corrections in the config block above), then rebuild "pour" from the corrected table so everything downstream uses the fixed locations.
  net_corrections <- site_corrections %>% filter(Site %in% sites$Site)
  if (nrow(net_corrections) > 0) {
    sites <- rows_update(sites, net_corrections, by = "Site", unmatched = "ignore")
    pour <- st_as_sf(sites, coords = c("Lon", "Lat"), crs = 4326) %>%
      st_transform(crs = cfg$epsg)
    cat(net_name, ": applied", nrow(net_corrections), "manual coordinate correction(s) for:",
        paste(net_corrections$Site, collapse = ", "), "\n")
  }
  ## NOTE: check that every site snapped to the correct stream, not a neighboring one or not on a stream at all. If a site lands on the wrong stream or no stream, options are to a) add a fix to site_corrections above, or b) if the stream doesn't exist, you may need to adjust flow_threshold if small tribs are missing/spurious, and rerun.
  # To easily adjust lat/lon, hover over this map (top-left corner shows live lat/lon as you move the mouse). Pick a lat/lon that lands EXACTLY on the correct stream, well away from a confluence, use the hover to read and add correct coordinates to site_corrections in the config block above, then rerun.
  print((mapview(streams1) + mapview(pour)) %>% leafem::addMouseCoordinates())

  outlet_path <- file.path("temp", paste0(cfg$label, "_pour.shp"))
  st_write(pour, outlet_path, delete_layer = TRUE)

  snap_path <- file.path("temp", paste0(cfg$label, "_pour_snap.shp"))
  wbt_snap_pour_points(pour_pts = outlet_path, flow_accum = step1$faccum,
                        snap_dist = snap_dist, output = snap_path, wd = wd)

  pour_snap <- st_read(snap_path, quiet = TRUE); st_crs(pour_snap) <- cfg$epsg
  if (nrow(pour_snap) != nrow(sites)) {
    stop(net_name, ": wbt_snap_pour_points() returned ", nrow(pour_snap), " points but ",
         nrow(sites), " sites were sent in -- it silently dropped point(s) that didn't ",
         "snap within snap_dist (", snap_dist, "m). The Site<->row assignment below is not ",
         "safe until this is resolved (increase snap_dist, or check for a bad Lat/Lon).")
  }
  pour_snap$Site <- sites$Site 
  print((mapview(streams1) + mapview(pour) + mapview(pour_snap, color = "red")) %>%
          leafem::addMouseCoordinates())

  ## This chunk delinates each site's FULL cumulative watershed, delineated with its own dedicated pour point via wbt_watershed() call -- NOT one batched multi-point call since it can't handle nested watersheds correctly. One call per site is slower but is the only way to get correct cumulative areas.
  delineate_one_shed <- function(site_name) {
    pt_path <- file.path("temp", paste0(cfg$label, "_", site_name, "_pour.shp"))
    st_write(pour_snap[pour_snap$Site == site_name, ], pt_path, delete_layer = TRUE)
    shed_path_i <- file.path("temp", paste0(cfg$label, "_", site_name, "_shed.tif"))
    wbt_watershed(d8_pntr = step1$fdir, pour_pts = pt_path, output = shed_path_i, wd = wd)
    poly <- st_as_stars(raster(shed_path_i)) %>% st_as_sf(merge = TRUE)
    if (nrow(poly) == 0) {
      cat(net_name, "-", site_name, ": watershed came back empty (point may not",
          "have snapped onto a real flow path) -- skipping.\n")
      return(NULL)
    }
    # A single-point watershed can create small raster artifacts. This creates a union to one feature so every site gets exactly one row below, with its total area intact.
    st_sf(Site = site_name, geometry = st_union(st_geometry(poly)))
  }
  sheds <- map(sites$Site, delineate_one_shed) %>% compact() %>% bind_rows()

  missing_sheds <- setdiff(sites$Site, sheds$Site)
  if (length(missing_sheds) > 0) {
    cat(net_name, "-", length(missing_sheds), "site(s) got no watershed at all:\n")
    print(missing_sheds)
  }
  
  area_check <- sheds %>%
    mutate(Area_m2_delineated = as.numeric(st_area(.))) %>%
    st_drop_geometry() %>%
    left_join(sites %>% dplyr::select(Site), by = "Site")

    ## This uses the site with the largest DELINEATED area (sheds, computed just above) as outlet.
    outlet_site <- area_check$Site[which.max(area_check$Area_m2_delineated)]

  crop_extent <- sheds[sheds$Site == outlet_site, ]
  cat("Crop extent (outlet =", outlet_site, "): n polygons =", nrow(crop_extent),
      "| area (km2) =", round(as.numeric(sum(st_area(crop_extent))) / 1e6, 2), "\n")
  if (nrow(crop_extent) != 1) {
    stop(net_name, ": expected exactly 1 outlet watershed polygon, got ", nrow(crop_extent),
         " -- check that '", outlet_site, "' snapped onto the right stream (see the mapview plot above).")
  }

  ## SECOND PASS: crop DEM to the outlet's watershed and rerun hydrology for a clean stream network restricted to the actual watershed extent. Note the FINAL exported network is still clipped to the true, unbuffered crop_extent on the next line, same as before.
  cropped_dem <- raster::crop(dem, st_buffer(crop_extent, crop_buffer_m))
  cropped_dem_path <- file.path("temp", paste0(cfg$label, "_dem_cropped.tif"))
  writeRaster(cropped_dem, cropped_dem_path, overwrite = TRUE)

  step2 <- run_hydrology_steps(cropped_dem_path, wd, paste0(cfg$label, "_pass2"), flow_threshold)

  streams_final <- st_read(step2$streams_vec, quiet = TRUE)
  st_crs(streams_final) <- cfg$epsg
  streams_final <- streams_final[crop_extent, ]

  ## wbt_raster_streams_to_vector() writes its own "FID"-type attribute column into the streams shapefile as a byproduct of the raster -> vector conversion (not anything we need here). SSNbler's lines_to_lsn() (Part B, below) creates its own "fid" reach-ID column and errors, so we'll drop it here so Part B doesn't choke downstream.
  streams_final <- streams_final %>% select(-matches("^fid$", ignore.case = TRUE))

  n_before <- nrow(streams_final)
  streams_final <- prune_short_tributaries(streams_final, min_trib_length)
  cat(net_name, ": pruned", n_before - nrow(streams_final), "of", n_before,
      "stream segment(s) shorter than", min_trib_length, "m (dangling leaf tips only).\n")

  ## Safety check: make sure pruning didn't strand a site by removing the one segment nearest to it. Compares each site's already-known snapped location (pour_snap, from the pass-1 network, before pruning) against the pruned network -- if a site is now farther than snap_dist from every remaining edge, this warns loudly rather than let it silently fail to snap onto the pruned network a few steps down (or, worse, silently snap onto the wrong nearby reach).
  site_dists <- st_distance(pour_snap, streams_final) %>% apply(1, min) %>% as.numeric()
  stranded <- pour_snap$Site[site_dists > snap_dist]
  if (length(stranded) > 0) {
    warning(net_name, ": pruning with min_trib_length = ", min_trib_length,
            "m removed the segment nearest to site(s): ", paste(stranded, collapse = ", "),
            " -- lower min_trib_length (or check these sites' placement) before continuing.")
  }

  ## This samples flow accumulation onto the final stream vector, needed for Part B. This samples at both endpoints of each line and takes the max, since flow accumulation is always higher at a segment's downstream end regardless of how the line's vertices happen to be ordered.
  faccum_r <- raster(step2$faccum)
  endpoints <- streams_final %>%
    st_coordinates() %>%
    as.data.frame() %>%
    group_by(L1) %>%
    summarize(x_first = first(X), y_first = first(Y),
              x_last = last(X), y_last = last(Y), .groups = "drop")
  acc_first <- raster::extract(faccum_r, endpoints[, c("x_first", "y_first")])
  acc_last  <- raster::extract(faccum_r, endpoints[, c("x_last", "y_last")])
  streams_final$flowaccum <- pmax(acc_first, acc_last, na.rm = TRUE)

  ## Re-snap sites onto the cropped/final network so the exported site layer matches the exported stream layer exactly.
  final_outlet_path <- file.path("temp", paste0(cfg$label, "_pour_final.shp"))
  st_write(pour, final_outlet_path, delete_layer = TRUE)
  final_snap_path <- file.path("temp", paste0(cfg$label, "_pour_snap_final.shp"))
  wbt_snap_pour_points(pour_pts = final_outlet_path, flow_accum = step2$faccum,
                        snap_dist = snap_dist, output = final_snap_path, wd = wd)
  sites_snapped <- st_read(final_snap_path, quiet = TRUE)
  st_crs(sites_snapped) <- cfg$epsg
  sites_snapped$Site <- sites$Site  # same order as pour/final_outlet_path
  # strips FID-type column again
  sites_snapped <- sites_snapped %>% select(-matches("^fid$", ignore.case = TRUE))

  out_dir <- file.path("geo_output", net_name)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  st_write(streams_final, file.path(out_dir, "streams.gpkg"), delete_layer = TRUE)
  st_write(sites_snapped, file.path(out_dir, "sites_snapped.gpkg"), delete_layer = TRUE)
  st_write(sheds, file.path(out_dir, "subwatersheds.gpkg"), delete_layer = TRUE)
  write_csv(area_check, file.path(out_dir, "area_check.csv"))
  cat("\nWrote:", out_dir, "/streams.gpkg, sites_snapped.gpkg, subwatersheds.gpkg, area_check.csv\n")

  ## faccum_path and cell_area_m2 are returned so Part B'4.5's 's make_pred_points() can sample drainage area (flowaccum cell count x cell area) at the synthetic prediction points the same way flowaccum was just sampled onto streams_final above. Those points need a logArea value too, since every ssn_lm() fixed-effect formula in next script includes logArea, and predict() needs it at every prediction location, not just the real sites.
  list(streams = streams_final, sites = sites_snapped, sheds = sheds, area_check = area_check,
       faccum_path = step2$faccum, cell_area_m2 = prod(raster::res(faccum_r)))
}

## Time to use all the functions defined above! imap applies the delineate_networks function to each element- used here and elsewhere in the script to work through all the elements in the stream network efficiently. This will take a minute to run through each individual subcatchment. Rerun to here if points need to be corrected.
geo_results <- imap(networks[active_networks], ~ delineate_network(.y, .x, all_sites))

#### PART B -- Assemble .ssn objects ####

## Requires geo_results from Part A.
## Note that there is an iterative fix section in here (called out as its own section) so do not run this through without inspection

## B1. Build the landscape network (LSN) and check topology ----

edges <- imap(geo_results, function(g, name) {
  lines_to_lsn(
    streams = g$streams,
    lsn_path = networks[[name]]$lsn_path,
    snap_tolerance = 1,
    topo_tolerance = 10,
    check_topology = TRUE, ## check_topology = TRUE flags dangling nodes, unsnapped intersections, and complex confluences, and produces node_errors.gpkg which is used in corrective work flows below. 
    overwrite = TRUE,
    verbose = TRUE
  )
})

## Surfaces anything lines_to_lsn() flagged with details
for (name in names(edges)) {
  err_path <- file.path(networks[[name]]$lsn_path, "node_errors.gpkg")
  if (file.exists(err_path)) {
    node_errors <- st_read(err_path, quiet = TRUE)
    cat("\n", name, "- lines_to_lsn() flagged", nrow(node_errors),
        "topology issue(s) in node_errors.gpkg:\n")
    print(node_errors)
  }
}

## B1. START of iterative fix for complex confluences ----
## SSN2/SSNbler need a strictly dendritic (binary-tree) network -- every node can have at most 2 upstream edges flowing into 1 downstream edge. The functions below find and fix these. These need to be run once per complex confluence found to make the correction work!

find_complex_confluences <- function(edges_sf, coord_tol = 0.01) {
  ends <- edges_sf %>%
    st_coordinates() %>%
    as.data.frame() %>%
    group_by(L1) %>%
    summarize(x1 = first(X), y1 = first(Y), x2 = last(X), y2 = last(Y), .groups = "drop")
  pts <- bind_rows(
    ends %>% transmute(x = round(x1 / coord_tol) * coord_tol, y = round(y1 / coord_tol) * coord_tol),
    ends %>% transmute(x = round(x2 / coord_tol) * coord_tol, y = round(y2 / coord_tol) * coord_tol)
  )
  ## A normal confluence has exactly 3 edge-ends meeting (2 upstream + 1
  ## downstream) -- more than that is the "complex confluence" ssn_assemble()
  ## rejects.
  pts %>% count(x, y, name = "n_edge_ends") %>% filter(n_edge_ends > 3)
}
# DIAGNOSTIC: surface bad nodes
for (name in names(edges)) {
  bad_nodes <- find_complex_confluences(edges[[name]])
  if (nrow(bad_nodes) > 0) {
    cat("\n", name, "-", nrow(bad_nodes), "node(s) with more than 2 edges",
        "flowing in -- almost certainly what ssn_assemble() will choke on:\n")
    print(bad_nodes)
    bad_pts <- st_as_sf(bad_nodes, coords = c("x", "y"), crs = st_crs(edges[[name]]))
    print(mapview(edges[[name]]) + mapview(bad_pts, color = "red", cex = 8))
    cat("Options: (1) lower min_trib_length isn't enough if this is a real",
        "3-way confluence, not a short rasterization artifact; (2) split it",
        "with split_confluence() -- see the commented example below -- then",
        "rerun lines_to_lsn() (B1) onward; (3) as a last resort, nudge one",
        "incoming line's endpoint by hand in QGIS.\n")
  }
}
# If there are bad notes, save them to fix
bad_nodes <- find_complex_confluences(edges$nm)
node_xy <- as.numeric(bad_nodes[1, c("x", "y")])
# rerun split_confluence to fix problem nodes
geo_results$nm$streams <- split_confluence(geo_results$nm$streams, node_xy = node_xy)
# resaves edges with fix
edges <- imap(geo_results, function(g, name) {
  lines_to_lsn(streams = g$streams, lsn_path = networks[[name]]$lsn_path,
               snap_tolerance = 1, topo_tolerance = 10, check_topology = TRUE,
               overwrite = TRUE, verbose = TRUE)
})
# run "DIAGNOSTIC: surface bad nodes" for loop above NOW to see if any remain!
## B1. END OF iterative fix for complex confluences ----

## B2. Add tail-up weighting (afv in SSN2 terms) from flow accumulation ----

edges <- imap(edges, function(e, name) {
  afv_edges(
    edges = e,
    lsn_path = networks[[name]]$lsn_path,
    infl_col = "flowaccum",
    segpi_col = "flowPI",
    afv_col = "afvFlow",
    overwrite = TRUE,
    save_local = TRUE
  )
})

## B3. Calculate upstream distances for edges ----

edges <- imap(edges, function(e, name) {
  updist_edges(
    edges = e,
    lsn_path = networks[[name]]$lsn_path,
    calc_length = TRUE,
    length_col = "Length_m",
    overwrite = TRUE,
    save_local = TRUE,
    verbose = TRUE
  )
})

## B4. Map LSN, AFV, and upstream distances to sites ----

obs <- imap(geo_results, function(g, name) {
  sites_to_lsn(
    sites = g$sites,
    edges = edges[[name]],
    lsn_path = networks[[name]]$lsn_path,
    snap_tolerance = 100,  # generous starting point for GPS coords- check snapdist in the output and tighten if needed
    save_local = TRUE,
    file_name = "sites.gpkg",
    overwrite = TRUE
  )
})

obs <- imap(obs, function(o, name) {
  afv_sites(
    sites = list(obs = o),
    edges = edges[[name]],
    afv_col = "afvFlow",
    lsn_path = networks[[name]]$lsn_path,
    save_local = TRUE,
    overwrite = TRUE
  )$obs
})

obs <- imap(obs, function(o, name) {
  updist_sites(
    sites = list(obs = o),
    edges = edges[[name]],
    length_col = "Length_m",
    lsn_path = networks[[name]]$lsn_path,
    save_local = TRUE,
    overwrite = TRUE
  )$obs
})

## B4.5 Generate synthetic prediction points along the stream network ----
## Output from this is included in the .ssn object as a prediction dataset ("predpts") so it can be used in future scripts to drive network-realistic Monte Carlo analyses directly off SSN2's fitted covariance.

pred_spacing_m <- 250  # distance between synthetic prediction points along the network -- fewer/coarser points means a smaller covariance matrix to invert/decompose in future scripts; start coarse and only go finer if the MC analysis needs more points than this produces

# this function defines the workflow to generate synthetic prediction points
make_pred_points <- function(edges_sf, spacing_m, faccum_path, cell_area_m2) {
  pts <- st_line_sample(edges_sf, density = 1 / spacing_m, type = "regular")
  pts <- st_sf(geometry = pts)

  ## this defines any edge shorter than spacing_m to get 0 points from density-based sampling (rounds down), which st_line_sample() represents as an EMPTY geometry for that row.
  empty <- st_is_empty(pts)
  if (any(empty)) {
    short_edges <- edges_sf[empty, ]
    half_len <- as.numeric(st_length(short_edges)) / 2
    ## st_line_interpolate() requires a bare sfc, not a full sf data frame. st_geometry() pulls just the geometry list-column out.
    st_geometry(pts)[empty] <- st_line_interpolate(st_geometry(short_edges), half_len)
  }

  pts <- st_cast(pts, "POINT")
  pts$predID <- paste0("P", seq_len(nrow(pts)))

  ## This estimates drainage area at each synthetic point needed for area-dependent predictions in fucture scripts. This is approximated the same way flow accumulation was sampled onto streams_final$flowaccum above as raster cell count x cell area, rather than a full per-point watershed delineation (wbt_watershed() per site), which is too slow for hundreds of densely-spaced synthetic points. A small buffer + max() makes this robust to a sampled point landing a pixel or two off the rasterized stream cell it's meant to represent.
  faccum_r <- raster(faccum_path)
  cellsize <- mean(raster::res(faccum_r))
  flowaccum <- raster::extract(faccum_r, st_coordinates(pts),
                                buffer = 2 * cellsize, fun = max, na.rm = TRUE)
  ## This guards against any point still coming back NA/non-positive (e.g. a genuine edge case right at the raster boundary). fall back to the smallest valid drainage area found elsewhere rather than letting a single bad point produce a -Inf/NaN logArea that breaks ssn_lm()/ predict() for the whole prediction dataset.
  bad <- is.na(flowaccum) | flowaccum <= 0
  if (any(bad)) {
    fallback <- suppressWarnings(min(flowaccum[!bad], na.rm = TRUE))
    if (!is.finite(fallback)) fallback <- 1  # every point came back bad -- last resort
    cat("    ", sum(bad), "prediction point(s) got no valid flow accumulation --",
        "using the smallest valid value found (", round(fallback), "cells) instead.\n")
    flowaccum[bad] <- fallback
  }
  pts$Area_m2 <- flowaccum * cell_area_m2
  pts$logArea <- log(pts$Area_m2)

  pts
}

# this applies the prediction function - takes a moment to run
preds <- imap(edges, function(e, name) {
  make_pred_points(e, pred_spacing_m, geo_results[[name]]$faccum_path,
                    geo_results[[name]]$cell_area_m2)
})

# this saves the predictions to a landscape network 
preds <- imap(preds, function(p, name) {
  sites_to_lsn(
    sites = p,
    edges = edges[[name]],
    lsn_path = networks[[name]]$lsn_path,
    snap_tolerance = 1,  # these points were generated ON the edges, so a tight tolerance is fine/expected here, unlike the generous 100m used for real site coords
    save_local = TRUE,
    file_name = "predpts.gpkg",
    overwrite = TRUE
  )
})

# this maps on afv
preds <- imap(preds, function(p, name) {
  afv_sites(
    sites = list(preds = p),
    edges = edges[[name]],
    afv_col = "afvFlow",
    lsn_path = networks[[name]]$lsn_path,
    save_local = TRUE,
    overwrite = TRUE
  )$preds
})

# this calculates upstream distances 
preds <- imap(preds, function(p, name) {
  updist_sites(
    sites = list(preds = p),
    edges = edges[[name]],
    length_col = "Length_m",
    lsn_path = networks[[name]]$lsn_path,
    save_local = TRUE,
    overwrite = TRUE
  )$preds
})

# this tells you how many synthetic points were produced given the spacing specified
print(map_int(preds, nrow))

## B5. Assemble the final .ssn objects from both real and synthetic datasets ----

ssn_objects <- imap(edges, function(e, name) {
  ssn_assemble(
    edges = e,
    lsn_path = networks[[name]]$lsn_path,
    obs_sites = obs[[name]],
    preds_list = list(predpts = preds[[name]]),
    ssn_path = networks[[name]]$ssn_path,
    import = TRUE,
    check = TRUE,
    afv_col = "afvFlow",
    overwrite = TRUE
  )
})

summary(ssn_objects$nm)

## B6. Build distance matrices ----
## Distance matrices among real and synthetic points are needed in downstream scripts that use tail-up/tail-down covariance, but ssn_assemble() above does not build them automatically, and they aren't restored by readRDS() either (they're written as files under each .ssn folder's distance/ subdirectory). This sections builds them once here for use in downstream scripts.

walk(ssn_objects, function(obj) {
  ssn_create_distmat(ssn.object = obj, predpts = "predpts",
                     among_predpts = TRUE, overwrite = TRUE)
})

## B7. Verify network topology routing ----
## Verifying topology/routing here to catch a broken/misrouted network before further use. Note that a broken/misrouted network can otherwise look like a genuine absence of spatial structure in future scripts, so it is critical to carefully verify it here first. 
## Three checks:
# 1. How many distinct NetworkIDs exist? Should be 1 per stream network unless you know of a real second, disconnected drainage in that watershed.
# 2. DistanceUpstream per site, sorted -- should be ~0 at the outlet and largest at the furthest headwater. Eyeball whether this matches what you know about each site's real position in the network.
# 3. check_pair(): pick two sites you already know the real hydrologic relationship for (e.g. the outlet and a known upstream site) and verify SSN2 agrees, with the actual routed distance printed so you can sanity-check its magnitude against a real map.
# 4. Visual check: two maps, one colored by NetworkID (fragmentation check), one by upstream distance (routing gut check)

## 1. How many distinct NetworkIDs exist? &
## 2. DistanceUpstream per site, sorted
# define check function
inspect_ssn_topology <- function(ssn.object, proj) {
  obs_geom <- ssn_get_netgeom(ssn.object$obs, reformat = TRUE)
  obs_geom$Site <- ssn.object$obs$Site

  cat("\n", proj, "-- distinct NetworkID(s) among obs sites:",
      paste(sort(unique(obs_geom$NetworkID)), collapse = ", "), "\n")
  if (length(unique(obs_geom$NetworkID)) > 1) {
    cat("  ! more than one NetworkID -- sites on different NetworkIDs are",
        "not flow-connected OR flow-unconnected to each other (SSN2 treats",
        "them as unrelated networks).\n")
  }

  cat(proj, "-- DistanceUpstream by site (~0 expected at the outlet,",
      "largest at the furthest headwater):\n")
  print(obs_geom[order(obs_geom$DistanceUpstream),
                 c("Site", "NetworkID", "DistanceUpstream")])

  obs_geom
}
# run check
ssn_topology <- imap(ssn_objects, inspect_ssn_topology)

## 3. check_pair(): 
# Per ssn_get_stream_distmat()'s convention, columns are the FROM site and rows are the TO site; for a flow-connected pair the downstream-only distance is > 0 in exactly one direction and 0 in the other, while a flow-unconnected pair (shares a downstream confluence but no direct flow) is > 0 in BOTH directions.
# define check function
check_pair <- function(proj, site_a, site_b) {
  geom <- ssn_topology[[proj]]
  net_a <- geom$NetworkID[geom$Site == site_a]
  net_b <- geom$NetworkID[geom$Site == site_b]
  if (length(net_a) == 0 || length(net_b) == 0) {
    cat("Couldn't find", site_a, "and/or", site_b, "in", proj, "-- check spelling.\n")
    return(invisible(NULL))
  }
  if (net_a != net_b) {
    cat(site_a, "and", site_b, "are on DIFFERENT NetworkIDs (", net_a, "vs", net_b,
        ") -- not comparable in this .ssn object at all.\n")
    return(invisible(NULL))
  }
  dist_obs <- ssn_get_stream_distmat(ssn_objects[[proj]])
  dmat <- dist_obs[[paste0("dist.net", net_a)]]
  pid_a <- as.character(geom$pid[geom$Site == site_a])
  pid_b <- as.character(geom$pid[geom$Site == site_b])
  d_ab <- dmat[pid_a, pid_b]
  d_ba <- dmat[pid_b, pid_a]
  flow_connected <- (min(d_ab, d_ba) == 0) && (max(d_ab, d_ba) > 0)
  cat(site_a, "->", site_b, ":", round(d_ab, 1), "m  |  ",
      site_b, "->", site_a, ":", round(d_ba, 1), "m  |  ",
      "flow-connected:", flow_connected,
      if (!flow_connected && max(d_ab, d_ba) > 0) " (flow-unconnected -- shares a downstream confluence, no direct flow)" else "",
      "\n")
}
# run check (Example: replace with a pair you already know the real relationship)
check_pair("nm", "USF12", "USF25")

## 4. Visual check: two maps, one colored by NetworkID, one by upstream distance. Compare against what you know of the real network shape and site positions.
walk2(ssn_objects, names(ssn_objects), function(ssn.object, proj) {
    edges_geom <- ssn_get_netgeom(ssn.object$edges, reformat = TRUE)
    ssn.object$edges$NetworkID <- edges_geom$NetworkID
    ssn.object$edges$upDist <- edges_geom$DistanceUpstream

    cat("\n", proj, "-- topology maps (check these against the real network shape):\n")
    print(mapview(ssn.object$edges, zcol = "NetworkID", layer.name = paste(proj, "NetworkID")) +
            mapview(ssn.object$obs, color = "black", col.regions = "yellow",
                     legend = FALSE, cex = 5, layer.name = paste(proj, "sites")))
    print(mapview(ssn.object$edges, zcol = "upDist", layer.name = paste(proj, "upstream distance (m)")) +
            mapview(ssn.object$obs, color = "black", col.regions = "yellow",
                     legend = FALSE, cex = 5, layer.name = paste(proj, "sites")))
  })

#### PART C -- Save results and clear temp data ####
## Save ssn objects so downstream scripts can either re-import from ssn_path or just readRDS() this directly.
saveRDS(ssn_objects, file.path(data_out_dir, "ssn_objects.rds"))

## Run after script 1 completes - deletes everything in temp/ 
unlink("temp", recursive = TRUE)
