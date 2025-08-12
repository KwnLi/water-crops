library(reticulate)
library(sf)
library(terra)
# gitdir <- "/Users/kevinl/Documents/GitHub/"
# devtools::load_all(paste0(gitdir,"rPCSE"))
# remotes::install_github("ClimateEcology/rPCSE")
library(rPCSE)

reticulate::use_condaenv("py3_pcse")
pcse <- reticulate::import("pcse")

#### load watershed data ####
datadir <- "/Users/kevinl/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/GIS/maxWFS/hu/"

hupa <- readRDS("data/pa-ws.rds")

#### wofost raw data dir ####
soilraw <- "data-raw/wofost/soil/"
weatherraw <- "data-raw/wofost/weather/"

# store grid ID combos
dir.create("data/gridIDs")
pa.gridIDs <- list()

for(i in seq_len(nrow(hupa))){
  hu.i <- hupa$huc12[i]
  datadir.i <- paste0(datadir,"hu",hu.i,"/data/")

  # watershed outline
  hu.vct <- terra::vect(paste0(datadir.i,"watershed.gpkg"))

  # land use raster to use as basis of projections
  lu.i <- terra::rast(paste0(datadir.i,"lu/fields.tif"))

  # soil data
  gssurgo.mukey <- terra::rast(paste0(soilraw,"mukey_",hu.i,".tif")) |>
    terra::project(crs(lu.i), method = "near") |>
    terra::resample(lu.i) |>  # match pixels
    terra::mask(hu.vct) # mask with watershed shape
  gssurgohz.dat <- read.csv(paste0(soilraw,"ssurgohz_",hu.i,".csv"))
  gssurgosub.dat <- read.csv(paste0(soilraw,"ssurgosub_",hu.i,".csv"))
  gssurgovalu1.dat <- read.csv(paste0(soilraw,"ssurgovalu1_",hu.i,".csv"))

  mukeysoil.dat <- dplyr::left_join(gssurgohz.dat, gssurgovalu1.dat)

  # weather data
  gridmet.id <- terra::rast(paste0(weatherraw,"weatherID_",hu.i,".tif")) |>
    terra::resample(gssurgo.mukey |>
                      terra::project(crs(terra::rast(paste0(weatherraw,"weatherID_",hu.i,".tif")))),
                    method = "near") |>  # match resolution
    terra::project(crs(gssurgo.mukey), method = "near") |>
    terra::resample(gssurgo.mukey, method = "near") |>
    terra::crop(gssurgo.mukey) |> # match extent
    terra::mask(hu.vct)
  gridmet.dat <- read.csv(paste0(weatherraw,"gridmet_",hu.i,".csv"))

  # elevation
  elev <- terra::rast(paste0(datadir.i,"elev.tif")) |>
    terra::project(crs(gssurgo.mukey), method = "near") |>
    terra::resample(gssurgo.mukey) |>  # match pixels
    terra::crop(gssurgo.mukey) # match extent

  #### match IDs ####

  # # OLD find unique combinations of mukey and weather id (use the most common gridID for each mukey)
  # gridIDs <- terra::zonal(gridmet.id, gssurgo.mukey, fun = rPCSE::Mode) |>
  #   dplyr::distinct() |>
  #   dplyr::rename(mukey = "gssurgocdl", grid.ids = "id")

  mukeygridid <- gssurgo.mukey
  mukeygridid$grid.ids <- terra::ifel(!is.na(gssurgo.mukey), gridmet.id, NA)

  # find unique combinations of mukey and weather id
  gridIDs <- terra::values(mukeygridid) |> as.data.frame() |>
    dplyr::group_by(gssurgocdl,grid.ids) |>
    dplyr::summarize(n=dplyr::n(),.groups = "drop") |>
    dplyr::rename(mukey = "gssurgocdl")

  #### weather data processing ####

  # find mean elevation of gridMET
  gridMET.elev <- terra::zonal(elev, gridmet.id, fun = "mean", na.rm = TRUE) |>
    dplyr::rename(ELEV = "elev3dep")

  # find mean longitude of gridMET
  gridMET.wgs84.lon <- terra::project(gridmet.id, "epsg:4326") |>
    terra::init("x") |> # calculate longitude
    terra::project(crs(gridmet.id), method = "near") |> # project back
    terra::resample(gridmet.id, method = "near")

  lon <- terra::zonal(gridMET.wgs84.lon, gridmet.id, fun = "mean") |>
    dplyr::rename(LON = "id.1")

  # find mean latitude of gridMET
  gridMET.wgs84.lat <- terra::project(gridmet.id, "epsg:4326") |>
    terra::init("y") |> # calculate latitude
    terra::project(crs(gridmet.id), method = "near") |> # project back
    terra::resample(gridmet.id, method = "near")

  lat <- terra::zonal(gridMET.wgs84.lat, gridmet.id, fun = "mean") |>
    dplyr::rename(LAT = "id.1")

  # write.csv(dplyr::left_join(lon,lat),"gridMETlatlon.csv",row.names = FALSE)

  weather <- gridmet.dat |> dplyr::filter(grid.ids %in% gridIDs$grid.ids) |>
    dplyr::mutate(rh = (rmin+rmax)/2, temp = (tmin+tmax)/2,
                  VAP = vap_from_relhum(rh = rh, temp = temp),  # calculate relative humidity
                  IRRAD = srad*60*60*24/1000,  # convert from srad units (W/m2) to IRRAD units (kJ/m2/day)
                  WIND = wind10to2(vs),
                  DAY = as.numeric(gsub("-","",date)),
                  RAIN = pr,
                  SNOWDEPTH = NaN) |>
    dplyr::rename(TMIN = tmin, TMAX = tmax) |>
    dplyr::select(grid.ids, DAY, IRRAD, TMIN, TMAX, VAP, WIND, RAIN, SNOWDEPTH) %>% # this dplyr pipe is necessary for split!
    split(.$grid.ids)

  # loop through all the grid IDs and make weather csvs
  if(!dir.exists("data/weather")) dir.create("data/weather", recursive = TRUE)
  for(j in seq_len(length(weather))){
    grid.j <- names(weather)[j]
    lat.j <- lat[lat$id == as.numeric(grid.j),"LAT"]
    lon.j <- lon[lon$id == as.numeric(grid.j),"LON"]
    elev.j <- gridMET.elev[gridMET.elev$id == as.numeric(grid.j),"ELEV"]
    weather.dat <- weather[[grid.j]]

    fname.j <- paste0("data/weather/id",grid.j)

    make_weather_csv(fname = fname.j,
                     dat=weather.dat, lon = lon.j, lat = lat.j, elev = elev.j,
                     country = "USA", station = grid.j, srce = "gridMET")

  }

  #### soil data processing ####
  if(!dir.exists("data/soil")) dir.create("data/soil", recursive = TRUE)

  nonzero.rootzone <- gridIDs |> dplyr::left_join(gssurgovalu1.dat) |>
    dplyr::filter(rootznemc > 0)

  for(k in seq_len(nrow(nonzero.rootzone))){
    mukey.k <- nonzero.rootzone$mukey[k]

    fname.k <- paste0("mukey",mukey.k)

    ssurgodat.k <- gssurgohz.dat |> dplyr::filter(mukey == mukey.k) |>
      dplyr::rename(ksat_r_rootzone = ksat_r) |>
      dplyr::left_join(gssurgosub.dat, by = "mukey") |>
      dplyr::rename(ksat_r_subsoil = ksat_r) |>
      dplyr::left_join(gssurgovalu1.dat, by = "mukey") |>
      as.list()

    make_soil_pcse(outpath = "data/soil", fname = fname.k,ssurgo.dat=ssurgodat.k)
  }

  # write out gridIDs
  write.csv(gridIDs, paste0("data/gridIDs/gridID_hu", hu.i,".csv"), row.names = FALSE)
  pa.gridIDs[[hu.i]] <- gridIDs
}


saveRDS(pa.gridIDs, "data/gridIDs.rds")



