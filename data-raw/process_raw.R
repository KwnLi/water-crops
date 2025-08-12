library(reticulate)
library(sf)
library(terra)
devtools::load_all()

reticulate::use_condaenv("py3_pcse")
pcse <- reticulate::import("pcse")

#### load watershed data ####
datadir <- "/Users/kevinl/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/GIS/maxWFS/hu/hu070200080301/data/"
hu <- "070200080301"

# watershed outline
hu.vct <- terra::vect(paste0(datadir,"watershed.gpkg"))

# soil data
gssurgo.mukey <- terra::rast("data-raw/wofost/soil/mukey_070200080301.tif") |>
  terra::project(crs(hu.vct), method = "near") |>
  terra::mask(hu.vct) # mask with watershed shape
gssurgohz.dat <- read.csv("data-raw/wofost/soil/ssurgohz_070200080301.csv")
gssurgosub.dat <- read.csv("data-raw/wofost/soil/ssurgosub_070200080301.csv")
gssurgovalu1.dat <- read.csv("data-raw/wofost/soil/ssurgovalu1_070200080301.csv")

mukeysoil.dat <- dplyr::left_join(gssurgohz.dat, gssurgovalu1.dat)

# weather data
gridmet.id <- terra::rast("data-raw/wofost/weather/weatherID_070200080301.tif") |>
  terra::resample(gssurgo.mukey |>
                    terra::project(crs(terra::rast("data-raw/wofost/weather/weatherID_070200080301.tif"))),
                  method = "near") |>  # match resolution
  terra::project(crs(hu.vct), method = "near") |>
  terra::resample(gssurgo.mukey, method = "near") |>
  terra::crop(gssurgo.mukey) |> # match extent
  terra::mask(hu.vct)
gridmet.dat <- read.csv("data-raw/wofost/weather/gridmet_070200080301.csv")

# elevation
elev <- terra::rast(paste0(datadir,"elev.tif")) |>
  terra::project(crs(hu.vct), method = "near") |>
  terra::crop(gssurgo.mukey) # match extent

#### match IDs ####

# find unique combinations of mukey and weather id (use the most common gridID for each mukey)
gridIDs <- terra::zonal(gridmet.id, gssurgo.mukey, fun = Mode) |>
  dplyr::distinct() |>
  dplyr::rename(mukey = "gssurgocdl", grid.ids = "id")

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
dir.create("data/weather", recursive = TRUE)
for(i in seq_len(length(weather))){
  grid.i <- names(weather)[i]
  lat.i <- lat[lat$id == as.numeric(grid.i),"LAT"]
  lon.i <- lon[lon$id == as.numeric(grid.i),"LON"]
  elev.i <- gridMET.elev[gridMET.elev$id == as.numeric(grid.i),"ELEV"]
  weather.dat <- weather[[grid.i]]

  fname.i <- paste0("data/weather/id",grid.i)

  make_weather_csv(fname = fname.i,
                   dat=weather.dat, lon = lon.i, lat = lat.i, elev = elev.i,
                   country = "USA", station = grid.i, srce = "gridMET")

}

#### soil data processing ####
dir.create("data/soil", recursive = TRUE)

nonzero.rootzone <- gridIDs |> dplyr::left_join(gssurgovalu1.dat) |>
  dplyr::filter(rootznemc > 0)

for(i in seq_len(nrow(nonzero.rootzone))){
  mukey.i <- nonzero.rootzone$mukey[i]

  fname.i <- paste0("mukey",mukey.i)

  ssurgodat.i <- gssurgohz.dat |> dplyr::filter(mukey == mukey.i) |>
    dplyr::rename(ksat_r_rootzone = ksat_r) |>
    dplyr::left_join(gssurgosub.dat, by = "mukey") |>
    dplyr::rename(ksat_r_subsoil = ksat_r) |>
    dplyr::left_join(gssurgovalu1.dat, by = "mukey") |>
    as.list()

  make_soil_pcse(outpath = "data/soil", fname = fname.i,ssurgo.dat=ssurgodat.i)
}
