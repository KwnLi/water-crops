library(reticulate)
library(sf)
library(terra)
# remotes::install_github("ClimateEcology/optimEcoServices")
# remotes::install_github("ClimateEcology/rPCSE")
library(rPCSE)
library(optimEcoServices)

reticulate::use_condaenv("py3_pcse")
pcse <- reticulate::import("pcse")

datadir <- "/Users/kevinl/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/GIS/maxWFS/"

# objective: identify all the mukey and grid ID combinations covered by corn and soy fields
optimEcoServices::biophys_pol |> dplyr::filter(CLASS_NAME %in% c("Corn","Soybeans"))
# Corn = 1
# Soybeans = 5

# make the agromanagement file for corn
make_agromanager(outpath = "data", fname = "corn",  # define output location and filename
                 years = c(2020,2021,2022),         # model run years
                 crop_name = "maize",               # crop name, must match available crop format
                 variety_name = "Grain_maize_201",  # crop variety name, must match available crops
                 crop_start_month = 5,              # start month (%m)
                 crop_start_day = 17,               # start day (%d)
                 crop_start_type = "sowing",        # crop start type
                 crop_end_month = 11,               # end month (%m)
                 crop_end_day = 2,                  # end day (%d)
                 crop_end_type = "maturity",        # crop end type
                 max_duration = 300                 # maximum duration crop model can run
)
agro.corn <- pcse$input$YAMLAgroManagementReader("data/corn.agro")

# make the agromanagement file for soybeans
make_agromanager(outpath = "data", fname = "soy",   # define output location and filename
                 years = c(2020,2021,2022),         # model run years
                 crop_name = "soybean",             # crop name, must match available crop format
                 variety_name = "Soybean_901",      # crop variety name, must match available crops
                 crop_start_month = 5,              # start month (%m)
                 crop_start_day = 30,               # start day (%d)
                 crop_start_type = "sowing",        # crop start type
                 crop_end_month = 10,               # end month (%m)
                 crop_end_day = 30,                 # end day (%d)
                 crop_end_type = "maturity",        # crop end type
                 max_duration = 300                 # maximum duration crop model can run
)
agro.soy <- pcse$input$YAMLAgroManagementReader("data/soy.agro")

# site parameters and data
siteparam <- pcse$input$WOFOST72SiteDataProvider(WAV=10)

hupa <- readRDS("data/pa-ws.rds")
gridIDs <- readRDS("data/gridIDs.rds")

# crop data provider
cropparam <- pcse$input$YAMLCropDataProvider()

# containers for outputs
wofost.out <- list()   # wofost results
cornsoy.id <- list()  # mukey ids of results + cell counts
for(i in seq_len(nrow(hupa))){
  hu.i <- hupa$huc12[i]
  hudir.i <- paste0(datadir,"hu/","hu",hu.i,"/data/")

  lu.i <- terra::rast(paste0(hudir.i,"lu/fields.tif"))
  lu.cornsoy <- terra::ifel(lu.i[["crop"]] %in% c(1,5),lu.i[["crop"]],NA)

  hu.mukey <- terra::rast(paste0("data-raw/wofost/soil/mukey_",hu.i,".tif")) |>
    terra::project(crs(lu.cornsoy), method = "near") |>
    terra::resample(lu.cornsoy) |>
    terra::crop(lu.cornsoy) |> terra::mask(lu.cornsoy)

  lu.cornsoy$mukey <- terra::ifel(!is.na(lu.cornsoy), hu.mukey, NA)

  hu.cornsoy <- terra::values(lu.cornsoy) |> as.data.frame() |>
    dplyr::left_join(gridIDs[[hu.i]]) |>
    dplyr::group_by(crop,mukey,grid.ids) |> dplyr::summarize(n=dplyr::n(),.groups = "drop") |>
    dplyr::filter(!is.na(crop))

  # run corn simulations in WOFOST
  hu.corn <- hu.cornsoy |> dplyr::filter(crop==1)
  hu.corn.wofost <- list()
  for(g in seq_len(nrow(hu.corn))){
    hu.corn.g <- hu.corn[g,]

    # load WOFOST soil file
    cornsoil <- pcse$input$PCSEFileReader(paste0("data/soil/mukey",hu.corn.g$mukey,".pcse"))

    # load weather file
    cornweather <- pcse$input$CSVWeatherDataProvider(
      here::here(paste0("data/weather/id",hu.corn.g$grid.ids,".csv")))

    # Bundle parameters
    cornparam <- pcse$base$ParameterProvider(sitedata=siteparam, soildata = cornsoil, cropdata = cropparam)

    # Load parameters into the model engine
    cornwofost <- pcse$models$Wofost72_WLP_CWB(cornparam, cornweather, agro.corn)

    # run WOFOST
    cornwofost$run_till_terminate()

    # get output
    corn.output <- cornwofost$get_output()

    hu.corn.wofost[[as.character(hu.corn.g$mukey)]] <- lapply(corn.output, dplyr::bind_cols) |> dplyr::bind_rows()
  }

  # run soy simulations in WOFOST
  hu.soy <- hu.cornsoy |> dplyr::filter(crop==5)
  hu.soy.wofost <- list()
  for(h in seq_len(nrow(hu.soy))){
    hu.soy.h <- hu.soy[h,]

    # load WOFOST soil file
    soysoil <- pcse$input$PCSEFileReader(paste0("data/soil/mukey",hu.soy.h$mukey,".pcse"))

    # load weather file
    soyweather <- pcse$input$CSVWeatherDataProvider(
      here::here(paste0("data/weather/id",hu.soy.h$grid.ids,".csv")))

    # Bundle parameters
    soyparam <- pcse$base$ParameterProvider(sitedata=siteparam, soildata = soysoil, cropdata = cropparam)

    # Load parameters into the model engine
    soywofost <- pcse$models$Wofost72_WLP_CWB(soyparam, soyweather, agro.soy)

    # run WOFOST
    soywofost$run_till_terminate()

    # get output
    soy.output <- soywofost$get_output()

    hu.soy.wofost[[as.character(hu.soy.h$mukey)]] <- lapply(soy.output, dplyr::bind_cols) |> dplyr::bind_rows()
  }

  # combine results for hu.i
  wofost.out[[hu.i]] <- dplyr::bind_rows(hu.corn.wofost, .id = "mukey") |>
    dplyr::mutate(crop = "corn") |>
    dplyr::bind_rows(
      dplyr::bind_rows(hu.soy.wofost, .id = "mukey") |>
        dplyr::mutate(crop = "soybeans")
    )

  # record hu.cornsoy
  cornsoy.id[[hu.i]] <- hu.cornsoy

}

# make data into dataframe
wofost.out <- dplyr::bind_rows(wofost.out, .id = "huc12")
cornsoy.id <- dplyr::bind_rows(cornsoy.id, .id = "huc12")
