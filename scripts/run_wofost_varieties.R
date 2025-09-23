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

# identify all the mukey and grid ID combinations covered by corn and soy fields
optimEcoServices::biophys_pol |> dplyr::filter(CLASS_NAME %in% c("Corn","Soybeans"))
# Corn = 1
# Soybeans = 5

# make the agromanagement file for corn and soy varieties
dir.create("data/wofost_varieties")

# corn varieties
corn_var <- c("Grain_maize_201",
              "Grain_maize_202", "Grain_maize_203",
              "Grain_maize_204", "Grain_maize_205",
              "Fodder_maize_nl", "Maize_VanHeemst_1988")

for(i in seq_along(corn_var)){
  cornvar.i <- corn_var[i]
  make_agromanager(outpath = "data/wofost_varieties",
                   fname = corn_var[i],               # define output location and filename
                   years = c(2020,2021,2022),         # model run years
                   crop_name = "maize",               # crop name, must match available crop format
                   variety_name = corn_var[i],        # crop variety name, must match available crops
                   crop_start_month = 5,              # start month (%m)
                   crop_start_day = 17,               # start day (%d)
                   crop_start_type = "sowing",        # crop start type
                   crop_end_month = 11,               # end month (%m)
                   crop_end_day = 2,                  # end day (%d)
                   crop_end_type = "maturity",        # crop end type
                   max_duration = 300                 # maximum duration crop model can run
  )
}

# make the agromanagement file for soybeans
soy_var <- c("Soybean_901", "Soybean_902", "Soybean_903", "Soybean_904",
             "Soybean_905", "Soybean_906", "Soybean_VanHeemst_1988")

for(i in seq_along(soy_var)){
  make_agromanager(outpath = "data/wofost_varieties",
                   fname = soy_var[i],   # define output location and filename
                   years = c(2020,2021,2022),         # model run years
                   crop_name = "soybean",             # crop name, must match available crop format
                   variety_name = soy_var[i],      # crop variety name, must match available crops
                   crop_start_month = 5,              # start month (%m)
                   crop_start_day = 30,               # start day (%d)
                   crop_start_type = "sowing",        # crop start type
                   crop_end_month = 10,               # end month (%m)
                   crop_end_day = 30,                 # end day (%d)
                   crop_end_type = "maturity",        # crop end type
                   max_duration = 300                 # maximum duration crop model can run
  )
}


# site parameters and data
siteparam <- pcse$input$WOFOST72SiteDataProvider(WAV=10)

hupa <- readRDS("data/pa-ws.rds")
# gridIDs <- readRDS("data/gridIDs.rds")

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
    terra::resample(lu.cornsoy)
    # terra::crop(lu.cornsoy) |> terra::mask(lu.cornsoy)

  gridmet.id <- terra::rast(paste0("data-raw/wofost/weather/weatherID_",hu.i,".tif")) |>
    terra::resample(lu.cornsoy |>
                      terra::project(crs(terra::rast(paste0("data-raw/wofost/weather/weatherID_",hu.i,".tif")))),
                    method = "near") |>  # match resolution
    terra::project(crs(lu.cornsoy), method = "near") |>
    terra::resample(lu.cornsoy, method = "near")

  lu.cornsoy$mukey <- terra::ifel(!is.na(lu.cornsoy), hu.mukey, NA)
  lu.cornsoy$grid.ids <- terra::ifel(!is.na(lu.cornsoy), gridmet.id, NA)

  hu.cornsoy <- terra::values(lu.cornsoy) |> as.data.frame() |>
    dplyr::group_by(crop,mukey,grid.ids) |> dplyr::summarize(n=dplyr::n(),.groups = "drop") |>
    dplyr::filter(!is.na(crop))

  # run corn simulations in WOFOST
  hu.corn <- hu.cornsoy |> dplyr::filter(crop==1)
  hu.corn.wofost <- list()
  for(g in seq_len(nrow(hu.corn))){
    hu.corn.g <- hu.corn[g,]

    # load WOFOST soil file
    cornsoil.file <- paste0("data/soil/mukey",hu.corn.g$mukey,".pcse")
    if(!file.exists(cornsoil.file)) next  # if this file doesn't exist, skip

    tryCatch({
      cornsoil <- pcse$input$PCSEFileReader(cornsoil.file)
    }, error = function(e){cat("ERROR :", hu.i, hu.corn.g$mukey, conditionMessage(e), "\n")})

    # load weather file
    cornweather <- pcse$input$CSVWeatherDataProvider(
      here::here(paste0("data/weather/id",hu.corn.g$grid.ids,".csv")))

    # Bundle parameters
    cornparam <- pcse$base$ParameterProvider(sitedata=siteparam, soildata = cornsoil, cropdata = cropparam)

    # loop through corn varieties
    hu.corn.var <- list()
    for(x in seq_along(corn_var)){
      agro.corn.x <- pcse$input$YAMLAgroManagementReader(paste0("data/wofost_varieties/",corn_var[x],".agro"))
      # Load parameters into the model engine
      cornwofost <- pcse$models$Wofost72_WLP_CWB(cornparam, cornweather, agro.corn.x)

      # potential production
      cornpp <- pcse$models$Wofost72_PP(cornparam, cornweather, agro.corn.x)

      # run WOFOST
      cornwofost$run_till_terminate()
      cornpp$run_till_terminate()

      # get output
      corn.output <- cornwofost$get_output()
      cornpp.output <- cornpp$get_output()

      # record output
      hu.corn.var[[corn_var[x]]] <- lapply(corn.output, dplyr::bind_cols) |>
        dplyr::bind_rows() |>
        dplyr::left_join(
          lapply(cornpp.output, dplyr::bind_cols) |> dplyr::bind_rows(),
          by = "day", suffix = c(".WL", ".PP")
        )
    }

    hu.corn.wofost[[paste0(hu.corn.g$mukey,"_",hu.corn.g$grid.ids)]] <- dplyr::bind_rows(hu.corn.var, .id = "variety")
  }

  # run soy simulations in WOFOST
  hu.soy <- hu.cornsoy |> dplyr::filter(crop==5)
  hu.soy.wofost <- list()
  for(h in seq_len(nrow(hu.soy))){
    hu.soy.h <- hu.soy[h,]

    # load WOFOST soil file
    soysoil.file <- paste0("data/soil/mukey",hu.soy.h$mukey,".pcse")
    if(!file.exists(soysoil.file)) next  # if this file doesn't exist, skip

    tryCatch({
      soysoil <- pcse$input$PCSEFileReader(soysoil.file)
    }, error = function(e){cat("ERROR :", hu.i, hu.soy.h$mukey, conditionMessage(e), "\n")})


    # load weather file
    soyweather <- pcse$input$CSVWeatherDataProvider(
      here::here(paste0("data/weather/id",hu.soy.h$grid.ids,".csv")))

    # Bundle parameters
    soyparam <- pcse$base$ParameterProvider(sitedata=siteparam, soildata = soysoil, cropdata = cropparam)

    # loop through soy varieties
    hu.soy.var <- list()
    for(y in seq_along(soy_var)){
      agro.soy.y <- pcse$input$YAMLAgroManagementReader(paste0("data/wofost_varieties/", soy_var[y], ".agro"))

      # Load parameters into the model engine
      soywofost <- pcse$models$Wofost72_WLP_CWB(soyparam, soyweather, agro.soy.y)

      # potential production
      soypp <- pcse$models$Wofost72_PP(soyparam, soyweather, agro.soy.y)

      # run WOFOST
      soywofost$run_till_terminate()
      soypp$run_till_terminate()

      # get output
      soy.output <- soywofost$get_output()
      soypp.output <- soypp$get_output()

      # record output
      hu.soy.var[[soy_var[y]]] <- lapply(soy.output, dplyr::bind_cols) |>
        dplyr::bind_rows() |>
        dplyr::left_join(
          lapply(soypp.output, dplyr::bind_cols) |> dplyr::bind_rows(),
          by = "day", suffix = c(".WL", ".PP")
        )
    }
    hu.soy.wofost[[paste0(hu.soy.h$mukey,"_",hu.soy.h$grid.ids)]] <- dplyr::bind_rows(hu.soy.var, .id = "variety")
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
wofost.outdf <- dplyr::bind_rows(wofost.out, .id = "huc12") |> dplyr::rename(mukey_gridid = mukey)
cornsoy.iddf <- dplyr::bind_rows(cornsoy.id, .id = "huc12")

saveRDS(wofost.outdf |> dplyr::filter(crop=="corn"), "PA_wofost_corn_var.rds")
saveRDS(wofost.outdf |> dplyr::filter(crop=="soybeans"), "PA_wofost_soy_var.rds")
saveRDS(cornsoy.iddf, "PA_cornsoyIDs_var.rds")
