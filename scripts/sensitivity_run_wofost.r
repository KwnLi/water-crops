library(reticulate)
# remotes::install_github("ClimateEcology/optimEcoServices")
# remotes::install_github("ClimateEcology/rPCSE")
library(rPCSE)

temp.weatherfile <- "data/sensitivityweather.csv"   # this is the name of the temp weather file that will be rewritten over the loops

soils <- Rwofost::wofost_soil() |> strsplit("\n|,\\s") |> unlist()

reticulate::use_condaenv("py3_pcse")
pcse <- reticulate::import("pcse")

# crop data provider
cropparam <- pcse$input$YAMLCropDataProvider()

# load corn variety
agro.corn <- pcse$input$YAMLAgroManagementReader("data/wofost_varieties/Grain_maize_201.agro")

# load soy variety
agro.soy <- pcse$input$YAMLAgroManagementReader("data/wofost_varieties/Soybean_901.agro")

# load in weathergrid
weathergrid <- readRDS("data/weathergrid.rds")
weathermns <- read.csv("data/weathergrid.csv")

# sensitivity parameters
# these are strings representing the numeric value of the 'bumps'
deltaT <- names(weathergrid)
deltaRAIN <- names(weathergrid[[1]])

# load middlemost grid to get the header data
midgrid <- weathermns$gridid[which(weathermns$distCenter==min(weathermns$distCenter))]
midgrid.head <- read.csv(paste0("data/weather/",midgrid),nrows=8, header = FALSE)

##### LOOP #####
corn.sens <- list()
soy.sens <- list()

for(g in seq_along(soils)){
  soilfile.g <- paste0("data/Rwofost_soil/",soils[g],".pcse")

  # load WOFOST soil file
  soil.g <- pcse$input$PCSEFileReader(soilfile.g)

  # site parameters and data
  siteparam.g <- pcse$input$WOFOST72SiteDataProvider(WAV=wofost_soil(soils[g])$WAV)

  # crop output for soil g
  corn.sens.g <- list()
  soy.sens.g <- list()

  # loop through delta T and delta Rain
  for(i in seq_along(deltaT)){
    deltaT.i <- deltaT[i]

    for(j in seq_along(deltaRAIN)){
      deltaRAIN.j <- deltaRAIN[j]

      weather.ij <- weathergrid[[deltaT.i]][[deltaRAIN.j]]

      # write out a temp weatherfile
      writeLines(midgrid.head[,1], temp.weatherfile)
      write.table(weather.ij, temp.weatherfile, sep = ",", eol = "\n", fileEncoding = "UTF-8",
                  na = "NaN",quote = FALSE,
                  col.names = TRUE, row.names = FALSE, append=TRUE)

      # load weatherfile
      pcseweather.ij <- pcse$input$CSVWeatherDataProvider(
        here::here(temp.weatherfile))

      # Bundle parameters
      param.ij <- pcse$base$ParameterProvider(sitedata=siteparam.g,
                                              soildata = soil.g,
                                              cropdata = cropparam)

      # Load parameters into the model engine
      # corn model
      cornwofost <- pcse$models$Wofost72_WLP_CWB(param.ij, pcseweather.ij, agro.corn)
      # corn potential production
      cornpp <- pcse$models$Wofost72_PP(param.ij, pcseweather.ij, agro.corn)
      # soy model
      soywofost <- pcse$models$Wofost72_WLP_CWB(param.ij, pcseweather.ij, agro.soy)
      # soy potential production
      soypp <- pcse$models$Wofost72_PP(param.ij, pcseweather.ij, agro.soy)

      # run WOFOST
      cornwofost$run_till_terminate()
      cornpp$run_till_terminate()
      soywofost$run_till_terminate()
      soypp$run_till_terminate()

      # get output
      corn.output <- cornwofost$get_output()
      cornpp.output <- cornpp$get_output()
      soy.output <- soywofost$get_output()
      soypp.output <- soypp$get_output()

      # record output
      corn.sens.g[[deltaT.i]][[deltaRAIN.j]] <- lapply(corn.output, dplyr::bind_cols) |>
        dplyr::bind_rows() |>
        dplyr::left_join(
          lapply(cornpp.output, dplyr::bind_cols) |> dplyr::bind_rows(),
          by = "day", suffix = c(".WL", ".PP")
        )
      soy.sens.g[[deltaT.i]][[deltaRAIN.j]] <- lapply(soy.output, dplyr::bind_cols) |>
        dplyr::bind_rows() |>
        dplyr::left_join(
          lapply(soypp.output, dplyr::bind_cols) |> dplyr::bind_rows(),
          by = "day", suffix = c(".WL", ".PP")
        )
    }
  }
  corn.sens[[soils[g]]] <- dplyr::bind_rows(lapply(corn.sens.g, dplyr::bind_rows, .id = "deltaRAIN"), .id = "deltaT")
  soy.sens[[soils[g]]] <- dplyr::bind_rows(lapply(soy.sens.g, dplyr::bind_rows, .id = "deltaRAIN"), .id = "deltaT")
}

corn.sens.df <- dplyr::bind_rows(corn.sens, .id = "soil")
soy.sens.df <- dplyr::bind_rows(soy.sens, .id = "soil")

saveRDS(corn.sens.df,"results-raw/cornsensresults.rds")
saveRDS(soy.sens.df,"results-raw/soysensresults.rds")
