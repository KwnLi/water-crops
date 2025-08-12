devtools::load_all()

reticulate::use_condaenv("py3_pcse")
pcse <- reticulate::import("pcse")

#### Run WOFOST using our data ####

# load soil file
soildir <- "data/soil/"
soilfiles <- list.files(soildir)

testsoil <- pcse$input$PCSEFileReader(paste0(soildir,soilfiles[1]))

# load weather file
weatherdir <- "data/weather/"
weatherfiles <- list.files(weatherdir)

testweather <- pcse$input$CSVWeatherDataProvider(here::here(paste0(weatherdir,weatherfiles[1])))

# site
testsite <- pcse$input$WOFOST72SiteDataProvider(WAV=10)

# crop
testcrop <- pcse$input$YAMLCropDataProvider()

# agromanagement
testagro <- pcse$input$YAMLAgroManagementReader("test.agro")

# parameters
testparam <- pcse$base$ParameterProvider(sitedata=testsite, soildata = testsoil, cropdata = testcrop)

testwofost <- pcse$models$Wofost72_WLP_CWB(testparam, testweather, testagro)

test.output <- testwofost$get_output()

test.outputdf <- lapply(test.output, dplyr::bind_cols) |> dplyr::bind_rows()

#### Run WOFOST using default parameters ####
