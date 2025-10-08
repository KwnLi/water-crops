library(reticulate)
library(sf)
library(terra)
# gitdir <- "/Users/kevinl/Documents/GitHub/"
# devtools::load_all(paste0(gitdir,"rPCSE"))
# remotes::install_github("ClimateEcology/rPCSE")
library(rPCSE)
library(Rwofost)
library(ggplot2)

reticulate::use_condaenv("py3_pcse")
pcse <- reticulate::import("pcse")

#### soil data processing ####
if(!dir.exists("data/Rwofost_soil")) dir.create("data/Rwofost_soil", recursive = TRUE)

# Rwofost soils
soils <- wofost_soil() |> strsplit("\n|,\\s") |> unlist()

for(k in seq_along(soils)){
  id.k <-soils[k]

  make_soil_pcse(outpath = "data/Rwofost_soil", fname = id.k,
                 Rwofost.dat=wofost_soil(id.k), Rwofost.id=id.k)
}

### review the available weather data
gridweather <- list.files("data/weather")

weathermns <- list()

for(j in seq_along(gridweather)){
  weather.j <- read.csv(paste0("data/weather/",gridweather[j]),skip = 8) |>
    dplyr::mutate(TMEAN = (TMIN+TMAX)/2) |>
    dplyr::summarize(Tmn = mean(TMEAN),
                     RAINmn = mean(RAIN),
                     .groups = "drop")

  weathermns[[gridweather[j]]] <- weather.j
}

weathermns <- weathermns |> dplyr::bind_rows(.id="gridid") |>
  dplyr::mutate(allmnT = mean(Tmn), allmnRAIN = mean(RAINmn)) |>
  dplyr::rowwise() |>
  dplyr::mutate(distCenter = sqrt((Tmn-allmnT)^2 + (RAINmn-allmnRAIN)^2)) |>
  dplyr::ungroup()

write.csv(weathermns, "data/weathergrid.csv",row.names = FALSE)

# middlemost grid
midgrid <- weathermns$gridid[which(weathermns$distCenter==min(weathermns$distCenter))]
midgrid.dat <- read.csv(paste0("data/weather/",midgrid),skip=8)

# bump up and down

Trange <- seq(-1,1, length.out=11)
Rrange <- seq(-0.2,0.2, length.out=11)
bumps <- list()
bumpsmns <- list()

for(i in seq_along(Trange)){
  deltT <- Trange[i]

  for(j in seq_along(Rrange)){
    deltR <- Rrange[j]

    bumpdat.ij <-
      midgrid.dat |> dplyr::mutate(TMIN = TMIN + deltT,
                                   TMAX = TMAX + deltT,
                                   RAIN = RAIN + deltR) |>
      dplyr::mutate(RAIN = ifelse(RAIN<0,0,RAIN))

    bumps[[as.character(deltT)]][[as.character(deltR)]] <- bumpdat.ij

    bumpsmns[[as.character(deltT)]][[as.character(deltR)]] <- bumpdat.ij |>
      dplyr::mutate(TMEAN = (TMAX+TMIN)/2) |>
      dplyr::summarize(Tmn = mean(TMEAN), RAINmn = mean(RAIN), .groups = "drop")
  }
  bumpsmns[[as.character(deltT)]] <- dplyr::bind_rows(bumpsmns[[as.character(deltT)]], .id="deltR")
}
bumpsmns.df <- dplyr::bind_rows(bumpsmns, .id="deltT")

write.csv(bumpsmns.df,"data/weathergridmns.csv", row.names = FALSE)

saveRDS(bumps, "data/weathergrid.rds")

# from quarter to double rain
Trange <- seq(-1,1, length.out=11)
Rrange <- seq(0.25,2, by = 0.25)
bumps2 <- list()
bumpsmns2 <- list()

for(i in seq_along(Trange)){
  deltT <- Trange[i]

  for(j in seq_along(Rrange)){
    deltR <- Rrange[j]

    bumpdat.ij <-
      midgrid.dat |> dplyr::mutate(TMIN = TMIN + deltT,
                                   TMAX = TMAX + deltT,
                                   RAIN = RAIN * deltR)

    bumps2[[as.character(deltT)]][[as.character(deltR)]] <- bumpdat.ij

    bumpsmns2[[as.character(deltT)]][[as.character(deltR)]] <- bumpdat.ij |>
      dplyr::mutate(TMEAN = (TMAX+TMIN)/2) |>
      dplyr::summarize(Tmn = mean(TMEAN), RAINmn = mean(RAIN), .groups = "drop")
  }
  bumpsmns2[[as.character(deltT)]] <- dplyr::bind_rows(bumpsmns2[[as.character(deltT)]], .id="deltR")
}
bumpsmns2.df <- dplyr::bind_rows(bumpsmns2, .id="deltT")
write.csv(bumpsmns2.df,"data/weathergridmns2.csv", row.names = FALSE)

saveRDS(bumps2, "data/weathergrid2.rds")

