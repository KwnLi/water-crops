# must be run from home directory /home/kevin on Sarah's server!

# remotes::install_github("ClimateEcology/EnvToolbox", lib="/home/kevin/Rlib")

scp -P 1111 "data/pa-ws.csv" kevin@50.173.75.66:/home/kevin/
ssh -YC kevin@50.173.75.66 -p 1111

library(EnvToolbox, lib.loc="/home/kevin/Rlib")
library(sf)
# library(terra)

# hu <- sf::st_read("/home/kevin/selected_hu12_v2.gpkg")
hu <- st_read("/home/kevin/selected_hu12_v2_candidates.gpkg")

pa.ws <- read.csv("pa-ws.csv") |> dplyr::mutate(huc12 = as.character(paste0("0",huc12)))

hupa <- hu |> dplyr::filter(huc12 %in% pa.ws$huc12)

# parent directory
parentdir <- "/home/kevin/gisdata/hu_437/wofost/"

dir.create(parentdir, recursive = TRUE)

# extract soil data
dir.create(paste0(parentdir, "soil"))

for(i in 1:nrow(hupa)){
  hu.i <- hupa[i,]
  hu.name <- hu.i$huc12

  # root zone variables
  hu.ssurgo.hz <- querySoils(dat=hu.i,
                        return.grids = TRUE,
                        index = "gssurgo",
                        ssurgo.table = "chorizon",
                        soil.layers = c('wthirdbar_r','wsatiated_r','wfifteenbar_r','ksat_r'), # wofost params
                        horizon.depth = 'rootzone',     # depth-weighted mean from the rootzone variable in valu1
                        return.SSURGO.table = TRUE,
                        verbose=FALSE)

  terra::writeRaster(hu.ssurgo.hz[[1]],
                     paste0(parentdir, "soil/mukey_", hu.name,".tif"),
                     overwrite = TRUE)

  write.csv(hu.ssurgo.hz[[2]],
            paste0(parentdir, "soil/ssurgohz_",hu.name,".csv"), row.names = FALSE)

  # subsoil variable (KSUB)
  hu.ssurgo.sub<- querySoils(dat=hu.i,
                             return.grids = TRUE,
                             index = "gssurgo",
                             ssurgo.table = "chorizon",
                             soil.layers = c('ksat_r'), # wofost params
                             horizon.depth = 'subsoil',     # depth-weighted mean from the rootzone to max depth
                             return.SSURGO.table = TRUE,
                             verbose=FALSE)

  write.csv(hu.ssurgo.sub[[2]],
            paste0(parentdir, "soil/ssurgosub_",hu.name,".csv"), row.names = FALSE)

  hu.ssurgo.valu1 <- querySoils(dat=hu.i,
                             return.grids = TRUE,
                             index = "gssurgo",
                             ssurgo.table = "valu1",
                             soil.layers = c('rootznemc'), # wofost params
                             return.SSURGO.table = TRUE,
                             verbose=FALSE)

  write.csv(hu.ssurgo.valu1[[2]], paste0(parentdir, "soil/ssurgovalu1_",hu.name,".csv"), row.names = FALSE)
}


# extract climate data

dir.create(paste0(parentdir, "weather"))

date.range <- readRDS(file = "/home/kevin/ClimateDateRangeTable.rds")

for(i in 1:nrow(hupa)){
  hu.i <- hupa[i,]
  hu.name <- hu.i$huc12

  hu.gridmet <- getGrid_climate(vct=hu.i,
                               from.when = "2020-01-01", to.when = "2023-12-31",
                               index = "gridMET", period = "daily",
                               by=by, verbose=FALSE, showPB = TRUE)

  saveRDS(hu.gridmet, paste0(parentdir, "weather/gridmet.rds"))

  terra::writeRaster(hu.gridmet[[1]], paste0(parentdir, "weather/weatherID_",
                                            hu.name,".tif"), overwrite = TRUE)

  write.csv(hu.gridmet[[2]],
            paste0(parentdir, "weather/gridmet_",hu.name,".csv"), row.names = FALSE)
}
