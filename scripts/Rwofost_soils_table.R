# Rwofost soils table
soils <- Rwofost::wofost_soil() |> strsplit("\n|,\\s") |> unlist()

Rwofostsoil <- list()
for(i in seq_along(soils)){
  soil.i <- soils[i]
  Rws.i <- Rwofost::wofost_soil(soil.i)
  Rwofostsoil[[soil.i]] <- data.frame(
    SMFCF = Rws.i$SMFCF,
    SM0 = Rws.i$SM0,
    SMW = Rws.i$SMW,
    CRAIRC = Rws.i$CRAIRC,
    SOPE = Rws.i$SOPE,
    KSUB = Rws.i$KSUB,
    RDMSOL = Rws.i$RDMSOL,
    WAV = Rws.i$WAV
  )
}

Rwofostsoil.df <- dplyr::bind_rows(Rwofostsoil, .id = "soilid")
write.csv(Rwofostsoil.df, "data/Rwofostsoil.csv", row.names = FALSE)
