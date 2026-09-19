suppressPackageStartupMessages({library(sf); library(dplyr); library(readr)})
sf_use_s2(TRUE)
raw <- "data/raw"
adm1 <- st_read(file.path(raw,"geoBoundaries-IDN-ADM1.geojson"), quiet=TRUE)
cat("ADM1 names:\n"); print(sort(adm1$shapeName))
