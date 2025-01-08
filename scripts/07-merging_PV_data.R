library(sf)
library(raster)
library(dplyr)
library(mapview)
library(MazamaSpatialUtils)
library(mapSpain)
library(units)
library(terra) # better for rasterizing

# some of the following processes are very heavy so it's highly recommended to run in supercluster

# spanish provinces borders
provinces <- mapSpain::esp_get_prov()
# excluding Canarias and Baleares
provinces <- provinces[!provinces$iso2.prov.name.es %in% c("Las Palmas", "Santa Cruz de Tenerife", "Baleares",
                                                           "Melilla", "Ceuta"), ]
# dissolve
provinces$campo <- 1
mask_spain <- provinces %>%
  dissolve(field='campo') %>%
  st_transform(crs = 3035)

gpkg_files <- list.files("spatial_data/siose_ar", pattern = "\\.gpkg$",
                         full.names = TRUE, recursive = TRUE)
                         # with recursive=T it searches in all subfolders

# extract the USOS layer
extract_solar <- function(file_path) {
  # read all layers
  layers <- st_layers(file_path)
  # find the one that matches _USOS in its name
  layer_name <- layers$name[grep("_T_USOS$", layers$name)]
  if (length(layer_name) == 0) {
    # if not found return a warning
    stop("'_T_USOS' layer not found")
  }
  # read layer
  layer <- st_read(file_path, layer = layer_name[1])
  # select only solar and wind polygons
  layer <- layer[layer$ID_USO_MAX %in% c(2442), ]
  return(layer)
}

# erase baleares [6], las palmas [30], tenerife [32], ceuta [45],
# melilla [46] and posesiones [53]
capas <- gpkg_files[-c(6,30,32,45,46,53)]

# better run on cluster because it takes a while!
capas2 <- lapply(capas, extract_layer)

# project all to etrs89 30n (25830)
reproject <- function(sf_layer) {
  reprojected_layer <- st_transform(sf_layer, crs = 3035)
  return(reprojected_layer)
}
capas3 <- lapply(capas2, reproject)

# join layers

# initialize object to store layers
sf_merged <- capas3[[1]]
# iterate through each layer
for (i in 2:length(capas3)) {
  if (ncol(capas3[[i]]) > 17) {
    capas3[[i]] <- capas3[[i]][,-1] # erase the id column
    sf_merged <- rbind(sf_merged, capas3[[i]]) # add layer
  } else {  # ncol = 17 ok
    sf_merged <- rbind(sf_merged, capas3[[i]])
  }
}
mapview(sf_merged)

write_sf(sf_merged, "spatial_data/solar_energy/siose_solar.shp")

#############################################
### merging of all layers ###
############################################

# siose
siose <- st_read("spatial_data/solar_energy/siose_solar.shp")
# andalucia
andalucia <- st_read("spatial_data/solar_energy/other_sources/andalucia/solar.gpkg")
# aragon
aragon <- st_read("spatial_data/solar_energy/other_sources/aragon/v_vall_aae.shp")
# cataluña
# catalunya <- st_read("spatial_data/other_sources/catalunya/solar_multipoly.gpkg")
# write_sf(catalunya, "spatial_data/other_sources/catalunya/solar_multipoly.shp")
catalunya <- st_read("spatial_data/solar_energy/other_sources/catalunya/solar_multipoly.shp")
catalunya2 <- st_read("spatial_data/solar_energy/other_sources/catalunya/solar_points.gpkg")
# navarra
navarra <- st_read("spatial_data/solar_energy/other_sources/navarra/solar_funcionando.shp")
# valencia
# valencia <- st_read("spatial_data/other_sources/valencia/solar.gpkg")
# write_sf(valencia, "spatial_data/other_sources/valencia/solar_multipoly.shp")
valencia <- st_read("spatial_data/solar_energy/other_sources/valencia/solar_multipoly.shp")
# mis punticos
punticos <- st_join(st_read("spatial_data/solar_energy/mis_punticos.shp"),
                    st_read("spatial_data/solar_energy/punticos_angela.shp"))
punticos <- st_join(punticos, st_read("spatial_data/solar_energy/punticos_elena.shp"))
punticos <- st_join(punticos, st_transform(st_read("spatial_data/solar_energy/punticos_ana.shp"),
                                           crs = st_crs(punticos)))
# SatlasPreTrain
satlas <- st_read('spatial_data/solar_energy/satlas_solar/2024-01_solar.shp') %>%
  st_transform(crs = 3035) %>%
  st_intersection(mask_spain)

# list all files
lista <- list(siose, andalucia, aragon, catalunya,
              catalunya2, navarra, valencia, punticos, satlas)

topolys <- function(x) {
  
  # reproject x
  x <- st_transform(x, crs = 3035)
  
  # if sf layer is points convert to polygons
  if (head(class(x$geom), 1) %in% c('sfc_MULTIPOINT', 'sfc_POINT')) {
    
    # extract coords and generate sf object
    coords <- as.data.frame(st_coordinates(x))
    puntos <- st_as_sf(coords, coords = c('X','Y'), crs = 3035)
    # generate 250m buffer around each point
    buff <- st_buffer(puntos, dist=250)
    return(buff)
  }
  
  # if POLYGONS generate super small buffer to fix invalid geometries
  else { 
    buff <- st_buffer(x, dist= 0.01) # 1cm
    return(buff)
  }
}

new_list <- lapply(lista, topolys)

# initialize object to store layers
allpolys <- new_list[[1]] %>% dplyr::select(geometry)
# iterate through each layer
for (i in 2:length(new_list)) {
  geometries <- new_list[[i]] %>% dplyr::select(geometry)
  allpolys <- rbind(allpolys, geometries)
}

plot(allpolys)

write_sf(allpolys, 'spatial_data/solar_energy/allsolar.shp')

toerase <- setdiff(ls(), # list all environment
                   c('capas', 'reproject', 'mask_spain')) # keep this
rm(list = toerase)

allsolar <- st_read('spatial_data/solar_energy/allsolar.shp') %>%
  st_make_valid()

# now we need to cut this polygons by urban layer (corine because less geometry issues)
# and by greenhouses uses (from SIOSE AR, SIGPAC uses)

# first we're extracting greenhouses from SIOSE

# extract the greenhouses polygons
extract_greenhouse <- function(file_path) {
  # read all layers
  layers <- st_layers(file_path)
  # find the one that matches POLIGONOS in its name
  layer_name <- layers$name[grep("POLIGONOS$", layers$name)]
  if (length(layer_name) == 0) {
    # if not found return a warning
    warning("'POLIGONOS' layer not found")
  }
  # read layer
  layer <- st_read(file_path, layer = layer_name[1])
  # select only greenhouse polygons
  layer <- layer[layer$USO_SIGPAC %in% c('IV'), ]
  return(layer)
}

capas2 <- lapply(capas, extract_greenhouse)
capas3 <- lapply(capas2, reproject)

# join layers
# initialize object to store layers
greenhouse <- capas3[[1]][,c('ID_POLYGON','geom')]
# iterate through each layer
for (i in 2:length(capas3)) {
  capas3[[i]] <- capas3[[i]][,c('ID_POLYGON','geom')]
  greenhouse <- rbind(greenhouse, capas3[[i]])
}

write_sf(greenhouse, "spatial_data/siose_ar/siose_IV.shp")
# greenhouse <- st_read("spatial_data/siose_ar/siose_IV.shp")

greenhouse$campo <- 1
greenhouse <- dissolve(greenhouse, field = 'campo')

rm(capas2, capas3)

# now we can cut solar by greenhouses

cut1 <- st_difference(allsolar, greenhouse)

# now we're extracting the urban layer from corine (2018)
# this step is crazy computationally so don't run in local if possible

clc_urban <- st_read('spatial_data/corine_2018/clc2018.gpkg') %>%
  st_transform(crs = st_crs(mask_spain)) %>%
  st_intersection(mask_spain) %>%
  filter(Code_18 %in% # urban uses
           c('111', '112', '122', '123', '124', '142')) 
# use 121 is industrial and includes some solar farms so cant erase
# 133 is construction sites and covers some previously unfinished farms

write_sf(clc_urban, 'spatial_data/corine_2018/urban_uses.shp')

clc_urban$campo <- 1
clc_urban <- dissolve(clc_urban, field = 'campo')

finalsolar <- st_difference(cut1, clc_urban)

write_sf(finalsolar, 'spatial_data/solar_energy/allsolar_cut.shp')

# calculate areas per pixel cell

# method one

# read raster base and solar final polys in terra format
base <- rast("results/iberica_preds.asc") %>%
  project("EPSG:3035")
base[] <- 1:ncell(base)
solar <- vect('spatial_data/solar_energy/allsolar_cut.shp')

# rasterize polygons
cont <- terra::rasterize(solar, base, background=0, touches=F, cover=T)
# with cover=T it estimates the fraction of cell that is covered by polygons
mask <- vect(mask_spain)
cropped <- crop(cont, mask)
solar.cont.r <- mask(cont, mask)

mapview(solar.cont.r, maxpixels=2025093) + mapview(solar) # super done

# to raster format, resample to original and save
raster(solar.cont.r) %>%
  projectRaster(raster("results/iberica_preds.asc")) %>%
  raster::writeRaster('spatial_data/solar_energy/solar_continuous.asc', format='ascii', overwrite=T)

# binarize solar presence
solar.bin.r <- solar.cont.r
solar.bin.r[] <- ifelse(values(solar.cont.r) > 0, 1, 0)
raster::writeRaster(solar.bin.r, 'spatial_data/solar_energy/solar_binary.asc', format='ascii', overwrite=T)

# method two

base <- rast("results/iberica_preds.asc") %>%
  project("EPSG:3035")
base[] <- 1:ncell(base)
solar <- vect('spatial_data/solar_energy/allsolar_cut_diss.shp') # diss and area from QGIS

cont <- terra::rasterize(solar, base, background=0, touches=F, field='area', fun=sum)
mask <- vect(mask_spain)
cropped <- crop(cont, mask)
solar.cont.r <- mask(cont, mask)

mapview(solar.cont.r, maxpixels=2025093) + mapview(solar) # super done

# to raster format, resample to original and save
raster(solar.cont.r) %>%
  projectRaster(raster("results/iberica_preds.asc")) %>%
  raster::writeRaster('spatial_data/solar_energy/solar_continuous_funsum.asc', format='ascii', overwrite=T)

# method one overestimates, method two underestimates

# method three: package exactextractr
library(exactextractr)

base <- raster("results/iberica_preds.asc") %>%
  projectRaster(crs=3035)
base[] <- 1:ncell(base)
solar <- st_read('spatial_data/solar_energy/allsolar_cut_diss.shp') # diss and area from QGIS

cont <- coverage_fraction(base, solar)
cont1 <- stack(cont)
cont2 <- calc(cont1, sum)
rm(cont, cont1)

cropped <- crop(cont2, extent(mask_spain))
solar.cont.r <- mask(cropped, mask_spain)

solar.cont.r %>%
  projectRaster(raster("results/iberica_preds.asc")) %>%
  raster::writeRaster('spatial_data/solar_energy/solar_continuous_exactextractr.asc', format='ascii', overwrite=T)
