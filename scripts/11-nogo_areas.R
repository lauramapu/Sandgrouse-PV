
rm(list=ls())
source('scripts/utils.R') # all used libraries and my custom functions

# spanish IP borders
spain <- getpeninsularspain()

# load ebird data (points) and convert to raster the same way we did for sampling observations
pts_ebird_points <- st_read('data/ebird_alchata_4326.shp')
ort_ebird_points <- st_read('data/ebird_orientalis_4326.shp')

pts_ebird <- raster('results/iberica_preds.asc')
# assign unique value to each pixel
pts_ebird[] <- 1:ncell(pts_ebird)
# raster::extract values to points
ncell <- raster::extract(pts_ebird, pts_ebird_points)
# assign zero to pts_ebird raster when ncell matches
pts_ebird[] <- ifelse(values(pts_ebird) %in% ncell, 1, 0)

ort_ebird <- raster('results/iberica_preds.asc')
# assign unique value to each pixel
ort_ebird[] <- 1:ncell(ort_ebird)
# raster::extract values to points
ncell <- raster::extract(ort_ebird, ort_ebird_points)
# assign zero to ort_ebird raster when ncell matches
ort_ebird[] <- ifelse(values(ort_ebird) %in% ncell, 1, 0)

writeRaster(pts_ebird, 'data/alchata_ebird.asc', format='ascii', overwrite=T)
writeRaster(ort_ebird, 'data/orientalis_ebird.asc', format='ascii', overwrite=T)

data <- as.data.frame(stack('results/iberica_r_presence.asc',
                            'results/ortega_r_presence.asc',
                            'data/alchata_ebird.asc',
                            'data/orientalis_ebird.asc',
                            'results/iberica_preds.asc',
                            'results/ortega_preds.asc',
                            'results/iberica_ssp126_preds.asc',
                            'results/ortega_ssp126_preds.asc',
                            'results/iberica_ssp585_preds.asc',
                            'results/ortega_ssp585_preds.asc',
                            'spatial_data/solar_energy/solar_binary.asc'), xy=T) %>%
  
  # binarize preds, for future ortega we use 0.6 because very low values
  mutate(
    bin_iberica = ifelse(iberica_preds >= 0.8, 1, 0),
    bin_ortega = ifelse(ortega_preds >= 0.8, 1, 0),
    bin_iberica_ssp126 = ifelse(iberica_ssp126_preds >= 0.6, 1, 0),
    bin_ortega_ssp126 = ifelse(ortega_ssp126_preds >= 0.6, 1, 0),
    bin_iberica_ssp585 = ifelse(iberica_ssp585_preds >= 0.8, 1, 0),
    bin_ortega_ssp585 = ifelse(ortega_ssp585_preds >= 0.6, 1, 0)
  )

# nogo cells

# first sum up presences

data <- data %>%
  mutate(presence = rowSums(dplyr::select(., iberica_r_presence, ortega_r_presence, alchata_ebird, orientalis_ebird)))

#data$nogo[data$nogo==0] <- NA

presence_r <- rasterFromXYZ(data[,c('x','y','presence')], res = res(pts_ebird), crs=crs(pts_ebird)); mapview(presence_r)
presence_p <- rasterToPolygons(presence_r, fun = function(x){x>0})

# now presence + present distris

data <- data %>%
  mutate(presence_present = rowSums(dplyr::select(., iberica_r_presence, ortega_r_presence, alchata_ebird, orientalis_ebird,
                                                  bin_iberica, bin_ortega)))

presence_present_r <- rasterFromXYZ(data[,c('x','y','presence_present')], res = res(pts_ebird), crs=crs(pts_ebird)); mapview(presence_present_r)
presence_present_p <- rasterToPolygons(presence_present_r, fun = function(x){x>0})

# now presence + present distris filtered by presence

data <- data %>%
  mutate(
    presence_present_filtered = ifelse(presence>0,
      rowSums(dplyr::select(., iberica_r_presence, ortega_r_presence, alchata_ebird, orientalis_ebird,
                            bin_iberica, bin_ortega), na.rm=T),  # sum the binary distributions
      0  # and else assign zero
    )
  )

presence_present_filtered_r <- rasterFromXYZ(data[,c('x','y','presence_present_filtered')], res = res(pts_ebird), crs=crs(pts_ebird)); mapview(presence_present_filtered_r)
presence_present_filtered_p <- rasterToPolygons(presence_present_filtered_r, fun = function(x){x>0})

# now all the presence and distris together

data <- data %>%
  mutate(
    alldistris = ifelse(presence>0,
                        rowSums(dplyr::select(., iberica_r_presence, ortega_r_presence, alchata_ebird, orientalis_ebird,
                                              bin_iberica, bin_ortega, bin_iberica_ssp126, bin_ortega_ssp126,
                                              bin_iberica_ssp585, bin_ortega_ssp585), na.rm=T),  # sum the binary distributions
                        0  # and else assign zero
    )
  )

alldistris_r <- rasterFromXYZ(data[,c('x','y','alldistris')], res = res(pts_ebird), crs=crs(pts_ebird)); mapview(alldistris_r)
alldistris_p <- rasterToPolygons(alldistris_r, fun = function(x){x>0})

# all distris not filtered by presence

data <- data %>%
  mutate(
    alldistris_nofilt = rowSums(dplyr::select(., iberica_r_presence, ortega_r_presence, alchata_ebird, orientalis_ebird,
                                              bin_iberica, bin_ortega, bin_iberica_ssp126, bin_ortega_ssp126,
                                              bin_iberica_ssp585, bin_ortega_ssp585), na.rm=T)  # sum the binary distributions
  )

alldistris_nofilt_r <- rasterFromXYZ(data[,c('x','y','alldistris_nofilt')], res = res(pts_ebird), crs=crs(pts_ebird)); mapview(alldistris_nofilt_r)
alldistris_nofilt_p <- rasterToPolygons(alldistris_nofilt_r, fun = function(x){x>0})

# plot everything
mapview(presence_p) + mapview(presence_present_filtered_p) + mapview(alldistris_p)

# 'plot whatever column you want' function

plot.any <- function(x) {
  r <- rasterFromXYZ(data[,c('x','y',x)], res = res(pts_ebird), crs=crs(pts_ebird))
  r <- r>0; mapview(r)
  # p <- rasterToPolygons(r, fun = function(x){x>0}); mapview(p)
}

plot.any('bin_iberica_ssp585')


## plot along with solar
pv <- rasterFromXYZ(data[,c('x','y','solar_binary')], res = res(pts_ebird), crs=crs(pts_ebird)) %>%
  rasterToPolygons(fun = function(x){x>0})

mapview(alldistris_nofilt_r) + mapview(presence_p)


mapview(alldistris_p) + mapview(alldistris_nofilt_r) 

writeRaster(alldistris_r, 'results/alldistris_r.asc', format='ascii', overwrite=T)
write_sf(alldistris_points, 'results/alldistris_points.shp')

############################################################

### plot nicely

# resampling to 5x5 to better vissualize

columns_max <- c('alldistris', 'alldistris_nofilt')
resampled_max <- lapply(columns_max, function(col) resample.max(data, col))
resampled_max_df <- do.call(cbind, resampled_max)
resampled_df <- resampled_max_df[, !duplicated(names(resampled_max_df))]

onlypresence <- resampled_df %>% filter(alldistris>0)
everything <- resampled_df %>% filter(alldistris_nofilt>0)

ccaa <- esp_get_ccaa(epsg = 4326, moveCAN = F) # with CCAA
ccaa <- ccaa[!ccaa$nuts2.name %in% c("Illes Balears", "Canarias"),]

plot1 <- ggplot(onlypresence, aes(x = x, y = y, fill = alldistris)) +
  geom_raster() +
  scale_fill_viridis(title('')) +
  labs(x = "Longitude", y = "Latitude", title = "Presence cells: sum of presences + suitabilities") +
  theme_bw() +
  geom_sf(data = ccaa, fill = "transparent", color = "black", inherit.aes = FALSE)
plot1 <- get_leyend(plot1, x = 0.70, y = 0.15)
plot1

plot2 <- ggplot(everything, aes(x = x, y = y, fill = alldistris_nofilt)) +
  geom_raster() +
  scale_fill_viridis(title('Sum')) +
  labs(x = "Longitude", y = "Latitude", title = "All cells: sum of presences + suitabilities") +
  theme_bw() +
  geom_sf(data = ccaa, fill = "transparent", color = "black", inherit.aes = FALSE)
plot2 <- get_leyend(plot2, x = 0.70, y = 0.15)
plot2

plot3 <- grid.arrange(plot1, plot2, ncol=2)
ggsave('results/nogo-areas.jpg', plot3, height=5, width=14)
