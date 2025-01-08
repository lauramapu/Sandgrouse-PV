rm(list=ls())
source('scripts/utils.R') # all used libraries and my custom functions

# pseudoabsences
# set of 10000 random points
# we already checked there was no difference between sets in the results

### load env_stack and presence points

env_stack <- list.files("spatial_data/env_asciis", pattern = "", full.names = TRUE)
env_stack <- raster::stack(env_stack)

# WE DO THIS FIRST WITH IBERICA

obs_iberica <- readRDS("objects/iberica_points.rds")

# load mask and generate pseudoabsences within mask (10000 points per set)
provinces <- mapSpain::esp_get_prov()
provinces <- provinces[!provinces$iso2.prov.name.es %in% c("Las Palmas", "Santa Cruz de Tenerife", "Baleares"), ]
provinces$campo <- 1
mask_spain <- dissolve(provinces, field='campo')

set.seed(1)
pseudoabsences <- as.data.frame(st_sample(mask_spain, size = 10000, type = "random"))

# convert to single spatial objects to later extract raster values
PA_to_points <- function(original_points) {
  coordinates_obj <- as.data.frame(st_coordinates(original_points$geometry)) 
  points_object <- st_as_sf(coordinates_obj, coords = c("X", "Y"), crs = 4326)
  return(points_object)
}
pseudoabsences_sf <- PA_to_points(pseudoabsences)

points <- rbind(obs_iberica, pseudoabsences_sf)

# we don't need the full set of env variables because landuses are categorical (1/0)
# so we select only the quantitative ones
stack_corr <- stack(env_stack$elev, env_stack$slope, env_stack$mean_temp,
                    env_stack$max_temp, env_stack$min_temp, env_stack$mean_temp_seas,
                    env_stack$prec, env_stack$prec_seas, env_stack$linear_infrastr)

# extract env values to both set of points
points_env <- as.data.frame(raster::extract(stack_corr, points))

# corr and vif test

# Calculate the correlation matrix
cor_matrix <- cor(points_env, use = "complete.obs")
# Transpose the data frame, needed for cluster
transposed_data <- t(points_env)
# Calculate Euclidean distances
dist_matrix <- dist(transposed_data)
# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "complete")  # You can use other linkage methods
# Plot the dendrogram
plot(hc, main = "Dendrograma", xlab = "Variables", sub = NULL)
# corrplot
corr_plot <- corrplot(cor_matrix, type = "full", method = "number")

# custom labels for visualization
custom_labels <- c("Elevation", "Slope", "Mean Temp.", "Max. Temp.", "Min. Temp.", 
                   "Temp. Seasonality", "Precipitation", "Precip. Seasonality", 
                   "Linear Infrastruct.")
colnames(cor_matrix) <- custom_labels
rownames(cor_matrix) <- custom_labels
corsave <- corrplot(cor_matrix, type = "full", method = "number", tl.srt = 45)

# high correlations in min and max temp, so erasing
points_env2 <- points_env[, -which(names(points_env) %in% c('max_temp', 'min_temp'))]

### VIF ANALYSIS
model <- lm(mean_temp ~ . -mean_temp, data = points_env2)
vifmodel <- as.data.frame(car::vif(model))
write.csv2(vifmodel, "results/vif_iberica.csv", row.names=T)
# no vif values over 5

# generate dataframes for modelling

# present
dataset <- as.data.frame(extract(env_stack, points))
dataset <- dataset %>% dplyr::select(-min_temp, -max_temp)
# generate presence/absence column as two-level factor
dataset$PA <- factor(c(rep(1, nrow(obs_iberica)),
                       rep(0, nrow(dataset) - nrow(obs_iberica))), levels = c(1, 0))
# bind coordinates
dataset <- cbind(dataset, st_coordinates(points))
# erase duplicates
duplicates <- duplicated(dataset[, -which(names(dataset) == c("PA","X","Y"))])
dataset <- subset(dataset, !duplicates)
# presences number
count(dataset[dataset$PA==1,])
# save
write.csv2(dataset, "data/iberica_PA_present.csv", row.names=F)

############################################

# ortega

env_stack <- list.files("spatial_data/env_asciis", pattern = "", full.names = TRUE)
env_stack <- raster::stack(env_stack)

obs_ortega <- readRDS("objects/ortega_points.rds")

points <- rbind(obs_ortega, pseudoabsences_sf)

points_env <- as.data.frame(extract(stack_corr, points))

cor_matrix <- cor(points_env, use = "complete.obs")
transposed_data <- t(points_env)
dist_matrix <- dist(transposed_data)
hc <- hclust(dist_matrix, method = "complete")
plot(hc, main = "Dendrograma", xlab = "Variables", sub = NULL)
corr_plot <- corrplot(cor_matrix, type = "full", method = "number")
points_env2 <- points_env[, -which(names(points_env) %in% c('max_temp', 'min_temp'))]
model <- lm(mean_temp ~ . -mean_temp, data = points_env2)
car::vif(model)
# same results

# present
dataset <- as.data.frame(extract(env_stack, points))
dataset <- dataset %>% dplyr::select(-min_temp, -max_temp)
dataset$PA <- factor(c(rep(1, nrow(obs_ortega)),
                       rep(0, nrow(dataset) - nrow(obs_ortega))), levels = c(1, 0))
# bind coordinates
dataset <- cbind(dataset, st_coordinates(points))
# erase duplicates
duplicates <- duplicated(dataset[, -which(names(dataset) == c("PA","X","Y"))])
dataset <- subset(dataset, !duplicates)
# presences number
count(dataset[dataset$PA==1,])
# save
write.csv2(dataset, "data/ortega_PA_present.csv", row.names=F)

################################
# preparation of env_df (same for both species)

# present
env_stack <- list.files("spatial_data/env_asciis", pattern = "", full.names = TRUE)
env_stack <- raster::stack(env_stack)
env_df <- na.omit(as.data.frame(env_stack, xy = T))
# add missing columns to match with the model's dataset 
env_df <- env_df %>%
  dplyr::select(-max_temp, -min_temp) %>%
  rename(X = x, Y = y) %>% # rename coords
  mutate(PA = 0) # create PA column
# create two levels for PA
env_df$PA <- factor(env_df$PA, levels = c("0", "1"))
write.csv2(env_df, "data/envdf_present.csv", row.names=F)

# ssp126
ssp126 <- list.files("spatial_data/env_asciis_fut", pattern = "126", full.names = TRUE)
ssp126 <- raster::stack(ssp126)
# erase present clim variables
env_stack_noclim <- env_stack[[setdiff(names(env_stack), 
                                       c("mean_temp", "mean_temp_seas", "min_temp",
                                         "max_temp", "prec", "prec_seas"))]]
# add future clim variables
env_stack_126 <- stack(env_stack_noclim, ssp126)
env_df <- na.omit(as.data.frame(env_stack_126, xy = T))
# change column names
columns <- names(env_df)
mod_names <- gsub("_126_4170", "", columns)
names(env_df) <- mod_names
# add missing columns to match with the model's dataset 
env_df <- env_df %>%
  dplyr::select(-max_temp, -min_temp) %>%
  rename(X = x, Y = y) %>% # rename coords
  mutate(PA = 0) # create PA column
# create two levels for PA
env_df$PA <- factor(env_df$PA, levels = c("0", "1"))
write.csv2(env_df, "data/envdf_ssp126.csv", row.names=F)

# ssp585
ssp585 <- list.files("spatial_data/env_asciis_fut", pattern = "585", full.names = TRUE)
ssp585 <- raster::stack(ssp585)
# add future clim variables
env_stack_585 <- stack(env_stack_noclim, ssp585)
env_df <- na.omit(as.data.frame(env_stack_585, xy = T))
# change column names
columns <- names(env_df)
mod_names <- gsub("_585_4170", "", columns)
names(env_df) <- mod_names
# add missing columns to match with the model's dataset 
env_df <- env_df %>%
  dplyr::select(-max_temp, -min_temp) %>%
  rename(X = x, Y = y) %>% # rename coords
  mutate(PA = 0) # create PA column
# create two levels for PA
env_df$PA <- factor(env_df$PA, levels = c("0", "1"))
write.csv2(env_df, "data/envdf_ssp585.csv", row.names=F)