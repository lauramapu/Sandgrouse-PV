
#####################################
####################################################
######################################################################
# NO RESAMPLING, ALL IN 1X1KM ####################################################
######################################################################
####################################################
#####################################

rm(list=ls())
source('scripts/utils.R') # all used libraries and my custom functions

# spanish IP borders
spain <- getpeninsularspain()

# construct database with all our data

full <- as.data.frame(raster('results/iberica_preds.asc'), xy=T)[,1:2]

filepaths <- grep('bin', list.files('results', pattern='preds.asc', full.names=T), value=T, invert=T)

for (filepath in filepaths) {
  r <- raster(filepath)
  x <- as.data.frame(r, xy=T)[,3]
  full <- cbind(full, x)
}

colnames(full)[3:8] <- c('ib_pres', 'ib_sust', 'ib_foss', 'or_pres', 'or_sust', 'or_foss')

# load energy infrastructures 

solar.cont <- raster('spatial_data/solar_energy/solar_continuous.asc') %>%
  projectRaster(raster('results/iberica_preds.asc'))
solar.bin <- raster('spatial_data/solar_energy/solar_binary.asc') %>%
  projectRaster(raster('results/iberica_preds.asc'), method='ngb')

# add solar to df
full$solar.cont <- solar.cont[]
full$solar.bin <- solar.bin[]

# labels for solar
# split zeros and non-zeros
values_non_zero <- full$solar.cont[full$solar.cont > 0]
# calculate quantiles of non-zeros
quantiles <- quantile(values_non_zero, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm=T)
# define custom breaks
breaks <- c(-Inf, 0, quantiles[2], quantiles[3], quantiles[4], quantiles[5])
# cut by custom breaks
full$labels_solar <- cut(full$solar.cont,
                         breaks = breaks,
                         labels = c('1','2','3','4','5'),
                         include.lowest = TRUE)
# check
table(full$labels_solar)

# load solar potentiality
pvout <- raster('spatial_data/solar_energy/SolarGIS/YearlyMonthlyTotals/PVOUT.tif') %>%
  projectRaster(raster('results/iberica_preds.asc')) %>%
  crop(extent(spain)) %>%
  mask(spain)
full$pvout <- pvout[]

quintiles <- quantile(full$pvout, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm=T)
full$labels_pvout <- cut(full$pvout,
                         breaks = quintiles,
                         labels = c('1','2','3','4','5'),
                         include.lowest = TRUE)

full$labels_pvout <- as.factor(full$labels_pvout)

# set fixed breaks and labels
breakz <- c(0,0.2,0.4,0.6,0.8,1)
labelz <- c(1,2,3,4,5)

# generate labels for every distribution
for (col in colnames(full)[3:8]) {
  new_col_name <- paste0('labs_', col)
  full[[new_col_name]] <- cut(full[[col]],
                              breaks = breakz,
                              labels = labelz,
                              include.lowest = TRUE)
}

# convert to factor and assign 5 factors (because some are missing in future scenarios)
labs_cols <- grep('^labs_', names(full), value = TRUE)
for (col in labs_cols) {
  full[[col]] <- factor(full[[col]], levels = c('1', '2', '3', '4', '5'))
}

# first we're doing the overlap of the initial observations with binary solar occupation

iberica_sf <- readRDS('objects/iberica_points.rds')
ortega_sf <- readRDS('objects/ortega_points.rds')

iberica <- raster('results/iberica_preds.asc')
# assign unique value to each pixel
iberica[] <- 1:ncell(iberica)
# raster::extract values to points
ncell <- raster::extract(iberica, iberica_sf)
# assign zero to iberica raster when ncell matches
iberica[] <- ifelse(values(iberica) %in% ncell, 1, 0)

ortega <- raster('results/ortega_preds.asc')
# assign unique value to each pixel
ortega[] <- 1:ncell(ortega)
# raster::extract values to points
ncell <- raster::extract(ortega, ortega_sf)
# assign zero to ortega raster when ncell matches
ortega[] <- ifelse(values(ortega) %in% ncell, 1, 0)

# construct database with presence of full species and a new column of conflict
full$iberica <- iberica[]
full$ortega <- ortega[]

full$ib_conflicts <- ifelse(is.na(full$solar.bin), 
                            NA, 
                            ifelse(full$iberica == 1 & full$solar.bin == 1, 3, 
                                   ifelse(full$iberica == 1 & full$solar.bin == 0, 2, 
                                          ifelse(full$iberica == 0 & full$solar.bin == 1, 1, 
                                                 ifelse(full$iberica == 0 & full$solar.bin == 0, 0, NA)))))
full$or_conflicts <- ifelse(is.na(full$solar.bin), 
                            NA, 
                            ifelse(full$ortega == 1 & full$solar.bin == 1, 3, 
                                   ifelse(full$ortega == 1 & full$solar.bin == 0, 2, 
                                          ifelse(full$ortega == 0 & full$solar.bin == 1, 1, 
                                                 ifelse(full$ortega == 0 & full$solar.bin == 0, 0, NA)))))

p1df <- full[,c('x','y','ib_conflicts')]
p1df <- p1df[!is.na(p1df$ib_conflicts), ]

p1 <- ggplot(p1df, aes(x = x, y = y, fill = as.factor(ib_conflicts))) +
  geom_raster() +
  scale_fill_manual(
    values = c('0' = 'transparent', '1' = 'darkgoldenrod2', '2' = '#4d7ce0', '3' = 'black'),
    na.value = 'transparent',
    name = ' ',
    labels = c('0' = '', '1' = 'PV', '2' = 'Obs.', '3' = 'Conflict')
  ) +
  theme(
    legend.box.background = element_rect(fill = 'transparent', color = NA),
    legend.background = element_rect(fill = 'transparent', color = NA),
    panel.background = element_rect(fill = 'transparent', color = NA),
    plot.background = element_rect(fill = 'transparent', color = NA)
  ) +
  labs(x = 'Longitude', y = 'Latitude',
       title = expression('a)' ~ italic('P. alchata') ~ 'Observations vs PV Infrastructures')) +
  theme_bw() +
  geom_sf(data = spain, fill = 'transparent', color = 'black', inherit.aes = FALSE)
p1 <- get_leyend(p1)
p1

table(full$ib_conflicts)

p2df <- full[,c('x','y','or_conflicts')]
p2df <- p2df[!is.na(p2df$or_conflicts), ]

p2 <- ggplot(p2df, aes(x = x, y = y, fill = as.factor(or_conflicts))) +
  geom_raster() +
  scale_fill_manual(
    values = c('0' = 'transparent', '1' = 'darkgoldenrod2', '2' = '#4d7ce0', '3' = 'black'),
    na.value = 'transparent',
    name = ' ',
    labels = c('0' = '', '1' = 'PV', '2' = 'Obs.', '3' = 'Conflict')
  ) +
  labs(x = 'Longitude', y = 'Latitude',
       title = expression('b)' ~ italic('P. orientalis') ~ 'Observations vs PV Infrastructures')) + 
  theme_bw() +
  geom_sf(data = spain, fill = 'transparent', color = 'black', inherit.aes = FALSE)
p2 <- get_leyend(p2)
p2

table(full$or_conflicts)

# maps for each distribution and continuous solar occupation
# same palette for all bivariated (yellow-blue-black)
mypal <- custom_bipal(min.xy = 'white',
                      max.y = '#e29e11',
                      max.x = '#4d7ce0',
                      max.xy = 'black',
                      dim=5)

biclass <- bi_class(full, x = labs_ib_pres, y = labels_solar, dim=5)
unique(biclass$bi_class)

p3 <- bimap(biclass,
            title = expression('c)' ~ italic('P. alchata') ~ 'Present Suitability vs PV Infrastructures'),
            xlab = 'Suitability',
            ylab = 'PV Infrastr.',
            mypal = mypal,
            dim=5)
p3

# save number of cells in each square of the bivariated legend in the same position as the legend
table(full$labels_solar, full$labs_ib_pres, useNA = 'no') %>%
  apply(2, rev) %>%
  write.csv('results/ib_pres.csv', row.names=T)

# # save as image
# png('results/test.png', height=150, width=300)
# p<-tableGrob(matriz)
# grid.arrange(p)
# dev.off()

biclass <- bi_class(full, x = labs_or_pres, y = labels_solar, dim=5)
unique(biclass$bi_class)

p4 <- bimap(biclass,
            title = expression('d)' ~ italic('P. orientalis') ~ 'Present Suitability vs PV Infrastructures'),
            xlab = 'Suitability',
            ylab = 'PV Infrastr.',
            mypal = mypal,
            dim=5)
p4

table(full$labels_solar, full$labs_or_pres, useNA = 'no') %>%
  apply(2, rev) %>%
  write.csv('results/or_pres.csv', row.names=T)

# load solar potentiality layer to overlap with future scenarios

biclass <- bi_class(full, x = labs_ib_sust, y = labels_pvout, dim=5)
unique(biclass$bi_class)

p5 <- bimap(biclass,
            title = expression('e)' ~ italic('P. alchata') ~ 'Sustainable Scenario Suitability vs PV Potentiality'),
            xlab = 'Suitability',
            ylab = 'PV Potentiality',
            mypal = mypal,
            dim=5)
p5

table(full$labels_pvout, full$labs_ib_sust, useNA = 'no') %>%
  apply(2, rev) %>%
  write.csv('results/ib_sust.csv', row.names=T)

biclass <- bi_class(full, x = labs_or_sust, y = labels_pvout, dim=5)
unique(biclass$bi_class)

p6 <- bimap(biclass,
            title = expression('f)' ~ italic('P. orientalis') ~ 'Sustainable Scenario Suitability vs PV Potentiality'),
            xlab = 'Suitability',
            ylab = 'PV Potentiality',
            mypal = mypal,
            dim=5)
p6

table(full$labels_pvout, full$labs_or_sust, useNA = 'no') %>%
  apply(2, rev) %>%
  write.csv('results/or_sust.csv', row.names=T)

biclass <- bi_class(full, x = labs_ib_foss, y = labels_pvout, dim=5)
unique(biclass$bi_class)

p7 <- bimap(biclass,
            title = expression('g)' ~ italic('P. alchata') ~ 'Fossil-fueled Scenario Suitability vs PV Potentiality'),
            xlab = 'Suitability',
            ylab = 'PV Potentiality',
            mypal = mypal,
            dim=5)
p7

table(full$labels_pvout, full$labs_ib_foss, useNA = 'no') %>%
  apply(2, rev) %>%
  write.csv('results/ib_foss.csv', row.names=T)

biclass <- bi_class(full, x = labs_or_foss, y = labels_pvout, dim=5)
unique(biclass$bi_class)

p8 <- bimap(biclass,
            title = expression('h)' ~ italic('P. orientalis') ~ 'Fossil-fueled Scenario Suitability vs PV Potentiality'),
            xlab = 'Suitability',
            ylab = 'PV Potentiality',
            mypal = mypal,
            dim=5)
p8

table(full$labels_pvout, full$labs_or_foss, useNA = 'no') %>%
  apply(2, rev) %>%
  write.csv('results/or_foss.csv', row.names=T)

conflicts <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol=2)
ggsave('results/conflicts_per_distr_v3.jpg', conflicts, height=20, width=14)

write.csv(full, 'data/all_data_1x1km.csv', row.names=F)

#####################################
####################################################
######################################################################
# RESAMPLING TO 5x5KM TO IMPROVE VISUALIZATION ######################################
######################################################################
####################################################
#####################################

rm(list=ls())
source('scripts/utils.R')

spain <- getpeninsularspain()

# construct database with all our data

full <- as.data.frame(raster('results/iberica_preds.asc'), xy=T)[,1:2]

filepaths <- grep('bin', list.files('results', pattern='preds.asc', full.names=T), value=T, invert=T)

for (filepath in filepaths) {
  r <- raster(filepath)
  x <- as.data.frame(r, xy=T)[,3]
  full <- cbind(full, x)
}

colnames(full)[3:8] <- c('ib_pres', 'ib_sust', 'ib_foss', 'or_pres', 'or_sust', 'or_foss')

# load energy infrastructures 

solar.cont <- raster('spatial_data/solar_energy/solar_continuous.asc') %>%
  projectRaster(raster('results/iberica_preds.asc'))
solar.bin <- raster('spatial_data/solar_energy/solar_binary.asc') %>%
  projectRaster(raster('results/iberica_preds.asc'), method='ngb')

# solar
full$solar.cont <- solar.cont[]
full$solar.bin <- solar.bin[]

# load solar potentiality
pvout <- raster('spatial_data/solar_energy/SolarGIS/YearlyMonthlyTotals/PVOUT.tif') %>%
  projectRaster(raster('results/iberica_preds.asc')) %>%
  crop(extent(spain)) %>%
  mask(spain)
full$pvout <- pvout[]

# species presence
iberica_sf <- readRDS('objects/iberica_points.rds')
ortega_sf <- readRDS('objects/ortega_points.rds')

# construct raster of same res and extent in which presence is 1 and else 0
# load base raster
iberica <- raster('results/iberica_preds.asc')
# assign unique value to each pixel
iberica[] <- 1:ncell(iberica)
# extract values to points
ncell <- raster::extract(iberica, iberica_sf)
# assign 1 to raster cell when ncell matches
iberica[] <- ifelse(values(iberica) %in% ncell, 1, 0)
# same for ortega
ortega <- raster('results/ortega_preds.asc')
ortega[] <- 1:ncell(ortega)
ncell <- raster::extract(ortega, ortega_sf)
ortega[] <- ifelse(values(ortega) %in% ncell, 1, 0)

# construct database with presence of full species and a new column of conflict
full$iberica <- iberica[]
full$ortega <- ortega[]

# resample

columns_max <- c('iberica', 'ortega', 'solar.cont', 'solar.bin')

columns_avg <- c('ib_pres', 'or_pres', 'ib_sust',
                 'or_sust', 'ib_foss', 'or_foss', 'pvout')

resampled_max <- lapply(columns_max, function(col) resample.max(full, col))
resampled_avg <- lapply(columns_avg, function(col) resample.avg(full, col))

resampled_max_df <- do.call(cbind, resampled_max)
resampled_avg_df <- do.call(cbind, resampled_avg)

resampled_df <- cbind(resampled_max_df, resampled_avg_df)
# erase duplicated xy
resampled_df <- resampled_df[, !duplicated(names(resampled_df))]

# labels for solar

# split zeros and non-zeros
values_non_zero <- resampled_df$solar.cont[resampled_df$solar.cont > 0]
#  calculate quantiles of non-zeros
quantiles <- quantile(values_non_zero, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm=T)
# define custom breaks
breaks <- c(-Inf, 0, quantiles[2], quantiles[3], quantiles[4], quantiles[5])
# cut by custom breaks
resampled_df$labels_solar <- cut(resampled_df$solar.cont,
                                 breaks = breaks,
                                 labels = c('1','2','3','4','5'),
                                 include.lowest = TRUE)
# check
table(resampled_df$labels_solar)

# labels for potentiality

quintiles <- quantile(resampled_df$pvout, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm=T)
resampled_df$labels_pvout <- cut(resampled_df$pvout,
                                 breaks = quintiles,
                                 labels = c('1','2','3','4','5'),
                                 include.lowest = TRUE)

resampled_df$labels_pvout <- as.factor(resampled_df$labels_pvout)
table(resampled_df$labels_pvout)

# set fixed breaks and labels
breakz <- seq(0, 1, by=0.2)
labelz <- c(1,2,3,4,5)

# generate labels for every distribution
for (col in columns_avg[1:6]) {
  new_col_name <- paste0('labs_', col)
  resampled_df[[new_col_name]] <- cut(resampled_df[[col]],
                                      breaks = breakz,
                                      labels = labelz,
                                      include.lowest = TRUE)
}

# convert to factor and assign 5 factors (because some are missing in future scenarios)
labs_cols <- grep('^labs_', names(resampled_df), value = TRUE)
for (col in labs_cols) {
  resampled_df[[col]] <- factor(resampled_df[[col]], levels = c('1', '2', '3', '4', '5'))
}

# overlaps of presence points with solar

resampled_df$ib_conflicts <- ifelse(is.na(resampled_df$solar.bin), 
                                    NA, 
                                    ifelse(resampled_df$iberica == 1 & resampled_df$solar.bin == 1, 3, 
                                           ifelse(resampled_df$iberica == 1 & resampled_df$solar.bin == 0, 2, 
                                                  ifelse(resampled_df$iberica == 0 & resampled_df$solar.bin == 1, 1, 
                                                         ifelse(resampled_df$iberica == 0 & resampled_df$solar.bin == 0, 0, NA)))))
resampled_df$or_conflicts <- ifelse(is.na(resampled_df$solar.bin), 
                                    NA, 
                                    ifelse(resampled_df$ortega == 1 & resampled_df$solar.bin == 1, 3, 
                                           ifelse(resampled_df$ortega == 1 & resampled_df$solar.bin == 0, 2, 
                                                  ifelse(resampled_df$ortega == 0 & resampled_df$solar.bin == 1, 1, 
                                                         ifelse(resampled_df$ortega == 0 & resampled_df$solar.bin == 0, 0, NA)))))

# same plots as before

p1df <- resampled_df[,c('x','y','ib_conflicts')]
p1df <- p1df[!is.na(p1df$ib_conflicts) & p1df$ib_conflicts != 0, ]

p1 <- ggplot(p1df, aes(x = x, y = y, fill = as.factor(ib_conflicts))) +
  geom_raster() +
  scale_fill_manual(
    values = c('0' = 'transparent', '1' = 'darkgoldenrod2', '2' = '#4d7ce0', '3' = 'black'),
    na.value = 'transparent',
    name = '',
    labels = c('1' = 'PV', '2' = 'Obs.', '3' = 'Conflict')
  ) +
  labs(x = 'Longitude',
       y = 'Latitude',
       title = expression(paste('a) ', italic('P. alchata '), 'Observations vs PV Infrastructures'))) +
  theme_bw() +
  geom_sf(data = spain, fill = 'transparent', color = 'black', inherit.aes = FALSE)
p1 <- get_leyend(p1)
p1

table(resampled_df$ib_conflicts)

p2df <- resampled_df[,c('x','y','or_conflicts')]
p2df <- p2df[!is.na(p2df$or_conflicts) & p2df$or_conflicts != 0, ]

p2 <- ggplot(p2df, aes(x = x, y = y, fill = as.factor(or_conflicts))) +
  geom_raster() +
  scale_fill_manual(
    values = c('0' = 'transparent', '1' = 'darkgoldenrod2', '2' = '#4d7ce0', '3' = 'black'),
    na.value = 'transparent',
    name = '',
    labels = c('0' = '', '1' = 'PV', '2' = 'Obs.', '3' = 'Conflict')
  ) +
  labs(x = 'Longitude',
       y = 'Latitude',
       title = expression(paste('b) ', italic('P. orientalis '), 'Observations vs PV Infrastructures'))) +
  theme_bw() +
  geom_sf(data = spain, fill = 'transparent', color = 'black', inherit.aes = FALSE)
p2 <- get_leyend(p2)
p2

table(resampled_df$or_conflicts)

# maps for each distribution and continuous solar occupation

biclass <- bi_class(resampled_df, x = labs_ib_pres, y = labels_solar, dim=5)
unique(biclass$bi_class)

mypal <- custom_bipal(min.xy = 'white',
                      max.y = '#e29e11',
                      max.x = '#4d7ce0',
                      max.xy = 'black',
                      dim=5)

p3 <- bimap(biclass,
            title = expression(paste('c) ', italic('P. alchata '), 'Present Suitability vs PV Infrastructures')),
            xlab = 'Suitability',
            ylab = 'PV Infrastr.',
            mypal = mypal,
            dim=5)
p3

table(resampled_df$labels_solar, resampled_df$labs_ib_pres, useNA = 'no') %>%
  apply(2, rev) %>%
  write.csv('results/ib_pres_5x5.csv', row.names=T)

biclass <- bi_class(resampled_df, x = labs_or_pres, y = labels_solar, dim=5)
unique(biclass$bi_class)

p4 <- bimap(biclass,
            title = expression(paste('d) ', italic('P. orientalis '), 'Present Suitability vs PV Infrastructures')),
            xlab = 'Suitability',
            ylab = 'PV Infrastr.',
            mypal = mypal,
            dim=5)
p4

table(resampled_df$labels_solar, resampled_df$labs_or_pres, useNA = 'no') %>%
  apply(2, rev) %>%
  write.csv('results/or_pres_5x5.csv', row.names=T)

# load solar potentiality layer to overlap with future scenarios

biclass <- bi_class(resampled_df, x = labs_ib_sust, y = labels_pvout, dim=5)
unique(biclass$bi_class)

p5 <- bimap(biclass,
            title = expression(paste('e) ', italic('P. alchata '), 'Sustainable Scenario Suitability vs PV Potentiality')),
            xlab = 'Suitability',
            ylab = 'PV Potentiality',
            mypal = mypal,
            dim=5)
p5

table(resampled_df$labels_pvout, resampled_df$labs_ib_sust, useNA = 'no') %>%
  apply(2, rev) %>%
  write.csv('results/ib_sust_5x5.csv', row.names=T)

biclass <- bi_class(resampled_df, x = labs_or_sust, y = labels_pvout, dim=5)
unique(biclass$bi_class)

p6 <- bimap(biclass,
            title = expression(paste('f) ', italic('P. orientalis '), 'Sustainable Scenario Suitability vs PV Potentiality')),
            xlab = 'Suitability',
            ylab = 'PV Potentiality',
            mypal = mypal,
            dim=5)
p6

table(resampled_df$labels_pvout, resampled_df$labs_or_sust, useNA = 'no') %>%
  apply(2, rev) %>%
  write.csv('results/or_sust_5x5.csv', row.names=T)

biclass <- bi_class(resampled_df, x = labs_ib_foss, y = labels_pvout, dim=5)
unique(biclass$bi_class)

p7 <- bimap(biclass,
            title = expression(paste('g) ', italic('P. alchata '), 'Fossil-fueled Scenario Suitability vs PV Potentiality')),
            xlab = 'Suitability',
            ylab = 'PV Potentiality',
            mypal = mypal,
            dim=5)
p7

table(resampled_df$labels_pvout, resampled_df$labs_ib_foss, useNA = 'no') %>%
  apply(2, rev) %>%
  write.csv('results/ib_foss_5x5.csv', row.names=T)

biclass <- bi_class(resampled_df, x = labs_or_foss, y = labels_pvout, dim=5)
unique(biclass$bi_class)

p8 <- bimap(biclass,
            title = expression(paste('h) ', italic('P. orientalis '),'Fossil-fueled Scenario Suitability vs PV Potentiality')),
            xlab = 'Suitability',
            ylab = 'PV Potentiality',
            mypal = mypal,
            dim=5)
p8

table(resampled_df$labels_pvout, resampled_df$labs_or_foss, useNA = 'no') %>%
  apply(2, rev) %>%
  write.csv('results/or_foss_5x5.csv', row.names=T)

conflicts <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol=2)
ggsave('results/conflicts_per_distr_5x5_v3.jpg', conflicts, height=20, width=14)

write.csv(resampled_df, 'data/all_data_5x5_v3.csv', row.names=F)

# overlap plots to show distribution of PV occupancy over the suitability distributions
# fixed 0.2 bins

rm(list=ls())
source('scripts/utils.R') # all used libraries and my custom functions

# spanish IP borders
spain <- getpeninsularspain()

full <- read.csv('data/all_data_1x1km.csv')

# select distri labels columns
labs_columns <- grep("^labs_", names(full), value = TRUE)

# change 1-5 to 0.2-1
full[labs_columns] <- lapply(full[labs_columns], function(col) {
  # Only modify non-NA values
  ifelse(!is.na(col), 
         ifelse(col == 1, 0.2,
                ifelse(col == 2, 0.4,
                       ifelse(col == 3, 0.6,
                              ifelse(col == 4, 0.8,
                                     ifelse(col == 5, 1, col))))), NA)
})

ib_pres <- na.omit(full) %>% group_by(labs_ib_pres) %>% summarize(n = n(), solar.bin = sum(solar.bin)) %>% rename(labs = labs_ib_pres)
ib_sust <- na.omit(full) %>% group_by(labs_ib_sust) %>% summarize(n = n(), solar.bin = sum(solar.bin)) %>% rename(labs = labs_ib_sust)
ib_foss <- na.omit(full) %>% group_by(labs_ib_foss) %>% summarize(n = n(), solar.bin = sum(solar.bin)) %>% rename(labs = labs_ib_foss)
or_pres <- na.omit(full) %>% group_by(labs_or_pres) %>% summarize(n = n(), solar.bin = sum(solar.bin)) %>% rename(labs = labs_or_pres)
or_sust <- na.omit(full) %>% group_by(labs_or_sust) %>% summarize(n = n(), solar.bin = sum(solar.bin)) %>% rename(labs = labs_or_sust)
or_foss <- na.omit(full) %>% group_by(labs_or_foss) %>% summarize(n = n(), solar.bin = sum(solar.bin)) %>% rename(labs = labs_or_foss)

dataframes <- list(ib_pres, ib_sust, ib_foss, or_pres, or_sust, or_foss)

# define maxs to plot
max_n <- dataframes %>%
  map_dbl(~ max(.x$n, na.rm = TRUE)) %>%
  max()

max_solar <- dataframes %>%
  map_dbl(~ max(.x$solar.bin, na.rm = TRUE)) %>%
  max()

# add missing rows to future scenarios 
# ib_sust <- rbind(ib_sust, c(1,0,0))
# ib_foss <- rbind(ib_foss, c(1,0,0))
or_sust <- rbind(or_sust, c(1,0,0))
or_foss <- rbind(or_foss, c(1,0,0))

g1 <- plot_counts(ib_pres, mytitle = expression(italic('P. alchata') ~ '- Present'))
g2 <- plot_counts(ib_sust, mytitle = expression(italic('P. alchata') ~ '- Sustainable Fut.'))
g3 <- plot_counts(ib_foss, mytitle = expression(italic('P. alchata') ~ '- Fossil-Fueled Fut.'))
g4 <- plot_counts(or_pres, mytitle = expression(italic('P. orientalis') ~ '- Present'))
g5 <- plot_counts(or_sust, mytitle = expression(italic('P. orientalis') ~ '- Sustainable Fut.'))
g6 <- plot_counts(or_foss, mytitle = expression(italic('P. orientalis') ~ '- Fossil-Fueled Fut.'))

g7 <- grid.arrange(g1,g2,g3,g4,g5,g6, ncol=6, nrow=1)
ggsave('results/cell_counts_fifths.jpg', g7, height = 5, width = 18, units='in')

# arrange data as unified dataframe

names <- c('ib_pres', 'ib_sust', 'ib_foss', 'or_pres', 'or_sust', 'or_foss')

alldata <- Map(function(df, name) {
  df %>% mutate(distri = name)
}, dataframes, names)

alldata <- bind_rows(alldata)
write.csv(alldata, 'results/overlaps_fifths.csv', row.names=F)
