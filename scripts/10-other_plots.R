# other plots for the paper

rm(list=ls())
source('scripts/utils.R') # all used libraries and my custom functions

# observations used in modeling for each species

provinces <- mapSpain::esp_get_prov()[!mapSpain::esp_get_prov()$iso2.prov.name.es %in%
                                        c("Las Palmas", "Santa Cruz de Tenerife", "Baleares", "Melilla", "Ceuta"), ] %>%
  st_transform(crs=4326)

iberica <- read.csv2('data/iberica_PA_present.csv') %>%
  dplyr::filter(PA == 1) %>%
  st_as_sf(coords=c('X','Y'), crs=4326)

ortega <- read.csv2('data/ortega_PA_present.csv') %>%
  dplyr::filter(PA == 1) %>%
  st_as_sf(coords=c('X','Y'), crs=4326) 

p1 <- ggplot() +
  geom_sf(data = provinces, fill = NA, color = "black") +
  geom_sf(data = iberica, color = "cyan4") +
  theme_bw() +
  labs(title = expression(paste('a) ', italic('P. alchata'))),
    caption = 'n = 375')

p2 <- ggplot() +
  geom_sf(data = provinces, fill = NA, color = "black") +
  geom_sf(data = ortega, color = "darkred") +
  theme_bw() +
  labs(title = expression(paste('b) ', italic('P. orientalis'))),
       caption = 'n = 324')

p3 <- grid.arrange(p1, p2, ncol = 2, nrow = 1)
ggsave("data/obs_map.jpg", p3,
       width = 8, height = 3, units = 'in')

# summary of preds

topcells <- read.csv('data/all_data_1x1km.csv')

stat_box_data <- function(y) {
  return(data.frame(
    y = max(y),
    label = paste('Max:', round(max(y), 2), 
                  '\nMin:', round(min(y), 2), 
                  '\nAvg:', round(mean(y), 2),
                  '\nMedian:', round(median(y), 2))
  ))
}

ib.long <- topcells[,c('ib_pres','ib_sust','ib_foss')] %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
ordered_levels <- c('ib_pres', 'ib_sust', 'ib_foss')
ib.long$variable <- factor(ib.long$variable, levels = ordered_levels)

p1 <- ggplot(ib.long, aes(x = variable, y = value)) +
  geom_boxplot() +
  stat_summary(fun.data = stat_box_data, geom = "text", 
               hjust = 1, vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(title = expression(paste(italic("P. alchata"))),
       x = "",
       y = "Prediction") +
  scale_x_discrete(labels = c('Present','Sustainable','Fossil-fueled')) + 
  coord_flip()
p1

or.long <- topcells[,c('or_pres','or_sust','or_foss')] %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
ordered_levels <- c('or_pres', 'or_sust', 'or_foss')
or.long$variable <- factor(or.long$variable, levels = ordered_levels)

p2 <- ggplot(or.long, aes(x = variable, y = value)) +
  geom_boxplot() +
  stat_summary(fun.data = stat_box_data, geom = "text", 
               hjust = 1, vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(title = expression(paste(italic("P. orientalis"))),
       x = "",
       y = "Prediction") +
  scale_x_discrete(labels = c('Present','Sustainable','Fossil-fueled')) + 
  coord_flip()
p2

p3 <- grid.arrange(p1, p2, ncol=2, nrow=1)
ggsave('results/prediction_boxplots.jpg',p3, width=10, height=6, units='in')

# variable importance - iberica

# load models filenames
file_list <- list.files("results/alchata", pattern = 'model', full.names = T)

# empty dfs to store accuracy and gini values for each model
accu <- data.frame(matrix(nrow = 27, ncol = 0))
# iterate through the model list to extract importances
for (i in 1:length(file_list)) {
  
  model <- readRDS(file_list[[i]])
  col_name <- paste0('fold', i)
  
  accu <- cbind(accu, model$importance[, 3])
  names(accu)[ncol(accu)] <- col_name
}
write.csv2(accu, "results/iberica_accu.csv", row.names=F)

### coloured boxplots of variable importance

# first with accuracy

# transpose df
accut <- as.data.frame(t(accu))
# change display names
display_names <- c("slope" = "Slope",
                   "mean_temp" = "Mean Temperature",
                   "prec" = "Acc. Precipitation",
                   "forest" = "Forest",
                   "mean_temp_seas" = "Mean. Temp. Seasonality",
                   "prec_seas" = "Acc. Prec. Seasonality",
                   "elev" = "Elevation",
                   "linear_infrastr" = "Linear Infrastruc.",
                   "crop_low" = "Low Int. Crops",
                   "ext_perm_crop" = "Extensive Perm. Crops",
                   "for_shr_agric" = "Forest-Shrub-Agriculture",
                   "shrub" = "Shrubs",
                   "crop_med" = "Medium Int. Crops",
                   "settlements" = "Settlements",
                   "int_perm_crop" = "Intensive Perm. Crops",
                   "for_shr_grass" = "Forest-Shrub-Grass",
                   "grass_low" = "Low Int. Grass",
                   "crop_high" = "High Int. Crops",
                   "for_shr_crop" = "Forest-Shrub-Crop",
                   "water" = "Water Bodies",
                   "grass_med" = "Medium Int. Grass",
                   "mosaic_low" = "Low Int. Agriculture Mosaic",
                   "mosaic_med" = "Medium Int. Agriculture Mosaic",
                   "grass_high" = "High Int. Grass",
                   "bare" = "Bare Soil",
                   "for_shr_bare" = "Forest-Shrub-Bare",
                   "mosaic_high" = "High Int. Agriculture Mosaic")
# define variable categories (for color legend)
categorias <- list(climaticas = c("Mean Temperature", "Mean. Temp. Seasonality",
                                  "Acc. Precipitation", "Acc. Prec. Seasonality"),
                   topograficas = c("Elevation", "Slope"),
                   perturbacion = c("Linear Infrastruc."))
# replace names
colnames(accut) <- display_names[match(colnames(accut), names(display_names))]
# vector with all variable names
todas_variables <- colnames(accut)
# category for land uses
usos <- setdiff(todas_variables, unlist(categorias))
# calculate absolute means of each variable
accut_mean_abs <- accut %>%
  summarise_all(~ mean(abs(.))) %>%
  gather(key = "Variable", value = "Mean_Abs") %>%
  arrange(Mean_Abs) # order ascending
# define category of each variable
accut_mean_abs$categoria <- case_when(
  accut_mean_abs$Variable %in% categorias$climaticas ~ "Climate",
  accut_mean_abs$Variable %in% categorias$topograficas ~ "Topography",
  accut_mean_abs$Variable %in% categorias$perturbacion ~ "Anthropic",
  accut_mean_abs$Variable %in% usos ~ "Land Uses"
)
# vector for variables
ordered_vars <- accut_mean_abs$Variable
# transform to long format
accu_long_ordered <- accut %>%
  gather(key = "Variable", value = "Value") %>%
  mutate(Variable = factor(Variable, levels = ordered_vars))
# create category variable in long format df
accu_long_ordered <- accu_long_ordered %>%
  mutate(categoria = case_when(
    Variable %in% categorias$climaticas ~ "Climate",
    Variable %in% categorias$topograficas ~ "Topography",
    Variable %in% categorias$perturbacion ~ "Anthropic",
    TRUE ~ "Land Uses"  
  ))

# define colors
colores <- c("Climate" = "#95b8f6", 
             "Topography" = "#f9d99a", 
             "Anthropic" = "#dcd9f8",
             "Land Uses" = "#fa5f49")

# create boxplots
plot1 <- ggplot(accu_long_ordered, aes(x = Variable, y = abs(Value), fill = categoria)) +
  geom_boxplot() +
  coord_flip() +  # turn labels in x axis
  theme_minimal() +  
  labs(title = "Variable importance - P. alchata",
       x = "Variable", y = "Mean Decrease in Accuracy", fill = "Category") +
  scale_fill_manual(values = colores) +  # assign defined colors
  theme(legend.position = "left")
print(plot1)

# variable importance - ortega

# load models filenames
file_list <- list.files("results/orientalis", pattern = 'model', full.names = T)

# empty dfs to store accuracy and gini values for each model
accu <- data.frame(matrix(nrow = 27, ncol = 0))
# iterate through the model list to extract importances
for (i in 1:length(file_list)) {
  
  model <- readRDS(file_list[[i]])
  col_name <- paste0('fold', i)
  
  accu <- cbind(accu, model$importance[, 3])
  names(accu)[ncol(accu)] <- col_name
}
write.csv2(accu, "results/ortega_accu.csv", row.names=F)

### coloured boxplots of variable importance

# first with accuracy

# transpose df
accut <- as.data.frame(t(accu))
# replace names
colnames(accut) <- display_names[match(colnames(accut), names(display_names))]
# vector with all variable names
todas_variables <- colnames(accut)
# category for land uses
usos <- setdiff(todas_variables, unlist(categorias))
# calculate absolute means of each variable
accut_mean_abs <- accut %>%
  summarise_all(~ mean(abs(.))) %>%
  gather(key = "Variable", value = "Mean_Abs") %>%
  arrange(Mean_Abs) # order ascending
# define category of each variable
accut_mean_abs$categoria <- case_when(
  accut_mean_abs$Variable %in% categorias$climaticas ~ "Climate",
  accut_mean_abs$Variable %in% categorias$topograficas ~ "Topography",
  accut_mean_abs$Variable %in% categorias$perturbacion ~ "Anthropic",
  accut_mean_abs$Variable %in% usos ~ "Land Uses"
)
# vector for variables
ordered_vars <- accut_mean_abs$Variable
# transform to long format
accu_long_ordered <- accut %>%
  gather(key = "Variable", value = "Value") %>%
  mutate(Variable = factor(Variable, levels = ordered_vars))
# create category variable in long format df
accu_long_ordered <- accu_long_ordered %>%
  mutate(categoria = case_when(
    Variable %in% categorias$climaticas ~ "Climate",
    Variable %in% categorias$topograficas ~ "Topography",
    Variable %in% categorias$perturbacion ~ "Anthropic",
    TRUE ~ "Land Uses"  
  ))

# create boxplots
plot2 <- ggplot(accu_long_ordered, aes(x = Variable, y = abs(Value), fill = categoria)) +
  geom_boxplot() +
  coord_flip() +  # turn labels in x axis
  theme_minimal() +  
  labs(title = "Variable importance - P. orientalis",
       x = "Variable", y = "Mean Decrease in Accuracy", fill = "Category") +
  scale_fill_manual(values = colores) +  # assign defined colors
  theme(legend.position = "left")
print(plot2)

plot3 <- grid.arrange(plot1, plot2, nrow = 2, ncol = 1)
ggsave("results/both_varimp.jpg", plot3, width = 8, height = 9, units = "in")

# # predictions histogram
# 
# # iberica
# 
# dataset <- na.omit(read.csv2("data/iberica_PA_present.csv"))
# dataset$PA <- as.factor(dataset$PA)
# dataset_points <- st_as_sf(dataset, coords = c("X", "Y"), crs = 4326)
# 
# points <- readRDS("objects/abundance_iberica.RData")
# preds <- read.csv2("results/pres_preds_iberica.csv")
# preds$mean <- rowMeans(preds)
# preds_r <- raster("results/iberica_preds.asc")
# 
# points$preds <- raster::extract(preds_r, points)
# mapview(points, zcol='preds')
# hist(points$preds)
# 
# hist(preds$mean)
# 
# dataset_preds <- raster::extract(preds_r, dataset_points)
# top10 <- quantile(dataset_preds, probs = 0.9, na.rm = TRUE)
# down10 <- quantile(dataset_preds, probs = 0.1, na.rm = TRUE)
# above10 <- subset(dataset_preds, dataset_preds>down10)
# 
# hist(dataset_preds)
# abline(v = down10)
# abline(v = mean(above10))
# abline(v = top10)
# 
# down10 <- quantile(points_preds, probs = 0.1, na.rm = TRUE)
# 
# hist(points$preds)
# n_distinct(points_preds)
# 
# base <- preds_r
# base[] <- 1:ncell(preds_r)
# points$cell_number <- raster::extract(base, points)
# points_grouped <- points %>%
#   group_by(cell_number) %>%
#   summarize(preds = mean(preds),
#             Ngangas = sum(Ngangas),
#             geometry = first(geometry))
# 
# plot(points_grouped$preds, points_grouped$Ngangas)
# reg <- lm(points_grouped$Ngangas ~ points_grouped$preds)
# summary(reg)
# abline(reg, col = "red", lwd = 2)
# reg <- lm(log(points_grouped$Ngangas) ~ points_grouped$preds)
# summary(reg)
# plot(points_grouped$preds, log(points_grouped$Ngangas))
# abline(reg, col = "blue", lwd = 2)
# 
# corr <- cor.test(points_grouped$preds, points_grouped$Ngangas)
# corr$estimate

# most frequent land uses

# all territory

dataset <- na.omit(read.csv2("data/envdf_present.csv"))
dataset$PA <- as.factor(dataset$PA)

str(dataset)

binary_cols <- sapply(dataset, function(x) all(x %in% c(0, 1)))
subset_binary <- dataset[, binary_cols]

frequency_ones <- colSums(subset_binary[,1:20])
frequency_ones

frequency_percent_envdf <- frequency_ones/nrow(dataset) * 100
frequency_percent_envdf

# iberica
dataset <- na.omit(read.csv2("data/iberica_PA_present.csv"))
dataset$PA <- as.factor(dataset$PA)

str(dataset)

subset_data <- dataset[dataset$PA == 1, ]

binary_cols <- sapply(subset_data, function(x) all(x %in% c(0, 1)))
subset_binary <- subset_data[, binary_cols]

frequency_ones <- colSums(subset_binary[,1:20])
frequency_ones

frequency_percent_iberica <- frequency_ones/nrow(dataset[dataset$PA == 1, ]) * 100
frequency_percent_iberica

# ortega
dataset <- na.omit(read.csv2("data/ortega_PA_present.csv"))
dataset$PA <- as.factor(dataset$PA)

str(dataset)

subset_data <- dataset[dataset$PA == 1, ]

binary_cols <- sapply(subset_data, function(x) all(x %in% c(0, 1)))
subset_binary <- subset_data[, binary_cols]

frequency_ones <- colSums(subset_binary[,1:20])
frequency_ones

frequency_percent_ortega <- frequency_ones/nrow(dataset[dataset$PA == 1, ]) * 100
frequency_percent_ortega

# join in columns
frequencies <- as.data.frame(list(
  IP = frequency_percent_envdf,
  P.alchata = frequency_percent_iberica,
  P.orientalis = frequency_percent_ortega))

write.csv2(frequencies, "results/landuses_frequencies.csv", row.names=T)
