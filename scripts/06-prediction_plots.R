
rm(list=ls())
source('scripts/utils.R') # all used libraries and my custom functions

# spanish IP borders
spain <- getpeninsularspain()

full <- as.data.frame(stack(raster('results/iberica_preds.asc'),
                            raster('results/ortega_preds.asc'),
                            raster('results/iberica_ssp126_preds.asc'),
                            raster('results/ortega_ssp126_preds.asc'),
                            raster('results/iberica_ssp585_preds.asc'),
                            raster('results/ortega_ssp585_preds.asc')), xy=T)

colnames(full)[3:8] <- c('ib_pres', 'or_pres', 'ib_sust', 'or_sust', 'ib_foss', 'or_foss')

p1 <- ggplot(full, aes(x = x, y = y, fill = ib_pres)) +
  geom_raster() +
  scale_fill_viridis_c(na.value='transparent',
                     limits = c(0, 1),  # Set the limits of the scale
                     oob = scales::squish) +  # Handle out-of-bounds values
  labs(x = 'Longitude', y = 'Latitude',
       title = expression(paste('a) ', italic("P. alchata"), " - Present"))) +
  theme_bw() +
  theme(legend.position = 'none') + 
  geom_sf(data = spain, fill = 'transparent', color = 'black', inherit.aes = FALSE)
p1
# deactivate theme to extract legend and save for paper plot
# ggsave('results/getcontlegend.jpg', p1, height=15, width=14)

p2 <- ggplot(full, aes(x = x, y = y, fill = or_pres)) +
  geom_raster() +
  scale_fill_viridis_c(na.value='transparent',
                       limits = c(0, 1),  # Set the limits of the scale
                       oob = scales::squish) +  # Handle out-of-bounds values
  labs(x = 'Longitude', y = 'Latitude',
       title = expression(paste('b) ', italic("P. orientalis"), " - Present"))) +
  theme_bw() +
  theme(legend.position = 'none') + 
  geom_sf(data = spain, fill = 'transparent', color = 'black', inherit.aes = FALSE)
p2

p3 <- ggplot(full, aes(x = x, y = y, fill = ib_sust)) +
  geom_raster() +
  scale_fill_viridis_c(na.value='transparent',
                       limits = c(0, 1),  # Set the limits of the scale
                       oob = scales::squish) +  # Handle out-of-bounds values
  labs(x = 'Longitude', y = 'Latitude',
       title = expression(paste('c) ', italic("P. alchata"), " - Sustainable Future"))) +
  theme_bw() +
  theme(legend.position = 'none') + 
  geom_sf(data = spain, fill = 'transparent', color = 'black', inherit.aes = FALSE)
p3

p4 <- ggplot(full, aes(x = x, y = y, fill = or_sust)) +
  geom_raster() +
  scale_fill_viridis_c(na.value='transparent',
                       limits = c(0, 1),  # Set the limits of the scale
                       oob = scales::squish) +  # Handle out-of-bounds values
  labs(x = 'Longitude', y = 'Latitude',
       title = expression(paste('d) ', italic("P. orientalis"), " - Sustainable Future"))) +
  theme_bw() +
  theme(legend.position = 'none') + 
  geom_sf(data = spain, fill = 'transparent', color = 'black', inherit.aes = FALSE)
p4

p5 <- ggplot(full, aes(x = x, y = y, fill = ib_foss)) +
  geom_raster() +
  scale_fill_viridis_c(na.value='transparent',
                       limits = c(0, 1),  # Set the limits of the scale
                       oob = scales::squish) +  # Handle out-of-bounds values
  labs(x = 'Longitude', y = 'Latitude',
       title = expression(paste('e) ', italic("P. alchata"), " - Fossil-fueled Future"))) +
  theme_bw() +
  theme(legend.position = 'none') + 
  geom_sf(data = spain, fill = 'transparent', color = 'black', inherit.aes = FALSE)
p5

p6 <- ggplot(full, aes(x = x, y = y, fill = or_foss)) +
  geom_raster() +
  scale_fill_viridis_c(na.value='transparent',
                       limits = c(0, 1),  # Set the limits of the scale
                       oob = scales::squish) +  # Handle out-of-bounds values
  labs(x = 'Longitude', y = 'Latitude',
       title = expression(paste('f) ', italic("P. orientalis"), " - Fossil-fueled Future"))) +
  theme_bw() +
  theme(legend.position = 'none') + 
  geom_sf(data = spain, fill = 'transparent', color = 'black', inherit.aes = FALSE)
p6

p7 <- grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
ggsave('results/all_cont_preds_separated.jpg', p7, height=15, width=14)
