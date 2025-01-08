rm(list=ls())
source('scripts/utils.R')

# show barplot of installed energy per CCAA

solar.cont <- raster('spatial_data/solar_energy/solar_continuous.asc')
crs(solar.cont)

spain <- mapSpain::esp_get_ccaa()[!mapSpain::esp_get_ccaa()$iso2.ccaa.name.es %in%
                                    c('Islas Baleares', 'Canarias', 'Ceuta', 'Melilla'), ] %>%
  st_transform(crs(solar.cont))

resultados <- exact_extract(solar.cont, spain, fun = "sum")

# add to map
spain$solar_sum <- resultados

# barplot
ggplot(spain, aes(x = iso2.ccaa.name.es, y = solar_sum)) +
  geom_bar(stat = "identity", fill = "steelblue") +  
  labs(title = "Solar count per CCAA",
       x = "CCAA",
       y = "Aprox. km2") +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # rotate names

# cor between pretor database and our data

# load all pretor ccaa and bind
filepaths <- list.files('spatial_data/solar_energy/PRETOR', pattern='csv', full.names=T)
files <- lapply(filepaths, read.csv2)
pretor <- do.call(rbind, files)

# # to check cor with the data we used to georreference use this xlsx instead (the old one)
# pretor <- read_xlsx('spatial_data/solar_energy/PRETOR/energy_all_provinces-19022024.xlsx')

# load municipes sf
nuts <- st_read('spatial_data/boundaries/SHP_ETRS89/recintos_municipales_inspire_peninbal_etrs89/recintos_municipales_inspire_peninbal_etrs89.shp')

# data cleaning, different municipe names between pretor and nuts
# find differences between names
unique_names <- setdiff(pretor$Municipio.de.la.Instalación, nuts$NAMEUNIT)
print(unique_names)

municipios <- pretor$Municipio.de.la.Instalación

# transform 'Ejido (El)' to 'El Ejido' 
# this line first finds any pattern with ^(.*?), then adds a space,
# then finds anything between parenthesis and stops searching
# then takes the second object and then the first object and put them together
primera_correccion <- sub("^(.*?) \\((.*?)\\)$", "\\2 \\1", municipios)

unique_names <- setdiff(primera_correccion, nuts$NAMEUNIT)
print(unique_names) # check, still many differences and several patterns

# we're using the function stringdist to find best coincidences in the nuts names
?stringdist

find_best_match <- function(target, source) {
  similarities <- stringdist::stringdist(target, source, method = "osa")
  best_match_index <- which.min(similarities)
  best_match <- source[best_match_index]
  return(best_match)
}

# vector to store new names
updated_names <- vector("character", length = length(primera_correccion))

# iterate through pretor names and find best coincidences in nuts
for (i in seq_along(primera_correccion)) {
  updated_names[i] <- find_best_match(primera_correccion[i], nuts$NAMEUNIT)
}

prueba <- setdiff(updated_names, nuts$NAMEUNIT)
print(prueba) # check
# bind new names
pretor$NAMEUNIT <- updated_names

# create PV database by municipality 
PV_temp <- pretor %>%
  filter(Grupo.Normativo %in%  c("b.1.1", "b.1.1.I", "b.1.1.I.1", "b.1.1.I.2", "b.1.1.II", "b.1.2"))

PV <- PV_temp %>% 
  group_by(NAMEUNIT, Tipo.de.Inscripción) %>% 
  summarize(MW = sum(Potencia.Instalada.KW/1000))

PV_def <- PV %>% filter(Tipo.de.Inscripción == "DEFINITIVA")
PV_prev <- PV %>% filter(Tipo.de.Inscripción == "PREVIA")

# join nuts (sf) with the corresponding MW
nuts <- left_join(nuts, PV_def, by='NAMEUNIT')
mapview(nuts, zcol='MW')

# now we need to extract the total areas calculated in our raster 
nuts <- st_transform(nuts, crs=crs(solar.cont))
nuts$solar.cont <- exact_extract(solar.cont, nuts, fun = "sum")

# fix the wrong assignation of MW to Cieza (Cantabria), which corresponds to Cieza (Murcia)
nuts[nuts$NAMEUNIT == 'Cieza' & nuts$CODNUT3 == 'ES130',]$MW <- 0

# replace NAs in MW with zeros
nuts$MW[is.na(nuts$MW)] <- 0

write_sf(nuts, 'spatial_data/solar_energy/PRETOR/nuts_with_mw_25092024.shp')
# nuts <- st_read('spatial_data/solar_energy/PRETOR/nuts_with_mw_25092024.shp')

# and cor

corrdata <- nuts[,c('MW','solar.cont')] %>% st_drop_geometry()
cor_matrix <- cor(corrdata, use = "complete.obs")
print(cor_matrix)
# super high cor (>0.9)

logdata <- as.data.frame(cbind(log10(corrdata[,1]+0.001), log10(corrdata[,2]+0.001)))

p1 <- ggplot(corrdata, aes(x = MW, y = solar.cont)) +
  geom_point() +
  geom_smooth(method='lm',se=F) +
  labs(x = "PRETOR MW per municipe", y = "km2 per municipe (our obs)",
       title = 'Original') +
  theme_minimal()

p2 <- ggplot(logdata, aes(x = V1, y = V2)) +
  geom_point() +
  geom_smooth(method='lm',se=F) +
  labs(x = "PRETOR MW per municipe", y = "km2 per municipe (our obs)",
       title = 'log10') +
  theme_minimal()

p3 <- grid.arrange(p1, p2, ncol=2)
ggsave('results/current_pretor_vs_ourdata.jpg', p3, width=7, height=5)
