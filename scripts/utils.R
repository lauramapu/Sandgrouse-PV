library(dplyr)
library(mapview)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggspatial) 
library(ggrepel)
library(sf)
library(units)
library(mapSpain)
library(cowplot)   # combine ggplot2 objects
library(gridExtra)
library(raster)
library(MazamaSpatialUtils)
library(biscale)
library(grDevices) # to create custom color palette
library(biscale)
library(terra)
library(purrr)
library(randomForest)
library(caret)
library(pdp)
library(dplyr)
library(reshape2) # to melt data for graphics
library(tictoc)
library(pROC) # roc curve, calculate thresholds
library(tidyr)
library(modEvA) # boyce index
library(car) # vif
library(readxl)
library(stringdist)
library(exactextractr)
library(ggpubr)

# function to list and cite all the above packages
versionandcitation <- function() {
  # list
  packages <- c('dplyr', 'mapview', 'ggplot2', 'RColorBrewer', 'viridis',
                'ggspatial', 'ggrepel', 'sf', 'units', 'mapSpain', 
                'cowplot', 'gridExtra', 'raster', 'MazamaSpatialUtils', 
                'biscale', 'grDevices', 'terra', 'purrr', 'randomForest', 
                'caret', 'pdp', 'reshape2', 'tictoc', 'pROC', 'tidyr',
                'modEvA', 'car', 'readxl', 'stringdist', 'exactextractr')
  
  # format the reference
  format_citation <- function(citation_info) {
    title <- citation_info$title
    authors <- sapply(citation_info$author, function(a) paste(a$given, a$family))
    year <- citation_info$year
    note <- citation_info$note
    url <- citation_info$url
    
    citation_text <- paste0(title, '. ', paste(authors, collapse = ', '), '. ', year, '. ', note, '. ', url)
    return(citation_text)
  }
  
  # versions and references
  package_info <- lapply(packages, function(pkg) {
    version <- as.character(packageVersion(pkg))
    citation_info <- citation(pkg)
    citation_text <- format_citation(citation_info)
    list(Package = pkg, Version = version, Citation = citation_text)
  })
  
  # to df
  package_info_df <- do.call(rbind, lapply(package_info, as.data.frame))
  rownames(package_info_df) <- NULL
  
  # print
  print(package_info_df)
  
  # save to txt
  write.table(package_info_df, file = 'R-packages.txt', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# function to load spanish borders
getpeninsularspain <- function() {
  provinces <- mapSpain::esp_get_prov()[!mapSpain::esp_get_prov()$iso2.prov.name.es %in%
                                          c('Las Palmas', 'Santa Cruz de Tenerife', 'Baleares', 'Melilla', 'Ceuta'), ]
  provinces$campo <- 1
  spain <- provinces %>%
    dissolve(field='campo') %>%
    st_transform(crs=4326)
}

# function to construct bivariated maps
bimap <- function(biclass_object, title, xlab, ylab, mypal, dim) {
  # the plot itself
  p <- ggplot() +
    geom_raster(data = biclass_object, aes(x = x, y = y, fill = bi_class)) +
    bi_scale_fill(pal = mypal, dim=dim, na.value='transparent') +
    labs(x = 'Longitude', y = 'Latitude', fill = '') +
    ggtitle(title) +
    theme_bw() + 
    theme( # legend modifications
      legend.position = 'none',  
      legend.background = element_rect(fill = 'transparent'),  # transparent background
      legend.key = element_rect(fill = 'transparent', color = NA)  # transparent background for legend keys
    ) +
    geom_sf(data = spain, fill = 'transparent', color = 'black', inherit.aes=F)
  p
  
  # bivariate legend 
  legend <- bi_legend(pal = mypal,
                      dim = dim,
                      xlab = xlab,
                      ylab = ylab,
                      size = 10,
                      arrows = T,
                      pad_width = 1.5) # distance between colors
  legend <- legend +
    theme_void() +
    theme(
      panel.background = element_rect(fill = 'transparent', color = NA),
      plot.background = element_rect(fill = 'transparent', color = NA),
      axis.title.x = element_text(size = 10, color = 'black'),
      axis.title.y = element_text(size = 10, color = 'black', angle = 90), # vertical text
      axis.text = element_blank(),  # hide exes values
      axis.ticks = element_blank()  # hide exes values
    )
  legend
  
  # combine plot and legend in same plot
  combined_plot <- ggdraw() +
    draw_plot(p) +
    draw_plot(legend, x = 0.68, y = 0.11, width = 0.25, height = 0.25)
  return(combined_plot)
}

# # function to assign quintile to each cell basing on the three distributions
# labs_species <- function(row) {
#   # how many columns have a 5
#   count_5 <- sum(row == '5', na.rm = TRUE)
#   count_all <- sum(as.numeric(row))
#   # assign value to new column basing on the count of 5s and total count
#   if (any(is.na(row))) {
#     return(NA)
#   } else if (count_5 == 3) {
#     return(5)
#   } else if (count_5 == 2 & count_all >= 13) {
#     return(4)
#   } else if (count_5 == 2 & count_all < 13) {
#     return(3)
#   } else if (count_5 == 1) {
#     return(2)
#   } else {
#     return(1)
#   }
# }

# get legend from plot, extract it and place it over
get_leyend <- function(plot, x = 0.68, y = 0.11, width = 0.25, height = 0.25) {
  # Function to extract legend
  get_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == 'guide-box')
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  # Extract legend
  legend <- get_legend(plot)
  
  # Remove legend from original plot
  plot_no_legend <- plot + theme(legend.position = 'none')
  
  # Create a new ggplot object for the legend with transparent background
  legend_plot <- ggplot() +
    theme_void() +
    theme(
      legend.position = 'center',
      legend.box.background = element_rect(fill = 'transparent', color = NA),
      legend.background = element_rect(fill = 'transparent', color = NA),
      panel.background = element_rect(fill = 'transparent', color = NA),
      plot.background = element_rect(fill = 'transparent', color = NA)
    ) +
    annotation_custom(legend)
  
  # Combine the plot and legend
  combined_plot <- ggdraw() +
    draw_plot(plot_no_legend) +
    draw_plot(legend_plot, x = x, y = y, width = width, height = height)
  
  return(combined_plot)
}

# # same function as previous but for both species together
# lab_global <- function(row) {
#   # how many columns have a 5
#   count_5 <- sum(row == '5', na.rm = TRUE)
#   
#   # assign value to new column basing on the count of 3s
#   if (any(is.na(row))) {
#     return(NA)
#   } else if (count_5 == 6) {
#     return(5)
#   } else if (count_5 == 5) {
#     return(4)
#   } else if (count_5 == 4) {
#     return(3)
#   } else if (count_5 == 3) {
#     return(2)
#   } else {
#     return(1)
#   }
# }

# custom color palette for biscale
custom_bipal <- function(min.xy, max.y, max.x, max.xy, dim) {
  
  # convert colors to RGB values
  tl <- col2rgb(min.xy)
  tr <- col2rgb(max.y)
  bl <- col2rgb(max.x)
  br <- col2rgb(max.xy)
  
  # function to interpolate colors
  interpolate_color <- function(tl, tr, bl, br, x, y) {
    top <- (1 - x) * tl + x * tr
    bottom <- (1 - x) * bl + x * br
    color <- (1 - y) * top + y * bottom
    return(rgb(color[1], color[2], color[3], maxColorValue = 255))
  }
  
  # initialize matrix to store colors
  color_matrix <- matrix(NA, nrow = dim, ncol = dim)
  
  # generate color grid
  for (i in seq_len(dim)) {
    for (j in seq_len(dim)) {
      x <- (i - 1) / (dim - 1)
      y <- (j - 1) / (dim - 1)
      color_matrix[i, j] <- interpolate_color(tl, tr, bl, br, x, y)
    }
  }
  
  # generate names vector based on dim
  custom_pal_names <- as.vector(outer(1:dim, 1:dim, function(x, y) paste(x, y, sep = '-')))
  
  # convert matrix to vector
  color_vector <- as.vector(t(color_matrix))
  names(color_vector) <- custom_pal_names
  
  return(color_vector)
}

# function to plot together total cells per tenth and solar infrastructures presence per tenth
plot_counts <- function(df, mytitle) {
  # suitable cells
  plot1 <- ggplot() +
    geom_line(data = na.omit(df), 
              aes(x = labs, y = n, 
                  group = 1)) +
    scale_y_continuous(name = 'Cell count', limits=c(0, max_n)) +
    labs(x = 'Suitability probability', 
         subtitle = mytitle) +
    theme_minimal() +
    theme(legend.position = 'none')
  # solar overlap
  plot2 <- ggplot() +
    geom_line(data = na.omit(df), 
              aes(x = labs, y = solar.bin, 
                  group = 1)) +
    scale_y_continuous(name = 'Cell count', limits=c(0, max_solar)) +
    labs(x = 'Solar presence per suit. prob.') +
    theme_minimal() +
    theme(legend.position = 'none')
  
  return(grid.arrange(plot1, plot2, ncol=1, nrow=2))
}

# calculate solar overlap in each quintile
calc_overlap <- function(distri) {
  both %>%
    group_by(category = !!sym(distri)) %>%
    summarise(
      total_count = n(),
      solar_count = sum(solar == 1)
    ) %>%
    mutate(
      percentage = (solar_count / total_count) * 100,
      column_name = distri
    ) %>%
    mutate(category = as.character(category))
}

# function to add labels to the overlaps table
add_labels <- function(df) {
  df %>%
    mutate(
      group_species = case_when(
        grepl('ib', column_name) ~ 'P. alchata',   
        grepl('or', column_name) ~ 'P. orientalis', 
        TRUE ~ 'Both'
      ),
      group_energy = case_when(
        grepl('pres', column_name) ~ 'Present',       
        grepl('sust', column_name) ~ 'Sustainable',  
        grepl('foss', column_name) ~ 'Fossil-fueled', 
        TRUE ~ 'All'
      )
    )
}

# functions to resample from 1x1 to 5x5 by mean and max

resample.avg <- function(x, colname) {
  df <- rasterFromXYZ(data.frame(x[,c('x','y')], x[,colname]),
                      res = c(0.008333333,0.008333333), crs = 4326) %>%
    aggregate(fact=5, fun='mean') %>%
    as.data.frame(xy=T)
  colnames(df) <- c('x','y',paste0(colname))
  return(df)
}

resample.max <- function(x, colname) {
  df <- rasterFromXYZ(data.frame(x[,c('x','y')], x[,colname]),
                      res = c(0.008333333,0.008333333), crs = 4326) %>%
    aggregate(fact=5, fun='max') %>%
    as.data.frame(xy=T)
  colnames(df) <- c('x','y',paste0(colname))
  return(df)
}