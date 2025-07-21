library(osmdata)
library(tidyverse)
library(sf)
library(ggOceanMaps)
library(ggpubr)

dt <- expand.grid(lon = c(-170, 20), lat = c(15, 25))
sampling_hawaii <- data.frame(lon = c(-158, -157.8), 
                 lat = c(22.45, 21.5), 
                 var = c("Station ALOHA", "Kane'ohe Bay")
)


map_hawaii = basemap(data = dt, bathymetry = TRUE, bathy.style = "rcb", legends = FALSE, lon.interval=5) + 
  geom_point(data = transform_coord(sampling_hawaii), aes(x = lon, y = lat), color = "red") + 
  geom_text(data = sampling_hawaii, aes(x = lon, y = lat, label = var), hjust=-0.1, vjust=-0.1) +
  coord_sf(xlim=c(-162, -154),
           ylim=c(18.5, 24),
           expand=FALSE) 
ggsave('map_hawaii.pdf')