#SulawesiMapper
setwd("C:/Users/Isaac/Syncthing-Docs/Lab Meeting/Checklist/Katopasa")

library('sf')
library('elevatr')
library('raster')
library('terra')
library('tidyterra')
library('ggplot2')
library('ggnewscale')
library('rnaturalearth')
library('osmdata')
library('ggspatial')

#### Annotations dataset #### ○

annotations <- data.frame(lat = -1.11952,	lon = 121.45214, label = "Desa Mire (210 m)", type = "\u229A") %>%
  rbind(list(lat = -1.1296, lon = 121.44535, label = "Low Camp (375 m)", type = "\u229A")) %>%
  rbind(list(lat = -1.14952, lon =	121.44768, label = "Mid Camp (750 m)", type = "\u229A")) %>%
  rbind(list(lat = -1.18458, lon = 121.44527, label = "High Camp (1366 m)", type = "\u229A")) %>%
  rbind(list(lat = -1.24025, lon = 121.42624, label = "Gunung Katopasa Summit (2825 m)", type = "\u25b2")) %>%
  rbind(list(lat = -1.21401, lon = 121.43414, label = "First Summit (2575 m)", type = "\u25b2"))

annotations$lon2 <- annotations$lon + 0.01

#### KATOPASA MAP ####
extent.df <- data.frame(
  lat = c(-1.02,-1.02,-1.25,-1.25),
  lon = c(121.35,121.52 ,121.35, 121.52)
) %>% st_as_sf(coords = c("lon","lat")) %>% st_set_crs(4326) 

#load in forestwatch data
<<<<<<< HEAD
# Forest Watch data granules (Top Left Corner at 0n, 120E) available at: 
# https://storage.googleapis.com/earthenginepartners-hansen/GFC-2020-v1.8/download.html

=======
>>>>>>> bc4f5e8fa668f6b1184142c146ac765bed51d321
fw <- terra::rast("Hansen_GFC-2020-v1.8_treecover2000_00N_120E.tif") %>%
  crop(extent.df)

loss <- terra::rast("Hansen_GFC-2020-v1.8_lossyear_00N_120E.tif") %>%
  crop(extent.df)

loss_abs <- app(loss, function(x){x==0})

loss_15_16 <- loss %>% 
  subst(from = 0:14, to = NA) %>%
  subst(from = 17:20, to = NA)


current_trees <- (fw*loss_abs) %>% subst(from = 0:1, to = NA)
  
#Get raster elevation data for basemap
elev_r <- get_elev_raster(extent.df, z = 12, src = "aws") %>% 
  rast()

elev_c <- as.contour(elev_r, levels = c(700,1400,2100,2800)) %>% crop(extent.df)
elev_c5 <- as.contour(elev_r, levels = c(1:28)*100) %>% crop(extent.df)
elev_c$level <- as.factor(elev_c$level)
names(elev_c) <- c("elevation")

terrain_slope <- terrain(elev_r, v="slope", unit="radians", neighbors=8)
terrain_aspect <- terrain(elev_r, v="aspect", unit="radians", neighbors=8)

hillshade <- shade(terrain_slope,terrain_aspect, angle = 45, direction = 25, normalize = TRUE)

ocean_mask <- get_elev_raster(extent.df, src = "gl3") %>% rast() %>% 
  subst(from = 0:5, to = NA) %>% resample(elev_r) %>% 
  crop(extent.df)

hillshade <- resample(hillshade, elev_r) %>% 
  crop(extent.df) %>% 
  mask(ocean_mask)
# Get vector river & roads from basemap

river <- opq(st_bbox(extent.df)) %>%
  add_osm_feature(key = 'waterway') %>% osmdata_sf()

river_v <- river$osm_lines %>% st_crop(extent.df)

ocean <- opq(st_bbox(extent.df)) %>%
  add_osm_feature(key = 'natural', value = 'coastline') %>% osmdata_sf()
ocean_v <- ocean$osm_polygons %>% st_crop(extent.df)

roads <- opq(st_bbox(extent.df)) %>%
  add_osm_feature(key = 'highway') %>% osmdata_sf()
roads_v <- roads$osm_lines %>% st_crop(extent.df)

# Convert dataset for plotting

df_sf <- st_as_sf(df %>% filter(!is.na(Longitude), group != "turtle"), coords=c("Longitude","Latitude")) %>% 
  st_set_crs(4326)

# set water color
water <- "#aad3ff"

map <- ggplot() +
  
  # The shadows part of the raster map
  geom_spatraster(
    data = hillshade
  ) +
  paletteer::scale_fill_paletteer_c(
    "grDevices::Light Grays",
    na.value = water,
    direction = 1,
    guide = "none"
  ) +
  
  geom_sf(data = elev_c5, inherit.aes = FALSE, lwd=0.2, color = "#00000044")+
  # New fill scale
  ggnewscale::new_scale_fill() +
  
  # The elevation digital raster
  # geom_spatraster(
  #   data = elev_r,
  #   alpha = 0.7
  # ) +
  # scale_fill_hypso_c(
  #   palette = "dem_print",
  #   na.value = "#6688ff",
  #   breaks = c(700,1400),
  #   limits = c(0,5000)
  # ) +
  # 
  # #The elevation bands
  # geom_spatraster(
  #   data = elev_r
  # ) +
  # scale_fill_binned(
  #   low = "#77553300",
  #   high = "#88665500",
  #   na.value = water,
  #   breaks = c(700,1400),
  #   limits = c(0,5000),
  #   guide = "none"
  # ) +
  # 
  ggnewscale::new_scale_fill() +
  geom_spatraster(data = current_trees) +
  scale_fill_gradient(low = "#ffffff00", high = "#88aa3366", 
                      na.value = "#00000000",
                      guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_spatraster(data = loss_15_16)+
  scale_fill_gradient(low = "#11111166", high = "#11111166", 
                      na.value = "#00000000",
                      guide = "none") +
  geom_sf(data = elev_c, inherit.aes = FALSE, lwd=0.3, color = "#000000")+
  geom_sf(data = ocean_v, inherit.aes = FALSE, lwd=1, color = water)+
  geom_sf(data = river_v, inherit.aes = FALSE, lwd=1, color = water)+
  geom_sf(data = roads_v, inherit.aes = FALSE, lwd=1, color = "#555")+
  # New fill scale
  ggnewscale::new_scale_fill() +
  geom_sf(data = df_sf, inherit.aes = FALSE, aes(shape = group, fill = group), size = 3, color = "black")+
  scale_shape_manual(values = c(24,21,25)) +
  scale_fill_manual(values =  info$colors%>%unlist)+
  #geom_sf(data = df_sf, inherit.aes = FALSE, aes(color = group, shape = group), size = 3)+
  
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  annotation_scale()+
  geom_label(data = annotations, aes(x = lon2, y = lat, label = label),
             hjust = 0, size = 3)+
  annotate("segment", x = annotations$lon, y = annotations$lat, 
           xend = annotations$lon2, yend = annotations$lat)+
  annotate("text", x = annotations$lon[1:4], y = annotations$lat[1:4], label = annotations$type[1:4], size = 15)+
  annotate("text", x = annotations$lon[5:6], y = annotations$lat[5:6], label = annotations$type[5:6])+
  labs(elevation = "elevation (m)", group = "taxonomic group")+
  guides(fill = "legend", shape = "legend") +
  theme(legend.justification.right = "bottom",
        legend.title = element_text("animal"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


#### SULAWESI MAP ####
katopasa.box <- ext(extent.df) %>% as.polygons() 
crs(katopasa.box)<-crs(extent.df)

sulawesi.df <- data.frame(
  lat = c(3,3,-6,-6),
  lon = c(118.5,126 ,118.5, 126)
) %>% st_as_sf(coords = c("lon","lat")) %>% st_set_crs(4326) 

sulawesi_map <- get_elev_raster(sulawesi.df, src = "gl3") %>% rast()

sw_slope <- terrain(sulawesi_map, v="slope", unit="radians", neighbors=8)
sw_aspect <- terrain(sulawesi_map, v="aspect", unit="radians", neighbors=8)

sw_hillshade <- shade(sw_slope,sw_aspect, angle = 45, direction = 25, normalize = TRUE)

sulawesi_inset <- ggplot() +
  
  # The shadows part of the raster map
  geom_spatraster(
    data = sw_hillshade
  ) +
  paletteer::scale_fill_paletteer_c(
    "grDevices::Light Grays",
    direction = 1,
    na.value = water,
    guide = "none"
  ) +
  geom_sf(data=katopasa.box, fill = "#88aa3366")+
  scale_x_continuous(expand = c(0, 0), 
                     breaks = c(120,122,124), 
                     labels = c("120°E","122°E","124°E")) +
  scale_y_continuous(expand = c(0, 0))
  

library('cowplot')

ggdraw(map)+ draw_plot(sulawesi_inset, x = 0.6, y = 0.6, height = 0.4, width = 0.4)

<<<<<<< HEAD
ggsave("Figures/KatopasaMap.png",width = 8, height = 8, units = "in")  
=======
ggsave("KatopasaMap.png",width = 8, height = 8, units = "in")  
>>>>>>> bc4f5e8fa668f6b1184142c146ac765bed51d321
