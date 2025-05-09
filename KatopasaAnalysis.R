#### WHAT'S IN A MOUNTAIN? ####
# This code helps you find out!

#Put your REAL working directory on the line below.
setwd("C:/Users/Isaac/Syncthing-Docs/Lab Meeting/Checklist/Katopasa")
source("checklistFunctions2.R")

# Load in these libraries
library("tidyverse")
library("tibble")
library("vegan")
library("reshape2")
library("googlesheets4")
library("leaflet")
library("googledrive")
library("cowplot")
library("elevatr")
library("divDyn")
library("terra")
library("sf")
library("readxl")

##SETUP - simply supply
#the NAME of your mountain
mtn <- "Katopasa"

#### STEP 1 - Data Import/Cleanup ####
info <- readRDS("folderIDs.rds") %>% filter(Mountain == mtn)
df <- read_xlsx("Katopasa_Specimens.xlsx") %>% dfSetup() %>% addGroups()

#### STEP 2 - Add groups/Look up Elevation####

#### OPTIONAL: Check elevatr accuracy ####

# Check elevational accuracy according to aws data from elevatr

#

elevs <- filter(df, !is.na(Latitude), !is.na(Uncertainty))[,c("Longitude","Latitude","Uncertainty","Elevation")]
elevs <- filter(df, !is.na(Latitude))[,c("Longitude","Latitude")]

elevs_sf <- st_as_sf(x = elevs, 
                     coords = c("Longitude","Latitude"),
                     crs = 4326)

elev_v <- get_elev_raster(elevs_sf, z = 14, src = "aws") %>% 
  rast()

buffered_points <- st_buffer(elevs_sf,elevs_sf$Uncertainty)
elevs$dem <- terra::extract(elev_v, elevs_sf, fun = max, ID = FALSE) %>%
  unlist
elevs$maxelev <- terra::extract(elev_v, buffered_points, fun = max, ID = FALSE) %>%
  unlist
elevs$minelev <- terra::extract(elev_v, buffered_points, fun = min, ID = FALSE) %>%
  unlist
elevs$dem_uncertainty <- (elevs$maxelev - elevs$minelev)/2

elevs$elev_diff <- abs(elevs$Elevation - elevs$dem)

# the collected elevation data correlate well with the data recovered from the DEM
ggplot(elevs)+
  geom_point(aes(x = Elevation, y= dem))

#lon-lat coordinate uncertainty versus the range of elevations in the DEM in that error radius
# lon-lat uncertainty is much greater than corresponding elevational uncertainty
ggplot(elevs)+
  geom_point(aes(x = Uncertainty, y= dem_uncertainty))

# neither uncertainty in the lat-lon coordinates nor the resulting range of 
# possibly correct elevations in that error radius correlate 
# with the actual differences between reported elevation and dem elevation
ggplot(elevs)+
  geom_point(aes(x = dem_uncertainty, y= elev_diff))

ggplot(elevs)+
  geom_point(aes(x = Uncertainty, y= elev_diff))
elev_v <- get_aws_points(elevs_sf, z = 14, units = "meters")

df$elevatr <- NA

df$elevatr[which(!is.na(df$Latitude))] <- elev_v[[1]]$elevation

df$elev_diff <- df$Elevation-df$elevatr



ggplot(df)+
  geom_density(aes(x = elev_diff))

ggplot(df)+
  geom_point(aes(x = Elevation, y= elevatr))

ggplot(df)+
  geom_point(aes(x = Uncertainty, y= abs(elev_diff)))

error_lm <- lm(Uncertainty ~ abs(elev_diff), df)

elev_lm <- lm(elevatr ~Elevation, df, na.action = na.omit)

#elevatr is biased to +1.6m vs our gps points
correctionfactor <- mean(df$elev_diff%>% na.exclude())

checkelevs <- filter(df, elev_diff < -50 | elev_diff >50) 
#there are some major discrepancies, but all are solidly below 700m

elev_lm <- lm(Elevation ~ dem, data = elevs)
elev_uncertainty <- lm(elev_diff ~ dem_uncertainty, data = elevs)
elev_lm <- lm(Elevation ~ elevatr, data = df)

#### Set up df with elevational bands and taxon groups ####
cutoffs <- c(700,1400)

df <- df %>% elevBands(bbs = cutoffs)



#### STEP 3 - Generate Plots ####
accumulationPlots <- accCurve(
  filter(df, group %in% c("frog","lizard","snake"), JAM_Number !=15267),
  info$colors[[1]]%>%unlist)

elevationPlot <- plotElev(
  filter(df, group %in% c("frog","lizard","snake"), JAM_Number !=15267),
  info$colors[[1]]%>%unlist)


alignedAccumulation <- align_plots(plotlist = list(accumulationPlots$taxonacc + theme(legend.position = "none"),
                                                   accumulationPlots$bandacc + theme(legend.position = "none")),
                                   align = "hv",
                                   axis = "blr")

accumulation <- plot_grid(alignedAccumulation[[1]], alignedAccumulation[[2]], 
                          rel_widths = c(3,2),
                          rel_heights = c(2,1))

ggsave(filename = paste0("Figures/", mtn, "_AccumulationPlot.png"),
       plot = accumulation,
       width = 8, height = 6, units = "in",
       bg = "white")

ggsave(filename = paste0("Figures/", mtn, "_ElevationPlot.png"),
       plot = elevationPlot$rangeThrough,
       width = 8, height = 7, units = "in",
       bg = "white")

drive_put(paste0("Figures", mtn, "_ElevationPlot.png"), path = as_id(info$folder)) 
drive_put(paste0("Figures", mtn, "_AccumulationPlot.png"), path = as_id(info$folder))

#### STEP 4 - Generate Tables ####

#Table 1
t1 <- speciesTable(filter(df, group %in% c("lizard","frog","snake")))
write_sheet(t1$publicationTable, ss = as_id(info$t1), sheet = 1)

#Table 2
t2 <- missingSpecies(filter(df,group %in% c("lizard","frog","snake"))) 
write_sheet(t2$PublicationTable, ss = as_id(info$t2), sheet = 1)



#### STEP 5 - Leaflet Interactive Map #####
elevationcolors <- colorFactor(c("#A8CCDE","#4E819A","#3B5374"), df$eband)

df$Habitat[which(is.na(df$Habitat))] <- "not recorded"

map <- df %>% 
  leaflet() %>%
  addProviderTiles(providers$Esri.WorldTopo) %>% # using ESRI World Topo for the background map tiles
  
  addCircleMarkers(label = paste(df$JAM_Number,"-",df$binom),
                   color = ~elevationcolors(df$eband),
                   popup = paste("<center>",df$JAM_Number, "-", df$binom, "<br>",
                                 "Collected ",df$Collection_Date, "</center>","<br>",
                                 "Elevation: ", df$Elevation, "m", "<br>",
                                 "Habitat: ", df$Habitat),
                   radius=4) %>% #this adds in markers. on mouse-over, it will say their JAM ID and the binomial
  
  setView(lat = df$Latitude[10], lng = df$Longitude[10],zoom = 11) %>%
  addScaleBar(position = "bottomleft",
              options= scaleBarOptions(metric = TRUE))

library(htmlwidgets)
saveWidget(map, file = "Figures/Katopasa_leaflet.html")
