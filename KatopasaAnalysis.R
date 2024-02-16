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


##SETUP - simply supply
# 2 - the NAME of your mountain
mtn <- "Katopasa"

#### STEP 1 - Data Import/Cleanup ####
info <- readRDS("folderIDs.rds") %>% filter(Mountain == mtn)
df <- read_sheet(info$datasheet) %>% dfSetup() %>% addGroups()

#### STEP 2 - Add groups/Look up Elevation####
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

ggsave(filename = paste0(mtn, "_AccumulationPlot.png"),
       plot = accumulation,
       width = 8, height = 6, units = "in",
       bg = "white")

ggsave(filename = paste0(mtn, "_ElevationPlot.png"),
       plot = elevationPlot$rangeThrough,
       width = 8, height = 6.5, units = "in",
       bg = "white")

drive_put(paste0(mtn, "_ElevationPlot.png"), path = as_id(info$folder)) #not working
drive_put(paste0(mtn, "_AccumulationPlot.png"), path = as_id(info$folder))

#### STEP 4 - Generate Tables ####

#Table 1
t1 <- speciesTable(filter(df, group %in% c("lizard","frog","snake")))
write_sheet(t1$publicationTable, ss = as_id(info$t1), sheet = 1)

#Table 2
t2 <- missingSpecies(filter(df,group %in% c("lizard","frog","snake"))) 
write_sheet(t2$PublicationTable, ss = as_id(info$t2), sheet = 1)



#### STEP 5 - Calculate beta diversity for elevational bands

beta <- betadiver(veganize(df, elevation = TRUE))

#### STEP 6 - Leaflet Interactive Map #####
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

