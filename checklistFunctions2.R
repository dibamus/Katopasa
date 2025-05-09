#Checklist functions that accommodate an arbitrary number of elevational bands
accCurve <- function(df, colors){
  
  df <- df[order(df$Collection_Date),] #reorder df by Collection Date
  
  #We'll be referencing elevational bands a lot here, so extract them to variables
  bands <- levels(df$eband)
  
  #figure out what dates & days to use for plotting & organization
  dates <- unique(na.omit(df$Collection_Date))
  banddays <- c(1:max(na.omit(df$dayinband))+1) #adding an extra day so the code for removing points after band sampling works
  
  
  countSp <- function(x,df = df,tx){(filter(df,Collection_Date==x,group==tx,firstocc) %>% dim)[1]}
  countSpEl <- function(x,df = df,el){(filter(df,dayinband == x,eband==el,firstinband) %>% dim)[1]}
  
  #make data frame of new species counts per taxon per day
  newsp <- cbind(date = dates,
                 frog =sapply(dates,countSp,df=df,tx = 'frog'),
                 lizard =sapply(dates,countSp,df=df,tx = 'lizard'),
                 snake =sapply(dates,countSp,df=df,tx = 'snake')) %>%
    as.data.frame(row.names = as.character(dates)) %>%
    mutate(cufrog = cumsum(frog),
           culizard = cumsum(lizard),
           cusnake = cumsum(snake)) %>%
    mutate(total = frog+lizard+snake) %>%
    mutate(cumulative = cumsum(total))%>%
    transform(date = as.POSIXct(date, tz = 'GMT',origin = "1970-01-01"))
  
  #make data frame for species counts in elevation bands per day in band
  bandNewSp <- cbind(day = banddays, #adding 1 day so NA replacements below work
                     mapply(bands, FUN = function(x){
                       sapply(banddays, countSpEl, df = df, el = x)%>%cumsum()
                     })) %>%
    as.data.frame
  
  #Remove points after the end of band sampling
  for (i in 1:length(bands)){
    bandNewSp[which(bandNewSp$day > max(na.omit(filter(df, eband == bands[i])$dayinband))),bands[i]] <- NA
  }
  bandNewSp$day <- bandNewSp$day - 1 #subtracting 1 day to make up for the hack on line 13
  
  #melt dataset for plotting
  bandAccTable <- melt(bandNewSp, id.vars = 'day') %>% rename(c('day' = 'day',
                                                                'band' = 'variable',
                                                                'Species_encountered' = 'value')) %>% 
    filter(!is.na(Species_encountered))
  
  bandAcc <- ggplot(bandAccTable, aes(x = day, y = Species_encountered, group = band, linetype = band))+
    geom_line(linewidth = 1.5)+
    ylab("Cumulative Species Encountered")+
    xlab("Days of Sampling in Elevational Band")+
    labs(title = "B")+
    theme_minimal() +
    scale_x_continuous(expand = c(0,0),
                       breaks = seq(1,max(bandAccTable$day),1),
                       labels = seq(1,max(bandAccTable$day),1),
                       limits = c(1,max(bandAccTable$day)),
                       minor_breaks = NULL) +
    scale_y_continuous(breaks = seq(0,max(bandAccTable$Species_encountered),5),
                       labels = seq(0,max(bandAccTable$Species_encountered),5),
                       position = "right",
                       minor_breaks = NULL,
                       limits = c(0,max(bandAccTable$Species_encountered) + 1),
                       expand = c(0,0)) +
    scale_linetype_manual(values = c("solid","longdash","dotdash","dotted"))
    theme_minimal()
  
  #For accumulation curves by taxon
  mx <- max(newsp$cumulative)
  
  taxonAccTable <- melt(newsp %>% 
                          select(c("cufrog","culizard","cusnake","cumulative","date")) %>%
                          rename(c("frog"="cufrog",
                                   "lizard" = "culizard",
                                   "snake" = "cusnake",)), 
                        id.vars = "date") %>%
    rename(c("taxon" = "variable","Species_encountered" = "value"))
  
  # ggplot
  taxonAcc <- ggplot(taxonAccTable, aes(x = date, y=Species_encountered, group = taxon, color = taxon, linetype = taxon)) +
    geom_line(linewidth = 1.5) +
    labs(title = "A")+
    ylab("Cumulative Species Encountered") +
    xlab("Date of Sampling") +
    scale_x_datetime(breaks = seq(min(newsp$date),max(newsp$date), by = 'day'),
                     minor_breaks = NULL,
                     labels = format(seq(min(newsp$date),max(newsp$date), by = 'day'), format ="%b %d"),
                     expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0,mx,5),
                       minor_breaks = NULL,
                       position = "right",
                       expand = c(0,1)) +
    
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = c("longdash", "dashed", "dotted","solid")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
          #legend.position = c(0,1),
          #legend.justification = c(0,1),
          #legend.background = element_rect(fill = "white", colour = "white"),
          plot.margin = margin(l = 20))
  
  return(list('acc' = newsp, 'taxonacc' = taxonAcc, 'bandacc' = bandAcc))
}

addGroups<- function(df){
  id <- df$'binom' # extract IDs
  #Lists of frog, snake, lizard, & turtle genera
  grp <- list(
  frog = c("Tadpole","Tadpoles","Bufo","Chalcorana","Duttaphrynus","Hylarana","Ingerophrynus","Kaloula","Limnonectes","Occidozyga","Oreophryne","Papurana","Polypedates","Rhacophorus","Tadpole","tadpole"),
  snake = c("Shed","Ahaetulla","Amphiesma","Boiga","Calamaria","Coelognathus","Calamorhabdium","Chrysopelea","Cylindrophis","Dendrelaphis","Elaphe","Enhydris","Hebius","Hypsiscopus","Indotyphlops","Gonyosoma","Oligodon","Ophiophagus","Psammodynastes","Ptyas","Rabdion","Rhabdophis","snake","Tropidolaemus","Typhlops","Xenochrophis","Xenopeltis"),
  lizard = c("Bronchocoela","Cyrtodactylus","Dibamus","Draco","Emoia","Eutropis","Gehyra","Gekko","Hemidactylus","Hemiphyllodactylus","Lamprolepis","Lipinia","Sphenomorphus","Tiliqua","Tytthoscincus","Varanus"),
  turtle = c("Cuora","Leucocephalon")
  )
  
  #find and replace ids with general ID (frog, snake, lizard, turtle)
  grpmatch <- function(x, reference = grp){
    x <- sub(" .*", "", x)
    group <- x
    for (i in names(reference)){
     if(x %in% reference[[i[1]]]){
       group <- i
     }
    }
    return(group)
  }
  
  df$group <-sapply(id, FUN = grpmatch)
  df$matched <- TRUE
  df$matched[which(!(df$group %in% names(grp)))] <- FALSE
  if(!all(df$matched)){cat("could not match JAM #s",df$JAM_Number[!df$matched])}
  return(df)
}

dfSetup <- function(df){
  colnames(df) <- sub(" ","_",colnames(df))
  colnames(df) <- gsub("[()]","",colnames(df))
  df%>% mutate(binom = Genus_species) %>% problemRows()
}

elevationLookup <-function(r){#uses the elevatr package to fill out gaps in elevation data
  require('elevatr')
  print(paste0("Estimating Elevation for JAM",df[r,]$JAM_Number))
  if(is.na(df[r,]$Longitude)){
    print(paste0("No coordinates - skipping JAM",df[r,]$JAM_Number))
    return(NA)
  }
  else{
    return(get_elev_point(data.frame(x = df[r,]$Longitude, y = df[r,]$Latitude),
                          prj = "EPSG:4326", 
                          src = "aws",
                          z = 15)@data$elevation,)
    print(paste0("Got coordinates for JAM",df[r,]$JAM_Number))}
}

elevBands <- function(df,bbs,cf=0){ #df and a vector of cutoff elevations between elevational bands
  #label and set up elevational bands
  absmax <- bbs[length(bbs)] + 10000
  bbs <- c(bbs,absmax)
  
  names(bbs)[1] <- paste0("<",bbs[1],"m")
  
  for(i in 2:length(bbs-1)){
    names(bbs)[i] <- paste0(bbs[i-1], "-", bbs[i], "m")
  }
  
  names(bbs)[length(bbs)] <- paste0(">",bbs[length(bbs)-1],"m")
  
  #get the df ready to write elevational bands to
  df <- df[order(df$Collection_Date),] %>% #reorder df by Collection Date
    mutate(firstocc = rep(FALSE,times = length('binom')),
           firstinband = rep(FALSE,times = length('binom')))
  
  #insert absent elevation data
  
  #insert elevation for species without elevation data
  
  # if there are any specimens with coordinates but no elevation data
  # insert elevations using elevatr
  if(length(which(is.na(df$Elevation) & !is.na(df$Longitude))) !=0){
  
  elevest <- df[which(is.na(df$Elevation) & !is.na(df$Longitude)),c("Longitude","Latitude")] %>%
    as.data.frame 
  colnames(elevest) <- c("x","y")
  elevest <- get_elev_point(elevest, prj = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
                   src = 'aws', 
                   z = 14) #14 (9.6m resolution) is the highest-resolution dataset available for Sulawesi
   
  
  df[which(is.na(df$Elevation) & !is.na(df$Longitude)),"Elevation"] <- elevest$elevation +
    cf #add any correction factor to the elevations
  }
  #assign bands to specimens
  df$eband <- cut(df$Elevation, 
                  breaks = c(0,bbs),
                  labels = c(names(bbs))) %>%
    as.factor()
  
  df$dayinband <- 0
  
  
  
  # Label by Day within elevational band
  for(i in 1:length(bbs)){
  
  map1 <- setNames(c(1:length(unique(df[which(df$eband == names(bbs)[i]),]$Collection_Date))),
                   unique(df[which(df$eband == names(bbs)[i]),]$Collection_Date))
  df[which(df$eband == names(bbs)[i]),]$dayinband <- map1[as.character(df[which(df$eband == names(bbs)[i]),]$Collection_Date)]
  }
  
  ####FIRST OCCURRENCES####
  #Generate first occurrence data (overall)
  df$firstocc[match(unique(df$'binom'),
                    df$'binom')] <- TRUE
  
  #Generate first occurrence data (within band)
  #This marks which specimens represent the first occurrence of a species within each elevational band
  #the syntax of this is 
  #df$firstinband[specimens_in_band[first occurrence of each species within the band]]
  #the match(unique()): 
  #unique() makes a list of all the unique elements
  #match() finds the first occurrence of each unique element within some list
  #I know it's ghastly
  #band 1
  for(i in 1:length(bbs)){
    df$firstinband[which(df$eband == names(bbs[i]))[match(unique(df$binom[which(df$eband == names(bbs[i]))]),
                                                                                             df$binom[which(df$eband == names(bbs[i]))])]] <- TRUE
  }

  
  return(df)
}

getrange <- function(x,y){
  range <- x
  range[which(x!=y)] <- (paste0(x[which(x!=y)],"-",y[which(x!=y)]))
  range[which(range == "Inf--Inf")] <- "*"
  return(range)
}

missingSpecies <- function(df) {
  require('vegan')
  
  diversityEstimates <- rbind(specpool(veganize(df), smallsample = TRUE),
                              mapply(levels(df$eband), 
                                            FUN = function(x){
                                              specpool(veganize(filter(df, eband == x)), 
                                                       smallsample = TRUE)%>% unlist
                                            }) %>%t(),
                              specpool(veganize(filter(df, group == "frog")), smallsample = TRUE),
                              specpool(veganize(filter(df, group == "lizard")), smallsample = TRUE),
                              specpool(veganize(filter(df, group == "snake")), smallsample = TRUE)) %>%
    round(digits = 2) %>%
    as.data.frame()
  
  rownames(diversityEstimates) <- c("Total Mountain",
                                    levels(df$eband),
                                    "Frogs",
                                    "Lizards",
                                    "Snakes")
  
  pubDiversityEstimates <- diversityEstimates[,c("Species","n","chao","chao.se","jack2")] %>%
    mutate ("Missing species" = paste0(round(chao - Species),"/",round(jack2 - Species)))
  colnames(pubDiversityEstimates) <- c("Observed Species","Number of days sampled",
                                       "Chao estimate","Chao standard error",
                                       "Second-order Jackknife estimate", "Estimated Missing Species (Chao/Jackknife)")
  
  pubDiversityEstimates <- pubDiversityEstimates%>% rownames_to_column("Subset") %>% 
    mutate("Specimens Found" = c(dim(df)[1],table(df$eband),table(df$group)), .after = Subset)
  
  return(list("MissingDiversity" = diversityEstimates,
              "PublicationTable" = pubDiversityEstimates))
}

plotElev <- function(df, colors){
  require('tidyverse')
  require('ggplot2')
  require('cowplot')
  require('divDyn')
  
  edf<-df %>% filter(!is.na(Elevation))
  edf$Elevation <- as.integer(edf$Elevation)
  
  divDynGroup<- function(taxon = 'total',df){
    if(taxon != 'total'){
      df<-df%>%filter(group == taxon)
    }
    df$eBin <- df$Elevation %>% round(digits = -1) %>% as.factor
    g <- divDyn(df,tax = 'binom', bin = "eBin")%>%select(c('eBin','divRT'))
    g$eBin <- as.numeric(g$eBin)
    g$group <- taxon
    g
  }
  dd <- lapply(c("total",unique(edf$group)), #data on range-through diversity
               FUN = divDynGroup, df = edf)%>%
    bind_rows()
  
  elevation_order <- sapply(levels(as.factor(edf$binom)), #species list
                            function(x){max(filter(edf, binom == x)$Elevation)}) %>% # get max elevation for a species
    order()
  
  edf$binom <- factor(edf$binom, levels = levels(as.factor(edf$binom))[elevation_order])
  min.max <- data.frame(binom = levels(as.factor(edf$binom)),
                        max.el = sapply(levels(as.factor(edf$binom)),
                                        function(x){max(filter(edf, binom == x)$Elevation)}),
                        min.el = sapply(levels(as.factor(edf$binom)),
                                        function(x){min(filter(edf, binom == x)$Elevation)})) %>%
    addGroups()
  
  #for plotting, it's useful to remember what the cutoff values for the elevational bands were.
  #extracting them as variables here
  
  bbs <- regmatches(levels(edf$eband),(gregexpr("[0-9]+",levels(edf$eband)))) %>% #find all the elevation values
    unlist() %>% #unlist the result
    unique() %>% #remove all the duplicated values
    as.numeric() %>% #make them numbers instead of strings
    sort() #just in case of some weirdness..

  
  #make data frame with just first occurrences per elevational band
  bpdf <- edf[which(edf$firstinband==TRUE),]
  
  #bar plot of species within bands
  bar <- ggplot(data = bpdf, aes(x = eband, group = group, fill = group)) +
    geom_bar(width = 0.5) + #leaving just a bit of whitespace
    scale_fill_manual(values = colors) +
    coord_flip() +
    xlab(NULL)+
    scale_y_continuous(expand = c(0, 1),
                       position = "right") +
    
    ylab("Species Tally") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(angle = 90, hjust = .5),
          plot.margin = unit(c(0,0,0,0), "cm"),
          plot.background = element_rect(fill='transparent', color=NA),)
  
  #scatterplot of all specimen elevations, organized by species  
  scatter <- ggplot() +
    geom_hline(yintercept = bbs, color = "#586A6A", linetype = "dashed") +
    
    geom_point(data = edf, aes(x = binom, y=Elevation, group = group, color = group, fill = group), size = 0.5)+
    geom_linerange(data = min.max, aes(x = binom, ymin = min.el, ymax =  max.el), colour = "#586A6A", linetype = "dotted") +
    geom_point(data = edf, aes(x = binom, y=Elevation, group = group, color = group, fill = group, pch = group), size = 2.5)+
    
  
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = c(17,19,25,4)) +
    ylab("Elevation (m)") +
    xlab("Species") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, face = "italic", vjust = 0.5, hjust=1),
          plot.margin = unit(c(0,0,0,0), "cm"),
          legend.position = "none",
          legend.justification = c(0,1),
          legend.background = element_rect(fill = "white", colour = "white"))
  
  #line graph of range-through diversity
  rtDiv <- ggplot(dd, aes(x = eBin, y = divRT, group = group, color = group, linetype = group)) +
    geom_line(linewidth = 1.5)+ 
    geom_vline(xintercept = bbs, color = "#586A6A", linetype = "dashed") +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = c("longdash", "dashed", "dotted","solid")) +
    theme_minimal()+
    xlab(NULL)+
    
    ylab("Species Richness") +
    theme_minimal() +
    scale_y_continuous(minor_breaks = seq(0, max(dd$divRT), 1))+
    theme(legend.position = "none",
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0,0,0,0.5), "cm"),
          plot.background = element_rect(fill='transparent', color=NA))
  
  facetrt <-plot_grid(rtDiv + 
              theme(axis.text.y = element_text(color = "black")),
            scatter + 
              coord_flip() + 
              scale_x_discrete(limits = rev) + 
              theme(axis.text.x = element_text(face = "plain", angle = 0,vjust = 0.5, hjust=0.5),
                    axis.text.y = element_text (face = "italic")) +
              facet_grid(rows = vars(group), scales = "free", space = "free"),
            nrow = 2, rel_heights = c(1,4), align = "v", axis = "rl")
  
  
  aligned <- align_plots(scatter, bar, align="hv", axis="br")
  combined <- plot_grid(aligned[[1]],aligned[[2]],rel_widths = c(4, 1))
  rangeThrough <- align_plots(scatter, rtDiv, align="h", axis="tblr")
  rangeThrough <- plot_grid(rangeThrough[[1]],rangeThrough[[2]],
            rel_widths = c(length(unique(df$binom)), max(dd$divRT)))
  return(list("combined" = combined, "scatter" = scatter, "bar" = bar, "rangeThrough" = facetrt, "rtDiv" = rtDiv))
}

problemRows <- function(df){
  df$Flag <- NA
  
  df$Flag[c(which(is.na(df$JAM_Number)),
            which(is.na(df$Collection_Date)),
          which(is.na(df$Genus_species)),
          which(is.na(df$Latitude)))] <- TRUE
  return(df)
}

quickScrub <- function(df){
  df <- df %>% dfSetup() %>% addGroups()
  
  df <- filter(df, is.na(Flag) & matched)
  
  filter(df, group != "turtle")
}

speciesTable <- function(df) {
  speciestable <- df%>%group_by(binom)%>% # this df organizes species  and reports
    summarise(Number = n(), #number of specimens
              minMASS  = min(Weight_g, na.rm = TRUE), #minimum mass
              maxMASS  = max(Weight_g, na.rm = TRUE), #maximum mass
              minSVL = min(SVL_mm, na.rm = TRUE), #minimum svl
              maxSVL = max(SVL_mm, na.rm = TRUE), #maximum svl
              minel = min(Elevation, na.rm = TRUE), #minimum elevation
              maxel = max(Elevation, na.rm = TRUE)) %>% #maximum elevation 
    addGroups()
  
  pubTable <- data.frame(group = speciestable$group,
                         Species = speciestable$binom,
                         n = speciestable$Number,
                         "SVL" = getrange(speciestable$minSVL,speciestable$maxSVL),
                         "Mass" = getrange(speciestable$minMASS,speciestable$maxMASS),
                         "Elevation" = getrange(speciestable$minel,speciestable$maxel)) %>%
    arrange(group,Species)
  pubTable$SVL <- str_replace(pubTable$SVL, "-",'–')
  pubTable$Mass <- str_replace(pubTable$Mass, "-",'–')
  pubTable$Elevation <- str_replace(pubTable$Elevation, "-",'–')
  
  return(list("speciesTable" = speciesTable,
              "publicationTable" = pubTable))
  
}

veganize <- function(df, elevation = FALSE){
  
  vdf<-filter(df,!is.na(df$Collection_Date))
  
  spp <- unique(vdf$binom) #species list
  
  if (elevation == FALSE){
    dat <- unique(vdf$Collection_Date) #collection date list
    
    sppcount <- function(x,y) { #count number of species x found on day Y
      dim(filter(vdf,binom == x, Collection_Date == y))[1]
    }
    
    f <- function(x) {sapply(dat,sppcount,x=x)} #counts for each species on day y
    
    
    lapply(spp, f) %>% #counts for each species on each day
      as.data.frame(row.names = format(dat, format ="%b %d"),
                    col.names = spp)
  }
  else{
    eb <- levels(vdf$eband) #collection date list
    
    sppcount <- function(x,y) { #count number of species x found in band y
      dim(filter(vdf,binom == x, eband == y))[1]
    }
    
    f <- function(x) {sapply(eb,sppcount,x=x)} #counts for each species on day y
    
    
    lapply(spp, f) %>% #counts for each species in each band
      as.data.frame(row.names = eb,
                    col.names = spp)
    
  }
  
}
