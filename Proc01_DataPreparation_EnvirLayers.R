#===============================================================
# Proc01 Data preparation - Environmental layers
#===============================================================

rm(list = ls())
Sys.setenv(language="EN")
# start <- Sys.time()

# 1) Working directory ----------------------------------------------------


setwd("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/AnnualModeling_SOC")



# 2) Libraries ------------------------------------------------------------


pckg <- c('terra',     
          'magrittr',
          'magrittr',
          'devtools',
          'raster',
          'parallel',
          'rassta'
)

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
lapply(pckg,usePackage)


# 3) Static predictors ----------------------------------------------------

coast <- vect("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/CoastLine/Kyst.shp")
#data.frame(coast)
#plot(coast)

# 3.1) Geomorphology ------------------------------------------------------

all <- list.files("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Geomorphology/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
)

all <- lapply(all,function(x){
  rast(x)
})


geomorphology <- rast(all) %>% crop(coast,mask=T)

# 3.2) Geology ------------------------------------------------------------

all <- list.files("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Geology/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
)

all <- lapply(all,function(x){
  rast(x)
})

geology <- rast(all) %>% crop(coast,mask=T)

dummygeo1 <- dummies(geology[[1]],preval=1)
dummygeo2 <- dummies(geology[[2]],preval=1)
dummygeo1 <- dummygeo1[[-nlyr(dummygeo1)]]
dummygeo2 <- dummygeo2[[-nlyr(dummygeo2)]]

# 3.3) Soil ---------------------------------------------------------------

# 3.3.1) Clay layers for different depth intervals ------------------------

all <- list.files("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Soil/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
)

all <- lapply(all,function(x){
  rast(x)
})

soil <- rast(all) %>% crop(coast,mask=T)


# 3.4) Coordinates --------------------------------------------------------

cov <- rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif")
names(cov) <- readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds")

coord <- as.data.frame(soil[[1]], xy=TRUE)
x <- data.frame(coord$x,coord$y,coord$x) 
x <- rast(x, type="xyz")
crs(x)  <- "epsg:25832"
x <- x %>% crop(coast,mask=T) %>% resample(soil[[1]])
names(x) <- "CoordX"
x2 <- x^2
names(x2) <- "CoordX2"
x3 <- x^3
names(x3) <- "CoordX3"

y <- data.frame(coord$x,coord$y,coord$y)
y <- rast(y, type="xyz")   
crs(y)  <- "epsg:25832"
y <- y %>% crop(coast,mask=T) %>% resample(soil[[1]])
names(y) <- "CoordY"
y2 <- y^2
names(y2) <- "CoordY2"
y3 <- y^3
names(y3) <- "CoordY3"


# 4) Static predictors raster stack ---------------------------------------
StatPreds <- c(x,x2,x3,y,y2,y3,soil,dummygeo1,dummygeo2,geomorphology)
names(StatPreds)

terra::writeRaster(StatPreds,
                   "O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds05082022.tif",
                   datatype="FLT2S",
                   overwrite=T)

saveRDS(names(StatPreds),
        "O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds")


# 5) Dynamic predictors ---------------------------------------------------


# 5.1) Organisms ----------------------------------------------------------

indicesL5 <- function(bandstack){
  ndvi <- (bandstack[[4]]-bandstack[[3]])/(bandstack[[4]]+bandstack[[3]])
  
  kndvi <- tanh(((bandstack[[4]]-bandstack[[3]])/(bandstack[[4]]+bandstack[[3]]))^2)
  
  evi <- ((bandstack[[4]]-bandstack[[3]])/(bandstack[[4]] + 6*bandstack[[3]] - 7.5*bandstack[[1]] + 1))*2.5
  
  bsi <- ((bandstack[[6]] + bandstack[[3]]) - (bandstack[[3]]+bandstack[[1]]))/((bandstack[[6]] + bandstack[[3]]) + (bandstack[[3]]+bandstack[[1]]))
  
  savi <- (bandstack[[4]]-bandstack[[3]])/((bandstack[[4]]+bandstack[[3]] + 0.5)*1.5)
  
  str <- ((1-bandstack[[6]])^2)/(2*bandstack[[6]])
  
  brightness <- bandstack[[1]]*0.2043 + bandstack[[2]]*0.4158 + bandstack[[3]]*0.5524 + 
    bandstack[[4]]*0.5741 + bandstack[[5]]*0.3124 + bandstack[[6]]*0.2303
  
  greenness <- bandstack[[1]]*(-0.1603) + bandstack[[2]]*(-0.2819) + bandstack[[3]]*(-0.4934) + 
    bandstack[[4]]*0.7940 + bandstack[[5]]*0.0002 + bandstack[[6]]*(-0.1446) 
  
  wetness <- bandstack[[1]]*0.0315 + bandstack[[2]]*0.2021 + bandstack[[3]]*0.3102 + 
    bandstack[[4]]*0.1594 + bandstack[[5]]*0.6806 + bandstack[[6]]*(-0.6109)
  
  msavi <- (2*bandstack[[4]] + 1- sqrt((2*bandstack[[4]] + 1)^2 - 8 * (bandstack[[4]] - bandstack[[3]])))/2
 
   indices <- c(ndvi,kndvi,evi,savi,msavi,bsi,str,brightness,greenness,wetness)
 
    return(indices)
}

# 5.1.1) Period 1 1984-1986 -----------------------------------------------

paths <- dir("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/Year1985_1986",full.names = T)

all <- list()
for (i in 1:length(paths)) {
  all[[i]] <-  list.files(paths[i],
                          pattern = "tiff$|tif$",
                          full.names = TRUE
  ) %>% rast() %>% crop(coast,mask=T)

}

# a) 1985 -----------------------------------------------------------------
ind1985med <- app(all[[1]][[1:6]], fun=indicesL5, cores =15)
ind1985p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =15)
ind1985p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =15)

# b) 1986 -----------------------------------------------------------------
ind1986med <- app(all[[2]][[1:6]], fun=indicesL5, cores =15)
ind1986p10 <- app(all[[2]][[7:12]], fun=indicesL5, cores =15)
ind1986p90 <- app(all[[2]][[13:18]], fun=indicesL5, cores =15)

# c) Period 1 layers (organisms) ------------------------------------------

organismsP1 <- c(ind1985med,ind1985p10,ind1985p90,
                 ind1986med,ind1986p10,ind1986p90)


# 5.1.2) Period 2 1996-1997 -----------------------------------------------

paths <- dir("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/Year1996_1997",full.names = T)

all <- list()
for (i in 1:length(paths)) {
  all[[i]] <-  list.files(paths[i],
                          pattern = "tiff$|tif$",
                          full.names = TRUE
  ) %>% rast() %>% crop(coast,mask=T)

}

# a) 1996 -----------------------------------------------------------------
ind1996med <- app(all[[1]][[1:6]], fun=indicesL5, cores =15)
ind1996p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =15)
ind1996p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =15)

# b) 1997 -----------------------------------------------------------------
ind1997med <- app(all[[2]][[1:6]], fun=indicesL5, cores =15)
ind1997p10 <- app(all[[2]][[7:12]], fun=indicesL5, cores =15)
ind1997p90 <- app(all[[2]][[13:18]], fun=indicesL5, cores =15)

# c) Period 2 layers (organisms) ------------------------------------------

organismsP2 <- c(ind1996med,ind1996p10,ind1996p90,
                 ind1997med,ind1997p10,ind1997p90)

# 5.1.3) Period 3 2008-2009 -----------------------------------------------

paths <- dir("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/Year2008_2009",full.names = T)

all <- list()
for (i in 1:length(paths)) {
  all[[i]] <-  list.files(paths[i],
                          pattern = "tiff$|tif$",
                          full.names = TRUE
  ) %>% rast() %>% crop(coast,mask=T)

}

# a) 2008 -----------------------------------------------------------------
ind2008med <- app(all[[1]][[1:6]], fun=indicesL5, cores =15)
ind2008p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =15)
ind2008p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =15)

# b) 2009 -----------------------------------------------------------------
ind2009med <- app(all[[2]][[1:6]], fun=indicesL5, cores =15)
ind2009p10 <- app(all[[2]][[7:12]], fun=indicesL5, cores =15)
ind2009p90 <- app(all[[2]][[13:18]], fun=indicesL5, cores =15)

# c) Period 3 layers (organisms) ------------------------------------------

organismsP3 <- c(ind2008med,ind2008p10,ind2008p90,
                 ind2009med,ind2009p10,ind2009p90)


#  5.1.4) Period 4 2018-2019 ----------------------------------------------

paths <- dir("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/Year2018_2019",full.names = T)

all <- list()
for (i in 1:length(paths)) {
  all[[i]] <-  list.files(paths[i],
                          pattern = "tiff$|tif$",
                          full.names = TRUE
  ) %>% rast() %>% crop(coast,mask=T)

}

# a) 2018 -----------------------------------------------------------------
ind2018med <- app(all[[1]][[1:6]], fun=indicesL5, cores =15)
ind2018p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =15)
ind2018p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =15)

# b) 2019 -----------------------------------------------------------------
ind2019med <- app(all[[2]][[1:6]], fun=indicesL5, cores =15)
ind2019p10 <- app(all[[2]][[7:12]], fun=indicesL5, cores =15)
ind2019p90 <- app(all[[2]][[13:18]], fun=indicesL5, cores =15)


# c) Land use information taken from field block map 2018 -------------------------

crops <- rast("O:/Tech_AGRO/Jord/Sebastian/Fields_2011-2020/CropDataDK2011_2020/CropData2011_2020_30m.tif")
names(crops)
crops1819 <- crops[[c(8,9)]] %>% crop(coast,mask=T) %>% resample(all[[1]],method="near")

requirements <- read.table(file = 'O:/Tech_AGRO/Jord/Sebastian/Fields_2011-2020/CropDataDK2011_2020/crop_drain_req.csv',
                           sep = ';',
                           header = TRUE
)

reclasser <- as.matrix(requirements[, c(1, 6)])

# Check if the levels of each raster layer within the raster stack 
# has a corresponding code in the reclassification table
crops18<-ratify(raster(crops1819[[1]]))
(catg <- levels(crops18)[[1]])

catg$ID[catg$ID%nin%reclasser[,1]]


crops19<-ratify(raster(crops1819[[2]]))
(catg <- levels(crops19)[[1]])
catg$ID[catg$ID%nin%reclasser[,1]]
check <- data.frame(Kode=catg$ID,PermGrass_RotGrass=NA)

# Updating reclassification table
reclasser <- as.matrix(rbind(reclasser,check))


beginCluster(detectCores()-2)
crops_reclass <- clusterR(stack(crops18,crops19), reclassify,
                        args = list(rcl = reclasser))

endCluster()

ref <- raster("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif")
dummy18 <- layerize(crops_reclass[[1]])
names(dummy18) <- c("PermGrass18","GrassRot18","Crops18","BareSoil18")
dummy18 <- raster::resample(dummy18,ref,method="ngb")

dummy19 <- layerize(crops_reclass[[2]])
names(dummy19) <- c("PermGrass19","GrassRot19","Crops19")
dummy19 <- raster::resample(dummy19,ref,method="ngb")

writeRaster(dummy18,"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/FieldBlock18.tif",
            datatype='INT2U',overwrite=T)
saveRDS(names(dummy18),"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesFieldBlock18.rds")

writeRaster(dummy19,"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/FieldBlock19.tif",
            datatype='INT2U',overwrite=T)
saveRDS(names(dummy19),"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesFieldBlock19.rds")


# d) Land use information taken from ESRI LULC maps -----------------------
esri18 <- rast("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Organisms/LANDUSE/ESRILULC2018.tif") %>% 
  terra::resample(geology[[1]],method="near")

esri19 <- rast("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Organisms/LANDUSE/ESRILULC2019.tif") %>% 
  terra::resample(geology[[1]],method="near")

esri18d <- segregate(esri18)
esri19d <- segregate(esri19)

esri18d <- esri18d[[c(2,4,9)]]
names(esri18d) <- c("Forest18","Crop18","Grass18")
esri19d <- esri19d[[c(2,4,8)]]
names(esri19d) <- c("Forest19","Crop19","Grass19")



# d) Period 4 layers (organisms) ------------------------------------------

organismsP4 <- c(ind2018med,ind2018p10,ind2018p90,
                 ind2019med,ind2019p10,ind2019p90,
                 esri18d,esri19d
                 )


# 5.2) Climate ------------------------------------------------------------

# 5.2.1) Period 1 layers (climate) ----------------------------------------

all <- list.files("O:/Tech_AGRO/Jord/Sebastian/BIOLAYERS_CHELSA1985_2021_30M/monitoringgridlayers/Period1/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
) %>% rast() %>% crop(coast,mask=T)

climateP1 <- all

# 5.2.2) Period 2 layers (climate) ----------------------------------------

all <- list.files("O:/Tech_AGRO/Jord/Sebastian/BIOLAYERS_CHELSA1985_2021_30M/monitoringgridlayers/Period2/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
) %>% rast() %>% crop(coast,mask=T)

climateP2 <- all

# 5.2.3) Period 3 layers (climate) ----------------------------------------

all <- list.files("O:/Tech_AGRO/Jord/Sebastian/BIOLAYERS_CHELSA1985_2021_30M/monitoringgridlayers/Period3/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
) %>% rast() %>% crop(coast,mask=T)

climateP3 <- all

# 5.2.4) Period 4 layers (climate) ----------------------------------------

all <- list.files("O:/Tech_AGRO/Jord/Sebastian/BIOLAYERS_CHELSA1985_2021_30M/monitoringgridlayers/Period4/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
) %>% rast() %>% crop(coast,mask=T)

climateP4 <- all


# 6) Dynamic raster stacks ------------------------------------------------

# 6.1) Period 1 -----------------------------------------------------------
start <- Sys.time()
organismsP1 <- resample(organismsP1,StatPreds)
climateP1 <- resample(climateP1,StatPreds)
covP1 <- c(organismsP1,climateP1)
terra::writeRaster(covP1,"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP1985_1986.tif",datatype="FLT4S",overwrite=T)
names(covP1)
names(covP1) <- c(paste0(c("ndvi","kndvi","evi",
                           "savi","msavi","bsi",
                           "str","brightness",
                           "greenness","wetness") %>% rep(6) %>% as.factor(),
                                            gl(3,10,labels=c("Med","P10","P90")),
                                            "_",
                                            gl(2,30,labels=c("1985","1986"))),
                  names(covP1)[61:98])
                  
names(covP1)
saveRDS(names(covP1),"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP1985_1986.rds")


# # 6.2) Period 2 -----------------------------------------------------------
#

organismsP2 <- resample(organismsP2,StatPreds)
climateP2 <- resample(climateP2,StatPreds)
covP2 <- c(organismsP2,climateP2)
terra::writeRaster(covP2,"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP1996_1997.tif",datatype="FLT4S",overwrite=T)

names(covP2) <- c(paste0(c("ndvi","kndvi","evi",
                           "savi","msavi","bsi",
                           "str","brightness",
                           "greenness","wetness") %>% rep(6) %>% as.factor(),
                         gl(3,10,labels=c("Med","P10","P90")),
                         "_",
                         gl(2,30,labels=c("1996","1997"))),
                  names(covP2)[61:98])
saveRDS(names(covP2),"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP1996_1997.rds")

# 6.3) Period 3 -----------------------------------------------------------

organismsP3 <- resample(organismsP3,StatPreds)
climateP3 <- resample(climateP3,StatPreds)
covP3 <- c(organismsP3,climateP3)
terra::writeRaster(covP3,"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP2008_2009.tif",datatype="FLT4S",overwrite=T)

names(covP3) <- c(paste0(c("ndvi","kndvi","evi",
                           "savi","msavi","bsi",
                           "str","brightness",
                           "greenness","wetness") %>% rep(6) %>% as.factor(),
                         gl(3,10,labels=c("Med","P10","P90")),
                         "_",
                         gl(2,30,labels=c("2008","2009"))),
                  names(covP3)[61:98])

saveRDS(names(covP3),"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP2008_2009.rds")

# 6.4) Period 4 -----------------------------------------------------------

organismsP4 <- resample(organismsP4,StatPreds)
climateP4 <- resample(climateP4,StatPreds)
covP4 <- c(organismsP4,climateP4)
terra::writeRaster(covP4,"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP2018_2019_.tif",datatype="FLT4S",overwrite=T)

names(covP4) <- c(paste0(c("ndvi","kndvi","evi",
                           "savi","msavi","bsi",
                           "str","brightness",
                           "greenness","wetness") %>% rep(6) %>% as.factor(),
                         gl(3,10,labels=c("Med","P10","P90")),
                         "_",
                         gl(2,30,labels=c("2018","2019"))),
                  names(covP4)[61:79])
saveRDS(names(covP4),"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP2018_2019.rds")
Sys.time()-start