#===============================================================
# Proc01 Data preparation - Environmental layers
#===============================================================
rm(list = ls())


# 1) Working directory ----------------------------------------------------


setwd("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/AnnualModeling_SOC")



# 2) Libraries ------------------------------------------------------------


pckg <- c('terra',     
          'magrittr',
          'magrittr',
          'devtools',
          'raster',
          'parallel'
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



# 4) Static predictors raster stack ---------------------------------------

StatPreds <- c(soil,geology,geomorphology)
names(StatPreds)


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
  indices <- c(ndvi,evi,savi,bsi,str,brightness,greenness,wetness)
  return(indices)
}

# 5.1.1) Period 1 1984-1986 -----------------------------------------------

paths <- dir("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/Year1984_1986",full.names = T)

all <- list()
for (i in 1:length(paths)) {
  all[[i]] <-  list.files(paths[i],
                          pattern = "tiff$|tif$",
                          full.names = TRUE
  ) %>% rast() %>% crop(coast,mask=T)
  
}


# a) 1984 ----------------------------------------------------------------- 
ind1984med <- app(all[[1]][[1:6]], fun=indicesL5, cores =5)
ind1984p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =5)
ind1984p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =5)

# b) 1985 -----------------------------------------------------------------
ind1985med <- app(all[[2]][[1:6]], fun=indicesL5, cores =5)
ind1985p10 <- app(all[[2]][[7:12]], fun=indicesL5, cores =5)
ind1985p90 <- app(all[[2]][[13:18]], fun=indicesL5, cores =5)

# c) 1986 -----------------------------------------------------------------
ind1986med <- app(all[[3]][[1:6]], fun=indicesL5, cores =5)
ind1986p10 <- app(all[[3]][[7:12]], fun=indicesL5, cores =5)
ind1986p90 <- app(all[[3]][[13:18]], fun=indicesL5, cores =5)


# d) 1984-1986 ------------------------------------------------------------

all <- list.files("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/Aggregatedperperiod/Period1/",
                          pattern = "tiff$|tif$",
                          full.names = TRUE
  ) %>% rast() %>% crop(coast,mask=T)
  

indP1med <- app(all[[1:6]], fun=indicesL5, cores =5)
indP1p10 <- app(all[[7:12]], fun=indicesL5, cores =5)
indP1p90 <- app(all[[13:18]], fun=indicesL5, cores =5)


# e) Period 1 layers (organisms) ------------------------------------------

organismsP1 <- c(ind1984med,ind1984p10,ind1994p90,ind1985med , ind1985p10,ind1985p90 ,ind1986med , ind1986p10, ind1986p90, indP1med, indP1p10, indP1p90)


# 5.1.2) Period 2 1994-1997 -----------------------------------------------

paths <- dir("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/Year1994_1997",full.names = T)

all <- list()
for (i in 1:length(paths)) {
  all[[i]] <-  list.files(paths[i],
                          pattern = "tiff$|tif$",
                          full.names = TRUE
  ) %>% rast() %>% crop(coast,mask=T)
  
}



# a) 1994 -----------------------------------------------------------------
ind1994med <- app(all[[1]][[1:6]], fun=indicesL5, cores =5)
ind1994p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =5)
ind1994p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =5)

# b) 1995 -----------------------------------------------------------------
ind1995med <- app(all[[1]][[1:6]], fun=indicesL5, cores =5)
ind1995p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =5)
ind1995p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =5)

# c) 1996 -----------------------------------------------------------------
ind1996med <- app(all[[1]][[1:6]], fun=indicesL5, cores =5)
ind1996p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =5)
ind1996p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =5)

# d) 1997 -----------------------------------------------------------------
ind1997med <- app(all[[1]][[1:6]], fun=indicesL5, cores =5)
ind1997p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =5)
ind1997p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =5)


# e) 1994-1997 ------------------------------------------------------------

all <- list.files("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/Aggregatedperperiod/Period2/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
) %>% rast() %>% crop(coast,mask=T)

indP2med <- app(all[[1:6]], fun=indicesL5, cores =5)
indP2p10 <- app(all[[7:12]], fun=indicesL5, cores =5)
indP2p90 <- app(all[[13:18]], fun=indicesL5, cores =5)


# f) Period 2 layers (organisms) ------------------------------------------

organismsP2 <- c(ind1994med,ind1994p10,ind1994p90,ind1995med , ind1995p10,ind1995p90 ,ind1996med , ind1996p10, ind1996p90, ind1997med, ind1997p10, ind1997p90,
                 indP2med,
                 indP2p10,
                 indP2p90)

# 5.1.3) Period 3 2006-2009 -----------------------------------------------

paths <- dir("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/Year2006_2009",full.names = T)

all <- list()
for (i in 1:length(paths)) {
  all[[i]] <-  list.files(paths[i],
                          pattern = "tiff$|tif$",
                          full.names = TRUE
  ) %>% rast() %>% crop(coast,mask=T)
  
}

# a) 2006 -----------------------------------------------------------------
ind2006med <- app(all[[1]][[1:6]], fun=indicesL5, cores =5)
ind2006p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =5)
ind2006p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =5)

# b) 2007 -----------------------------------------------------------------
ind2007med <- app(all[[1]][[1:6]], fun=indicesL5, cores =5)
ind2007p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =5)
ind2007p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =5)

# c) 2008 -----------------------------------------------------------------
ind2008med <- app(all[[1]][[1:6]], fun=indicesL5, cores =5)
ind2008p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =5)
ind2008p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =5)

# d) 2009 -----------------------------------------------------------------
ind2009med <- app(all[[1]][[1:6]], fun=indicesL5, cores =5)
ind2009p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =5)
ind2009p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =5)


# e) 2006-2009 ------------------------------------------------------------

all <- list.files("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/Aggregatedperperiod/Period3/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
) %>% rast() %>% crop(coast,mask=T)

indP3med <- app(all[[1:6]], fun=indicesL5, cores =5)
indP3p10 <- app(all[[7:12]], fun=indicesL5, cores =5)
indP3p90 <- app(all[[13:18]], fun=indicesL5, cores =5)


# f) Period 3 layers (organisms) ------------------------------------------

organismsP3 <- c(ind2006med,ind2006p10,ind2006p90,ind2007med , ind2007p10,ind2007p90 ,ind2008med , ind2008p10, ind2008p90, ind2009med, ind2009p10, ind2009p90,
                 indP3med,
                 indP3p10,
                 indP3p90)


#  5.1.4) Period 4 2016-2019 ----------------------------------------------

paths <- dir("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/Year2016_2019",full.names = T)

all <- list()
for (i in 1:length(paths)) {
  all[[i]] <-  list.files(paths[i],
                          pattern = "tiff$|tif$",
                          full.names = TRUE
  ) %>% rast() %>% crop(coast,mask=T)
  
}

# a) 2016 -----------------------------------------------------------------
ind2016med <- app(all[[1]][[1:6]], fun=indicesL5, cores =5)
ind2016p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =5)
ind2016p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =5)

# b) 2017 -----------------------------------------------------------------
ind2017med <- app(all[[1]][[1:6]], fun=indicesL5, cores =5)
ind2017p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =5)
ind2017p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =5)

# c) 2018 -----------------------------------------------------------------
ind2018med <- app(all[[1]][[1:6]], fun=indicesL5, cores =5)
ind2018p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =5)
ind2018p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =5)

# d) 2019 -----------------------------------------------------------------
ind2019med <- app(all[[1]][[1:6]], fun=indicesL5, cores =5)
ind2019p10 <- app(all[[1]][[7:12]], fun=indicesL5, cores =5)
ind2019p90 <- app(all[[1]][[13:18]], fun=indicesL5, cores =5)


# e) 2016 - 2019 ----------------------------------------------------------

all <- list.files("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/Aggregatedperperiod/Period4/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
) %>% rast() %>% crop(coast,mask=T)

indP4med <- app(all[[1:6]], fun=indicesL5, cores =5)
indP4p10 <- app(all[[7:12]], fun=indicesL5, cores =5)
indP4p90 <- app(all[[13:18]], fun=indicesL5, cores =5)


# f) Period 4 layers (organisms) ------------------------------------------

organismsP4 <- c(ind2016med,ind2016p10,ind2016p90,ind2017med , ind2017p10,ind2017p90 ,ind2018med , ind2018p10, ind2018p90, ind2019med, ind2019p10, ind2019p90,
                 indP4med,
                 indP4p10,
                 indP4p90)



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



# 6) Static and dynamic raster stack --------------------------------------


# 6.1) Period 1 -----------------------------------------------------------

covP1 <- c(StatPreds,organismsP1,climateP1)
terra::writeRaster(covP1,"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/CovP1984_1986.tif",datatype="FLT4S")
saveRDS(names(covP1),"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesCovP1984_1986.rds")

# 6.2) Period 2 -----------------------------------------------------------

covP2 <- c(StatPreds,organismsP2,climateP2)
terra::writeRaster(covP2,"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/CovP1994_1997.tif",datatype="FLT4S")
saveRDS(names(covP2),"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesCovP1994_1997.rds")

# 6.3) Period 3 -----------------------------------------------------------

covP3 <- c(StatPreds,organismsP3,climateP3)
terra::writeRaster(covP3,"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/CovP2006_2009.tif",datatype="FLT4S")
saveRDS(names(covP3),"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesCovP2006_2009.rds")

# 6.4) Period 4 -----------------------------------------------------------

covP4 <- c(StatPreds,organismsP4,climateP4)
terra::writeRaster(covP4,"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/CovP2016_2019.tif",datatype="FLT4S")
saveRDS(names(covP4),"O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesCovP2016_2019.rds")