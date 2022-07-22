#============================================================================
# Proc02 Data preparation - Point data set --------------------------------
#============================================================================
rm(list = ls())
Sys.setenv(language="EN")


# 1) Working directory ----------------------------------------------------


setwd("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/AnnualModeling_SOC")


# 2) Libraries ------------------------------------------------------------


pckg <- c('DescTools',     
          'magrittr',
          'tidyr',
          'readr',
          'dplyr',
          'terra'
)

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
lapply(pckg,usePackage)

# 3) Predictors loading ---------------------------------------------------

covP1 <- c(rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP1985_1986.tif"))
names(covP1) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP1985_1986.rds"))
                  
covP2 <- c(rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP1996_1997.tif"))
names(covP2) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP1996_1997.rds"))

covP3 <- c(rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP2008_2009.tif"))
names(covP3) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP2008_2009.rds"))

covP4 <- c(rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP2018_2019.tif"))
names(covP4) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP2018_2019.rds"))

# 4) Data loading ---------------------------------------------------------

data <- read_delim("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/Points/DataJoined.csv",
                   delim = ";",col_types = "ffddfciiddiciddddid")

databd1 <- data %>% 
  filter(Depth==1 & !is.na(OrgC) & !is.na(BDwhole)) %>% 
  dplyr::select(PointID,BDwhole)

databd2 <- data %>% 
  filter(Depth==2 & !is.na(OrgC) & !is.na(BDwhole)) %>% 
  dplyr::select(PointID,BDwhole)

# 4.1) Depth 00-25 cm -----------------------------------------------------

data1 <- data %>% 
  filter(Depth==1 & !is.na(OrgC)) %>% 
  right_join(databd1,by="PointID",keep=F) %>% 
  mutate(SOC=OrgC*BDwhole.y*25) %>% 
  dplyr::select(PointID,X,Y,Year,SOC) %>% 
  pivot_wider(names_from = Year, values_from = SOC, names_prefix="SOC_")

# 4.2) Depth 25-50 cm -----------------------------------------------------

data2 <- data %>% 
  filter(Depth==2 & !is.na(OrgC)) %>% 
  right_join(databd2,by="PointID",keep=F) %>% 
  mutate(SOC=OrgC*BDwhole.y*25) %>% 
  dplyr::select(PointID,X,Y,Year,SOC) %>% 
  pivot_wider(names_from = Year, values_from = SOC, names_prefix="SOC_")

# 5) Data extraction ------------------------------------------------------

# 5.1) SOC 00-25 cm -------------------------------------------------------

data_sp1 <- vect(data1, geom=c("X", "Y"),crs="epsg:25832")# %>%  project("epsg:25832")#Data frame as spatial points data frame
data_sp1

dataP1D1 <- cbind(data1[,4],terra::extract(covP1,data_sp1))
dataP1D1$ID <- NULL
dataP1D1 <- dataP1D1 %>% na.omit
write_csv(dataP1D1,"RegMat_P1D1.csv")

dataP2D1 <- cbind(data1[,5],terra::extract(covP2,data_sp1))
dataP2D1$ID <- NULL
dataP2D1 <- dataP2D1 %>% na.omit
write_csv(dataP2D1,"RegMat_P2D1.csv")

dataP3D1 <- cbind(data1[,6],terra::extract(covP3,data_sp1))
dataP3D1$ID <- NULL
dataP3D1 <- dataP3D1 %>% na.omit
write_csv(dataP3D1,"RegMat_P3D1.csv")

dataP4D1 <- cbind(data1[,7],terra::extract(covP4,data_sp1))
dataP4D1$ID <- NULL
dataP4D1 <- dataP4D1 %>% na.omit
write_csv(dataP4D1,"RegMat_P4D1.csv")

# 5.2) SOC 25-50 cm -------------------------------------------------------

data_sp2 <- vect(data2, geom=c("X", "Y"),crs="epsg:25832")# %>%  project("epsg:25832")#Data frame as spatial points data frame
data_sp2

dataP1D2 <- cbind(data2[,4],terra::extract(covP1,data_sp2))
dataP1D2$ID <- NULL
dataP1D2 <- dataP1D2 %>% na.omit
write_csv(dataP1D2,"RegMat_P1D2.csv")

dataP2D2 <- cbind(data2[,5],terra::extract(covP2,data_sp2))
dataP2D2$ID <- NULL
dataP2D2 <- dataP2D2 %>% na.omit
write_csv(dataP2D2,"RegMat_P2D2.csv")

dataP3D2 <- cbind(data2[,6],terra::extract(covP3,data_sp2))
dataP3D2$ID <- NULL
dataP3D2 <- dataP3D2 %>% na.omit
write_csv(dataP3D2,"RegMat_P3D2.csv")

dataP4D2 <- cbind(data2[,7],terra::extract(covP4,data_sp2))
dataP4D2$ID <- NULL
dataP4D2 <- dataP4D2 %>% na.omit
write_csv(dataP4D2,"RegMat_P4D2.csv")


