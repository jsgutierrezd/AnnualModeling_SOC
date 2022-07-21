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


set.seed(0)
r <- rast(nrows=10, ncols=10)
values(r) <- sample(3, ncell(r), replace=TRUE)
is.factor(r)
cls <- c("forest", "water", "urban")
# make the raster start at zero
x <- r - 1
levels(x) <- cls
names(x) <- "land cover"
is.factor(x)
x


# 3) Predictors loading ---------------------------------------------------

covP1 <- rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/CovP1984_1986.tif")
names(covP1) <- readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesCovP1984_1986.rds")
levels(covP1[[5]]) <- paste0("geology",1:11)
levels(covP1[[6]]) <- paste0("georeg",1:10)
plot(covP1[[5]])

covP2 <- rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/CovP1994_1997.tif")
names(covP2) <- readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesCovP1994_1997.rds")
levels(covP2[[5]]) <- paste0("geology",1:11)
levels(covP2[[6]]) <- paste0("georeg",1:10)

covP3 <- rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/CovP2006_2009.tif")
names(covP3) <- readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesCovP2006_2009.rds")
levels(covP3[[5]]) <- paste0("geology",1:11)
levels(covP3[[6]]) <- paste0("georeg",1:10)

covP4 <- rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/CovP2016_2019.tif")
names(covP4) <- readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesCovP2016_2019.rds")
levels(covP4[[5]]) <- paste0("geology",1:11)
levels(covP4[[6]]) <- paste0("georeg",1:10)

# 4) Data loading ---------------------------------------------------------

data <- read_delim("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/Points/DataJoined.csv",
                   delim = ";",col_types = "ffddfciiddiciddddid")

databd <- data %>% 
  filter(Depth==1 & !is.na(OrgC) & !is.na(BDwhole)) %>% 
  dplyr::select(PointID,BDwhole)

# 4.1) Depth 00-25 cm -----------------------------------------------------

data1 <- data %>% 
  filter(Depth==1 & !is.na(OrgC)) %>% 
  right_join(databd,by="PointID",keep=F) %>% 
  mutate(SOC=OrgC*BDwhole.y*25) %>% 
  dplyr::select(PointID,X,Y,Year,SOC) %>% 
  pivot_wider(names_from = Year, values_from = SOC, names_prefix="SOC_")

# 4.2) Depth 25-50 cm -----------------------------------------------------

data2 <- data %>% 
  filter(Depth==2 & !is.na(OrgC)) %>% 
  right_join(databd,by="PointID",keep=F) %>% 
  mutate(SOC=OrgC*BDwhole.y*25) %>% 
  dplyr::select(PointID,X,Y,Year,SOC) %>% 
  pivot_wider(names_from = Year, values_from = SOC, names_prefix="SOC_")

# 5) Data extraction ------------------------------------------------------

# 5.1) SOC 00-25 cm -------------------------------------------------------

data_sp <- vect(data1, geom=c("X", "Y"),crs="epsg:25832")# %>%  project("epsg:25832")#Data frame as spatial points data frame
data_sp

dataP1D1 <- cbind(data1[,4],terra::extract(covP1,data_sp))
dataP1D1$ID <- NULL
dataP1D1 <- dataP1D1 %>% na.omit
write_csv(dataP1D1,"RegMat_P1D1.csv")

dataP2D1 <- cbind(data1[,5],terra::extract(covP2,data_sp))
dataP2D1$ID <- NULL
dataP2D1 <- dataP2D1 %>% na.omit
write_csv(dataP2D1,"RegMat_P2D1.csv")

dataP3D1 <- cbind(data1[,6],terra::extract(covP3,data_sp))
dataP3D1$ID <- NULL
dataP3D1 <- dataP3D1 %>% na.omit
write_csv(dataP3D1,"RegMat_P3D1.csv")

dataP4D1 <- cbind(data1[,7],terra::extract(covP4,data_sp))
dataP4D1$ID <- NULL
dataP4D1 <- dataP4D1 %>% na.omit
write_csv(dataP4D1,"RegMat_P4D1.csv")

# 5.2) SOC 25-50 cm -------------------------------------------------------

dataP1D2 <- cbind(data2[,4],terra::extract(covP1,data_sp))
dataP1D2$ID <- NULL
dataP1D2 <- dataP1D2 %>% na.omit
write_csv(dataP1D2,"RegMat_P1D2.csv")

dataP2D2 <- cbind(data2[,5],terra::extract(covP2,data_sp))
dataP2D2$ID <- NULL
dataP2D2 <- dataP2D2 %>% na.omit
write_csv(dataP2D2,"RegMat_P2D2.csv")

dataP3D2 <- cbind(data2[,6],terra::extract(covP3,data_sp))
dataP3D2$ID <- NULL
dataP3D2 <- dataP3D2 %>% na.omit
write_csv(dataP3D2,"RegMat_P3D2.csv")

dataP4D2 <- cbind(data2[,7],terra::extract(covP4,data_sp))
dataP4D2$ID <- NULL
dataP4D2 <- dataP4D2 %>% na.omit
write_csv(dataP4D2,"RegMat_P4D2.csv")


