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

covP1 <- rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/CovP1984_1986.tif")
names(COVP1) <- readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesCovP1984_1986.rds")

covP2 <- rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/CovP1994_1997.tif")
names(COVP2) <- readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesCovP1994_1997.rds")

covP3 <- rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/CovP2006_2009.tif")
names(COVP3) <- readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesCovP2006_2009.rds")

covP4 <- rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/CovP2016_2019.tif")
names(COVP4) <- readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesCovP2016_2019.rds")

# 4) Data loading ---------------------------------------------------------

# 4.1) Depth 00-25 cm -----------------------------------------------------

data1 <- read_delim("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/Points/PointsOC00_25cm.csv",
                    delim = ",",col_types = "fdddddddddddddddddddddddddd") %>% 
  select(PointID,X,Y,SOC_1986,SOC_1997,SOC_2009,SOC_2019)

head(data1)

# 4.2) Depth 25-50 cm -----------------------------------------------------

data2 <- read_delim("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/Points/PointsOC25_50cm.csv",
                    delim = ",",col_types = "fdddddddddddddddddddddddddd") %>% 
  select(Point,x,y,SOC1986,SOC1997,SOC2009,SOC2019)

names(data2) <- c("PointID",
                  "X",
                  "Y",
                  "SOC_1986",
                  "SOC_1997",
                  "SOC_2009",
                  "SOC_2019")
head(data2)

