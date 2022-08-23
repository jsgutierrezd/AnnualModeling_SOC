#============================================================================
# Proc03d_DataModelingP4D1 ------------------------------------------------
#============================================================================
rm(list = ls())
Sys.setenv(language="EN")


# 1) Working directory ----------------------------------------------------


setwd("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/AnnualModeling_SOC")


# 2) Libraries ------------------------------------------------------------

pckg <- c('raster',     
          'terra',
          'ggridges',
          'rasterVis',
          'RColorBrewer',
          'viridis',
          'hrbrthemes'
)

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
lapply(pckg,usePackage)


# 3) Data loading ---------------------------------------------------------

P1D1 <- raster("Outputs/Layers/P1D1/ModelP1D1QrfMedian_080822.tif")*0.1
names(P1D1) <- "SOCstock1986"
P2D1 <- raster("Outputs/Layers/P2D1/ModelP2D1QrfMedian_010822.tif")*0.1
names(P2D1) <- "SOCstock1997"
P3D1 <- raster("Outputs/Layers/P3D1/ModelP3D1QrfMedian_010822.tif")*0.1
names(P3D1) <- "SOCstock2009"
P4D1 <- raster("Outputs/Layers/P4D1/ModelP4D1QrfMedian_010822.tif")*0.1
names(P4D1) <- "SOCstock2019"
P.D1 <- stack(P1D1,P2D1,P3D1,P4D1)
myTheme <- rasterTheme(region = brewer.pal(9,"YlOrRd"))
x11()
levelplot(P.D1,par.settings = myTheme,main="SOC stock (Kg/m2)")


# 4) Ridgeline plot -------------------------------------------------------

ggplot(lincoln_weather, aes(x = `Mean Temperature [F]`, y = `Month`, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Temp. [F]", option = "C") +
  labs(title = 'Temperatures in Lincoln NE in 2016') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )

