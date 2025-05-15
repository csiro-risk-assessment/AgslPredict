## ITN and IRS plots for years 2002, 2008, 2014 and 2020

ITN2002 <- rast("//fs1-cbr.nexus.csiro.au/{d61-pmb-publications}/work/AgslDistribution/AgslPredict_instantiation/AgslPredict_instantiated_202504/vector_intervention/ITN_2002.tif")
ITN2008 <- rast("//fs1-cbr.nexus.csiro.au/{d61-pmb-publications}/work/AgslDistribution/AgslPredict_instantiation/AgslPredict_instantiated_202504/vector_intervention/ITN_2008.tif")
ITN2014 <- rast("//fs1-cbr.nexus.csiro.au/{d61-pmb-publications}/work/AgslDistribution/AgslPredict_instantiation/AgslPredict_instantiated_202504/vector_intervention/ITN_2014.tif")
ITN2020 <- rast("//fs1-cbr.nexus.csiro.au/{d61-pmb-publications}/work/AgslDistribution/AgslPredict_instantiation/AgslPredict_instantiated_202504/vector_intervention/ITN_2020.tif")

# projection to WGS84
ITN2002 <- project(ITN2002, "epsg:4326")
ITN2008 <- project(ITN2008, "epsg:4326")
ITN2014 <- project(ITN2014, "epsg:4326")
ITN2020 <- project(ITN2020, "epsg:4326")

ITN.list <- rast(ITN2002, nlyrs = 4, names = paste0("ITN_",c("2002","2008","2014","2020"),"_use_mean"))
ITN.list[[1]] <- ITN2002
ITN.list[[2]] <- ITN2008
ITN.list[[3]] <- ITN2014
ITN.list[[4]] <- ITN2020

IRS2002 <- rast("//fs1-cbr.nexus.csiro.au/{d61-pmb-publications}/work/AgslDistribution/AgslPredict_instantiation/AgslPredict_instantiated_202504/vector_intervention/IRS_2002.tif")
IRS2008 <- rast("//fs1-cbr.nexus.csiro.au/{d61-pmb-publications}/work/AgslDistribution/AgslPredict_instantiation/AgslPredict_instantiated_202504/vector_intervention/IRS_2008.tif")
IRS2014 <- rast("//fs1-cbr.nexus.csiro.au/{d61-pmb-publications}/work/AgslDistribution/AgslPredict_instantiation/AgslPredict_instantiated_202504/vector_intervention/IRS_2014.tif")
IRS2020 <- rast("//fs1-cbr.nexus.csiro.au/{d61-pmb-publications}/work/AgslDistribution/AgslPredict_instantiation/AgslPredict_instantiated_202504/vector_intervention/IRS_2020.tif")

# projection to WGS84
IRS2002 <- project(IRS2002, "epsg:4326")
IRS2008 <- project(IRS2008, "epsg:4326")
IRS2014 <- project(IRS2014, "epsg:4326")
IRS2020 <- project(IRS2020, "epsg:4326")

IRS.list <- rast(IRS2002, nlyrs = 4, names = paste0("IRS_",c("2002","2008","2014","2020"),"_use_mean"))
IRS.list[[1]] <- IRS2002
IRS.list[[2]] <- IRS2008
IRS.list[[3]] <- IRS2014
IRS.list[[4]] <- IRS2020

# plots ------------------------------------------------------------------------

library(rasterVis)

RColorBrewer::display.brewer.all(n=7, type="seq", select=NULL, exact.n=TRUE,
                                                     colorblindFriendly=TRUE)
RColorBrewer::display.brewer.all(n=7, type="div", select=NULL, exact.n=TRUE,
                                 colorblindFriendly=TRUE)

## annual ----------------------------------------------------------------------

cols <- rev(c(RColorBrewer::brewer.pal(9, "GnBu")[9:1]))
ann.theme.itn <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)
cols <- rev(c(RColorBrewer::brewer.pal(9, "BuPu")[9:1]))
ann.theme.irs <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)

library(marmap)

png("Yearly_ITN.png", 2.5*480, 2.5*480*(0.8), pointsize = 72)
levelplot(ITN.list[[1:4]],
          main = expression("ITN use mean"),
          names.attr = paste0("Year ", c("2002","2008","2014","2020")),
          layout = c(2, 2),
          par.settings = ann.theme.itn
)
dev.off()




png("Yearly_IRS.png", 2.5*480, 2.5*480*(0.8), pointsize = 72)
levelplot(IRS.list[[1:4]],
          main = expression("IRS use mean"),
          names.attr = paste0("Year ", c("2002","2008","2014","2020")),
          layout = c(2, 2),
          par.settings = ann.theme.irs
)
dev.off()

