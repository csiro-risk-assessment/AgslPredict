## ----include = FALSE--------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  echo = TRUE,
  comment = "#>"
)

## ----ra, eval = FALSE-------------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict", "01-RelativeAbundanceData/")

## ----gc,eval=FALSE----------------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict", "covariates_spatial/active.csv")

## ----mcvar,eval=FALSE-------------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "meteo_for_relative_abundance/meteo_download_RelativeAbundance.R")

## ----dat1, eval = FALSE-----------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict", "covariates_spatial")

## ----rac,eval=FALSE---------------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict", "01-RelativeAbundanceData/RelativeAbundanceDataCovariates.R")

## ----raa, eval = FALSE------------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "02-RelativeAbundanceAnalysis/relative_abundance_analysis_1m.R")

## ----rab, eval = FALSE------------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "02-RelativeAbundanceAnalysis/multinomialGLMfit1m.rds")

## ----adat, eval = FALSE-----------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "03-AbundanceData/AbundancePrep.R")

## ----apdat, eval = FALSE----------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "meteo_for_abundance/meteo_abundance_download.R")

## ----spac, eval = FALSE-----------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "03-AbundanceData/AbundanceCovarPrep.R")

## ----aanalysis, eval = FALSE------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "04-AbundanceAnalysis/agsl_negbinom.R")

## ----dat3s, eval = FALSE----------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "meteo_prediction_raw/download_from_power.R")

## ----pred1, eval = FALSE----------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "05-Prediction/prediction_tot_relabund.R")

## ----pred2, eval = FALSE----------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "05-Prediction/pXPrediction_function.R")

## ----pred2a, eval = FALSE---------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "05-Prediction/tempMask.R")

## ----pred3, eval = FALSE----------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "05-Prediction/prediction_plots.R")

## ----pred4, eval = FALSE----------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "05-Prediction/analysisSpec.R")

## ----metprep,eval=FALSE-----------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict",
#              "covar_plots/calc_metavg.R")

## ----metplots,eval=FALSE----------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict", "covar_plots/met_plots.R")

## ----spatplots,eval=FALSE---------------------------------------------------------------------------------------------------------------------------
#  system.file(package = "AgslPredict", "covar_plots/spatial_plots.R")

