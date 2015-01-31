################################################################################
# 
# predict.downscale.R
# Version 1.1
# 30/01/2015
#
# Updates:
#   30/01/2015: Thomas model added
#
# Prediction at other scales
# alias predict
# need class (downscale)
#
################################################################################

predict.downscale <- function(downscale.mod, areas, extent, tolerance = 1e-6) {
  areas.stand <- areas #^ 2 / extent
  params <- as.list(downscale.mod$pars)
  model <- downscale.mod$model
  predict.function <- getFunction(paste("Predict", model, sep = ""))
  
  if ((model == "Nachman") | (model == "PL") | (model == "Logis") | 
        (model == "Poisson") | (model == "NB") | (model == "GNB") | 
        (model == "INB")){    
    AOO <- exp(predict.function(par = params, area = areas.stand))
  }
  
  if (model == "FNB") {
    AOO <- exp(predict.function(par = params, 
                                area = areas.stand, 
                                extent = extent))
  }
  
  if (model == "Thomas") {
    AOO <- exp(predict.function(par = params,
                                tolerance = tolerance,
                                area = areas.stand, 
                                extent = extent))
  }
  expected <- data.frame("Cell.area" = areas, "Occupancy" = AOO)
  return(expected)
}