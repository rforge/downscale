#' Model area of occupancy against grain size for downscaling
#' 
#' @param AOO Observed area of occupancy
#' @param Area Grain size
#' @param Model Downscaling model (see Details)

################################################################################
# 
# downscale.R
# Version 1.1
# 30/01/2015
#
# Updates:
#   30/01/2015: Thomas model added 
#
# Model area of occupancy against grain size for downscaling. 
#
# Args:
#   Area: Vector of grain sizes for obsvered area of occupancies
#   Observed: Vector of observed area of occupancies
#   Model = "FNB"
#
# Returns:
#   optim.pars: list of parameters estimated from optimisation procedure
#
################################################################################

downscale <- function(occupancy, area, model, extent){
  input.data <- DataInput(occupancy = occupancy, area = area, extent = extent)
  model <- model
  
  if ((model == "Nachman") | (model == "PL") | (model == "Logis") | 
        (model == "Poisson") | (model == "NB") | (model == "GNB") | 
        (model == "INB")){
    optim.pars <- suppressWarnings(OptimiseParameters(area =
                                      input.data[!is.na(input.data[, "Occ"]),
                                                   "Cell.area"], 
                                      Observed = 
                                        input.data[!is.na(input.data[, "Occ"]),
                                                    "Occ"],
                                      model = model))
  }

  if (model == "FNB") { 
    optim.pars <- suppressWarnings(
      OptimiseParametersFNB(area =
                              input.data[!is.na(input.data[,"Occ"]),
                                         "Cell.area"], 
                            observed = 
                              input.data[!is.na(input.data[,"Occ"]),
                                         "Occ"],
                            extent = extent,
                            model = model))
  }

  if (model == "Thomas") { 
    optim.pars <- suppressWarnings(
      OptimiseParametersThomas(area =
                                 input.data[!is.na(input.data[,"Occ"]),
                                            "Cell.area"], 
                               observed = 
                                 input.data[!is.na(input.data[,"Occ"]),
                                            "Occ"],
                               extent = extent,
                               model = model,
                               tolerance = 1e-6))
  }
  output <- list("model" = model, "pars" = unlist(optim.pars))
  return(output)
}
