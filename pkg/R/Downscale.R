#' Model area of occupancy against grain size for downscaling
#'
#' @param AOO Observed area of occupancy
#' @param Area Grain size
#' @param Model Downscaling model (see Details)

################################################################################
# 
# downscale.R
# Version 1.0
# 28/01/2015
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

downscale <- function(AOO, Area, model, A0 = Null){
  model <- model
  if (model != "FNB") { 
    optim.pars <- OptimiseParameters(Area = Area, 
                                     Observed = AOO,
                                     model = model)
  }
  if (model == "FNB") { 
    optim.pars <- OptimiseParametersFNB(Area = Area, 
                                        Observed = AOO,
                                        A0 = A0,
                                        model = model)
  }
  return(optim.pars)
}
