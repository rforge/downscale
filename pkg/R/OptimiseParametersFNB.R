################################################################################
# 
# OptimiseParametersFNB.R
# Version 1.0
# 28/01/2015
#
# Optimisation procedure of finding parameters that best fit observed data for
# the Finite Negative Binomial model
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

OptimiseParametersFNB <- function(Area, Observed, A0, model = "FNB") {
  # Retrive residual function, downscaling function and starting parameters
  # for model of choice
  resid.fun <- getFunction(paste("Resid", model, sep = ""))
  pred.fun <- getFunction(paste("Predict", model, sep = ""))  
  starting.pars <- get(paste("Params", model, sep = ""))
  
  # Optimisation procedure
  optimisation <- minpack.lm::nls.lm(par = starting.pars,
                                     fn = resid.fun,
                                     A = Area,
                                     A0 = A0,
                                     observed = log(Observed),
                                     lower = c("W" = -Inf, "k" = 0),
                                     upper = c("W" = Inf, 
                                               "k" = min(Area, na.rm = TRUE) * 
                                                 100),
                                     control = minpack.lm::nls.lm.control(
                                       maxiter = 1000))
  optim.pars <- as.list(coef(optimisation))
  return(optim.pars)
}
