################################################################################
# 
# ResidFunctions.R
# Version 1.1
# 30/01/2015
#
# Updates:
#   30/01/2015: Thomas model added
#
# Functions for calculating the residuals between observed and expected area of
# occupancy used in the optimisation procedure for parameter fitting to coarse
# scale data.
#
# This file contains the 8 functions that are simple geometric extrapolations of
# the area-occupancy relationship at coarse grain sizes:
#   Nachman   Nachman model
#   PL        Power Law model
#   Logis     Logistic model
#   Poisson   Poisson model
#   NB        Negative binomial model
#   GNB       Generalised negative binomial model
#   INB       Improved negative binomial model
#   FNB       Finite negative binomial model
#   Thomas    Thomas model
#
################################################################################

### Nachman
ResidNachman <- function(par, observed, area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Nachman model
  #
  # Args:
  #   par: dataframe containing parameters C and z of the Nachman model
  #   area: Grain size (km2) to be predicted
  resids <- observed - PredictNachman(par, area)
  return(resids)
}

### Power Law
ResidPL <- function(par, observed, area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Power Law model
  #
  # Args:
  #   par: dataframe containing parameters C and z of the Power Law model
  #   area: Grain size (km2) to be predicted
  resids <- observed - PredictPL(par, area)
  return(resids)
}

### Logistic
ResidLogis <- function(par, observed, area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Logistic model
  #
  # Args:
  #   par: dataframe containing parameters C and z of the Logistic model
  #   area: Grain size (km2) to be predicted
  resids <- observed - PredictLogis(par, area)
  return(resids)
}

### Poisson
ResidPoisson <- function(par, observed, area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Poisson model
  #
  # Args:
  #   par: dataframe containing parameter lambda of the Poisson model
  #   area: Grain size (km2) to be predicted
  resids <- observed - PredictPoisson(par, area)
  return(resids)
}

### Negative binomial
ResidNB <- function(par, observed, area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Negative binomial model
  #
  # Args:
  #   par: dataframe containing parameters C and k of the Negative
  #        Binomial model
  #   area: Grain size (km2) to be predicted
  resids <- observed - PredictNB(par, area)
  return(resids)
}

### Generalised negative binomial model
ResidGNB <- function(par, observed, area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Generalised negative binomial model
  #
  # Args:
  #   par: dataframe containing parameters C, z and k of the Generalised 
  #        Negative Binomial model
  #   area: Grain size (km2) to be predicted
  resids <- observed - PredictGNB(par, area)
  return(resids)
}

### Improved negative binomial model
ResidINB <- function(par, observed, area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Improved negative binomial model
  #
  # Args:
  #   par: dataframe containing parameters C and b of the Improved 
  #        Negative Binomial model
  #   area: Grain size (km2) to be predicted
  resids <- observed - PredictINB(par, area)
  return(resids)
}

### Finite negative binomial model
ResidFNB <- function(par, observed, area, extent) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Finite negative binomial model
  #
  # Args:
  #   par: dataframe containing parameters W and k of the Finite
  #        Negative Binomial model
  #   area: Grain size (km2) to be predicted
  #   extent: Total area (km2)
  resids <- observed - PredictFNB(par, area, extent)
  return(resids)
}

### Thomas model
ResidThomas <- function(par, observed, area, extent, tolerance = 1e-6) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Thomas model
  #
  # Args:
  #   par: dataframe containing parameters rho, mu and sigman the Thomas model
  #   area: Grain size (km2) to be predicted
  #   extent: Total area (km2)
  #   tolerance: Tolerance to be given in intergration - the smaller the number
  #   the greater the accuracy but longer the processing time
  resids <- observed - 
    PredictThomas(par, area, extent, tolerance = tolerance)
  return(resids)
}
