################################################################################
# 
# ResidFunctions.R
# Version 1.0
# 28/01/2015
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
#
################################################################################

### Nachman
ResidNachman <- function(par, observed, Area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Nachman model
  #
  # Args:
  #   par: dataframe containing parameters C and z of the Nachman model
  #   A: Grain size (km2) to be predicted
  resids <- observed - PredictNachman(par, Area)
  return(resids)
}

### Power Law
ResidPL <- function(par, observed, Area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Power Law model
  #
  # Args:
  #   par: dataframe containing parameters C and z of the Power Law model
  #   A: Grain size (km2) to be predicted
  resids <- observed - PredictPL(par, Area)
  return(resids)
}

### Logistic
ResidLogis <- function(par, observed, Area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Logistic model
  #
  # Args:
  #   par: dataframe containing parameters C and z of the Logistic model
  #   A: Grain size (km2) to be predicted
  resids <- observed - PredictLogis(par, Area)
  return(resids)
}

### Poisson
ResidPoisson <- function(par, observed, Area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Poisson model
  #
  # Args:
  #   par: dataframe containing parameter lambda of the Poisson model
  #   A: Grain size (km2) to be predicted
  resids <- observed - PredictPoisson(par, Area)
  return(resids)
}

### Negative binomial
ResidNB <- function(par, observed, Area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Negative binomial model
  #
  # Args:
  #   par: dataframe containing parameters C and k of the Negative
  #        Binomial model
  #   A: Grain size (km2) to be predicted
  resids <- observed - PredictNB(par, Area)
  return(resids)
}

### Generalised negative binomial model
ResidGNB <- function(par, observed, Area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Generalised negative binomial model
  #
  # Args:
  #   par: dataframe containing parameters C, z and k of the Generalised 
  #        Negative Binomial model
  #   A: Grain size (km2) to be predicted
  resids <- observed - PredictGNB(par, Area)
  return(resids)
}

### Improved negative binomial model
ResidINB <- function(par, observed, Area) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Improved negative binomial model
  #
  # Args:
  #   par: dataframe containing parameters C and b of the Improved 
  #        Negative Binomial model
  #   A: Grain size (km2) to be predicted
  resids <- observed - PredictINB(par, Area)
  return(resids)
}

### Finite negative binomial model
ResidFNB <- function(par, observed, Area, A0) {
  # Calculates residual between observed and expected area of occupancy for 
  # grain size A using the Finite negative binomial model
  #
  # Args:
  #   par: dataframe containing parameters W and k of the Finite
  #        Negative Binomial model
  #   A: Grain size (km2) to be predicted
  #   A0: Total area (km2)
  resids <- observed - PredictFNB(par, Area, A0)
  return(resids)
}

