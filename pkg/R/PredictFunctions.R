################################################################################
# 
# PredictFunctions.R
# Version 1.1
# 28/01/2015
#
# Functions for downscaled prediction of area of occupancy for a given 
# grain size (A) given the parameters for that model. Model parameters may
# be the outputs from model optimisation to coarse-scale data
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
# requires(Rmpfr)
#   # Note: The Finite Negative Binomial model requires multiple precision 
#     floating-point numbers in package Rmpfr for calculation (see description)
#
################################################################################

### Nachman
PredictNachman <- function(par, Area) {
  # Predicts area of occupancy for grain size A using the Nachman model
  #
  # Args:
  #   par: dataframe containing parameters C and z of the Nachman model
  #   A: Grain size (km2) to be predicted
  AOO <- 1 - exp(-par$C * Area ^ par$z)
  return(AOO)
}

### Power Law
PredictPL <- function(par, Area) {
  # Predicts area of occupancy for grain size A using the Power law model
  #
  # Args:
  #   par: dataframe containing parameters C and z of the Power law model
  #   A: Grain size (km2) to be predicted
  AOO <- par$C * Area^par$z
  return(AOO)
}

### Logistic
PredictLogis <- function(par, Area) {
  # Predicts area of occupancy for grain size A using the Logistic model
  #
  # Args:
  #   par: dataframe containing parameters C and z of the Logistic model
  #   A: Grain size (km2) to be predicted
  AOO <- (par$C * (Area^par$z)) / (1 + (par$C * (Area^par$z)))
  return(AOO)
}

### Poisson
PredictPoisson <- function(par, Area) {
  # Predicts area of occupancy for grain size A using the Poisson model
  #
  # Args:
  #   par: dataframe containing parameter lambda of the Poisson model
  #   A: Grain size (km2) to be predicted
  AOO <- 1 - (exp(-par$lambda * Area))
  return(AOO)
}

### Negative binomial model
PredictNB <- function(par, Area) {
  # Predicts area of occupancy for grain size A using the Negative 
  # Binomial model
  #
  # Args:
  #   par: dataframe containing parameters C and k of the Negative
  #        Binomial model
  #   A: Grain size (km2) to be predicted
  AOO <- 1 - (1 + (par$C * Area) / par$k)^-par$k
  return(AOO)
}

### Generalised negative binomial model
PredictGNB <- function(par, Area) {
  # Predicts area of occupancy for grain size A using the Generalised Negative 
  # Binomial model
  #
  # Args:
  #   par: dataframe containing parameters C, z and k of the Generalised 
  #        Negative Binomial model
  #   A: Grain size (km2) to be predicted
  AOO <- 1 - (1 + (par$C * Area^par$z) / par$k)^-par$k
  return(AOO)
}

### Improved negative binomial model
PredictINB <- function(par, Area) {
  # Predicts area of occupancy for grain size A using the Improved Negative 
  # Binomial model
  #
  # Args:
  #   par: dataframe containing parameters C and b of the Improved 
  #        Negative Binomial model
  #   A: Grain size (km2) to be predicted
  AOO <- 1 - ((par$C* Area^(par$b - 1))^
                ((par$r * Area) /(1 - par$C * Area^(par$b - 1))))
  return(AOO)
}

### Finite negative binomial model
PredictFNB <- function(par, Area, A0){
  # Predicts area of occupancy for grain size A using the Finite Negative 
  # Binomial model. The function multiplies many gamma functions and so
  # integers may become larger than possible in R. Therefore  we use multiple
  # precision floatinf point numbers (the 'mpfr' function in package 'Rmpfr')
  # is used to make calculations possible.
  #
  # Args:
  #   par: dataframe containing parameters W and k of the Finite
  #        Negative Binomial model
  #   A: Grain size (km2) to be predicted
  #   A0: Total area (km2)
  gamma1 <- par$W + ((A0 * par$k) / Area) - par$k
  gamma2 <- (A0 * par$k) / Area
  gamma3 <- par$W + ((A0 * par$k) / Area)
  gamma4 <- ((A0 * par$k) / Area) - par$k
  AOO <- as.numeric(1 - (
    (gamma(Rmpfr::mpfr(gamma1, 128)) * gamma(Rmpfr::mpfr(gamma2, 128))) /
      (gamma(Rmpfr::mpfr(gamma3, 128)) * gamma(Rmpfr::mpfr(gamma4, 128)))))
  return(AOO)
}
