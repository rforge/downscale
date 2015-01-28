################################################################################
# 
# StartingParams.R
# Version 1.0
# 28/01/2015
#
# Dataframes containing the starting parameters for the optimisation procedure
# fitting model parameters to observed coarse-scale occupancy data:
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

### Nachman model
ParamsNachman <- data.frame("C" = 0.01, "z" = 0.01)

### Power Law model
ParamsPL <- data.frame("C" = 0.01, "z" = 0.01)

### Logistic model
ParamsLogis <- data.frame("C" = 0.01, "z" = 0.01)

### Poisson model
ParamsPoisson <- data.frame("lambda" = 0.001)

### Negative binomial model
ParamsNB <- data.frame("C" = 0.01, "k" = 0.01)

### Generalised negative binomial model
ParamsGNB <- data.frame("C" = 0.1, "z" = 0.001, "k" = 1)

### Improved negative binomial model
ParamsINB <- data.frame("C" = 0.01, "r" = 1e-5, "b" = 1)

### Finite negative binomial model
ParamsFNB <- data.frame("W" = 10, "k" = 10)
