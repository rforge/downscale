#' Model area of occupancy against grain size for downscaling
#' 
#' @aliases downscale
#' 
#' @description The function fits the log observed area of occupancy against 
#'   grain size for coarse-scale data (typically atlas data) for nine commonly
#'   used downscaling models (see Azaele et al. 2012 and Barwell et al. 2014).
#'   The parameters of the fitted models may then be used to estimate the area
#'   of occupancy at finer grain sizes than the observed data using 
#'   \code{\link{predict.downscale}}.
#'   
#' @param occupancy vector of observed area of occupancies in squared units
#'   (e.g. km^{2}).
#' @param area vector of grain sizes (in same units as occupancy).
#' @param model selected downscaling model, chosen from one of \code{"Nachman"},
#'   \code{"PL"}, \code{"Logis"}, \code{"Poisson"}, \code{"NB"}, \code{"GNB"},
#'   \code{"INB"}, \code{"FNB"}, \code{"Thomas"}. See \code{Details} below for
#'   model descriptions.
#' @param extent total area in same units as occupancy (only required for 
#'   \code{FNB} and \code{Thomas} models).
#' @param tolerance only applicable for the \code{Thomas} model. The tolerance 
#'   used during integration in the Thomas model during optimisation of 
#'   parameters. Lower numbers allow for greater accuracy but require longer 
#'   processing times (default = \code{1e-6}).
#'   
#' @details Nine downscaling models are available. \code{area} is the grain size
#'   and \code{extent} the total area in the same units:
#'  \tabular{llll}{
#'    \code{"Nachman"} \tab  \tab Nachman model \tab 
#'      \emph{log(1 - exp(-C * area ^ z))}\cr
#'    \code{"PL"} \tab  \tab  Power law model \tab 
#'      \emph{log(C * area ^ z)}\cr
#'    \code{"Logis"} \tab  \tab  Logistic model \tab 
#'      \emph{log((C * (area ^ z)) / (1 + (C * (area ^ z))))}\cr
#'    \code{"Poisson"} \tab  \tab  Poisson model \tab 
#'      \emph{log(1 - (exp(-lambda * area)))}\cr
#'    \code{"NB"} \tab  \tab  Negative binomial model \tab 
#'      \emph{log(1 - (1 + (C * area) / k) ^ -k)}\cr
#'    \code{"GNB"} \tab  \tab  Generalised negative binomial model \tab 
#'      \emph{log(1 - (1 + (C * area ^ z) / k) ^ -k)}\cr
#'    \code{"INB"} \tab  \tab  Improved negative binomial model \tab 
#'      \emph{log(1 - ((C * area ^ (b - 1)) ^ ((r * area) / 
#'        (1 - C * area ^ (b - 1)))))}\cr
#'    \code{"FNB"} \tab  \tab  Finite negative binomial model \tab 
#'      \emph{log(1 - 
#'        ((gamma(W + ((extent * k) / area) - k) * gamma(extent * k) / area) /
#'        (gamma(W + ((extent * k) / area)) * gamma(((extent * k) / area) - k)))
#'        }\cr
#'    \code{"Thomas"} \tab  \tab  Thomas model \tab see below\cr
#'    }
#'  
#'  The finite negative binomial model (\code{"FNB"}) incorporates several gamma
#'  functions. This may result in integers larger than is possible to store in 
#'  R. Therefore multiple precision floating point numbers (\code{\link{mpfr}}
#'  function in package \pkg{Rmpfr}) are used to make calculations possible.
#'  
#'  For the Thomas model LOTS OF BORING STUFF HERE
#'   
#' @return \code{downscale} returns an object of class "downscale" containing
#'   three objects:
#'   \item{model}{Downscaling model selected.} 
#'   \item{pars}{Estimated parameters for the downscaling model.}
#'    \item{observed}{Data frame containing two columns: 
#'    \tabular{lll}{
#'      \code{Cell.area} \tab  \tab Grain sizes for which occupancy have been
#'        observed\cr 
#'      \code{Occupancy} \tab  \tab Observed area of occupancy for each grain 
#'        size\cr
#'      }
#'    }
#'   
#' @author Charles Marsh <\email{charliem2003@@gmail.com}> with input from
#'   Louise Barwell.
#'   
#' @references Azaele, S., Cornell, S.J., & Kunin, W.E. (2012). Downscaling
#'   species occupancy from coarse spatial scales. \emph{Ecological
#'   Applications} 22, 1004-1014.
#' @references Barwell, L.J., Azaele, S., Kunin, W.E., & Isaac, N.J.B. (2014).
#'   Can coarse-grain patterns in insect atlas data predict local occupancy?
#'   \emph{Diversity and Distributions} 20, 895-907.
#'   
#' @seealso The function output may be used as the input for
#'   \code{\link{predict.downscale}} for extrapolating downscaling functions to
#'   smaller grain sizes using the estimated parameters from the downscale
#'   output.
#'   
#' @example R/Examples/Examples.R
#'
#' @export
#' 
#' @exportClass downscale

################################################################################
# 
# downscale.R
# Version 1.3
# 02/02/2015
#
# Updates:
#   03/02/2015: output defined as class 'downscale'
#   03/02/2015: observed data included in output
#   02/02/2015: help file updated
#   02/02/2015: error check for model name and extent
#   30/01/2015: Thomas model added 
#
# Model area of occupancy against grain size for downscaling. 
#
# Args:
#   occupancy: Vector of observed area of occupancies
#   area: Vector of grain sizes for observed area of occupancies (eg km2)
#   model: function to use
#   extent: total area (same units as area)- required only for FNB and Thomas
#           models
#   tolerance: tolerance for integration of Thomas model. Lower numbers allow 
#              for greater accuracy but require longer processing times
#
# Returns:
#   list of three objects of class 'downscale'
#     $model    Downscaling model used
#     $pars     List of parameters estimated from optimisation procedure
#     $observed Dataframe of grain sizes and observed occupancies used for 
#               modelling
#
################################################################################

downscale <- function(occupancy, area, model, extent, tolerance = 1e-6) {
  input.data <- DataInput(occupancy = occupancy, area = area, extent = extent)
  model <- model
  # error checking - model name is correct
  if (model %in% c("Nachman", "PL", "Logis", "Poisson", "NB", "GNB", "INB",
                   "FNB", "Thomas") == FALSE) {
    stop("Model name invalid", call. = FALSE)
  }
  
  # error checking - extent is bigger than largest grain size
  if ((model == "FNB") | (model == "Thomas")){
    if (max(area, na.rm = TRUE) > extent) {
      stop("Largest grain size is greater than total extent", call. = FALSE)
    }
  }
  
  if ((model == "Nachman") | (model == "PL") | (model == "Logis") | 
        (model == "Poisson") | (model == "NB") | (model == "GNB") | 
        (model == "INB")){
    optim.pars <- suppressWarnings(OptimiseParameters(area =
                                      input.data[!is.na(input.data[, "Occ"]),
                                                   "Cell.area"], 
                                      observed = 
                                        input.data[!is.na(input.data[, "Occ"]),
                                                    "Occ"],
                                      model = model))
  }

  if (model == "Logis") { 
    optim.pars <- suppressWarnings(
      OptimiseParametersLogis(area =
                              input.data[!is.na(input.data[,"Occ"]),
                                         "Cell.area"], 
                            observed = 
                              input.data[!is.na(input.data[,"Occ"]),
                                         "Occ"],
                            model = model))
  }
  
  if (model == "GNB") { 
    optim.pars <- suppressWarnings(
      OptimiseParametersGNB(area =
                                input.data[!is.na(input.data[,"Occ"]),
                                           "Cell.area"], 
                              observed = 
                                input.data[!is.na(input.data[,"Occ"]),
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
                               tolerance = tolerance))
  }
  observed <- data.frame("Cell.area" = input.data[,"Cell.area"],
                         "Occupancy" = input.data[,"Occ"])
  output <- list("model" = model,
                 "pars" = unlist(optim.pars),
                 "observed" = observed)
  class(output) <- "downscale"
  return(output)
}
