#' @title Predict occupancy at fine grain sizes
#' @name predict.downscale
#' @aliases predict
#' 
#' @description Predict area of occupancy at fine grain sizes using parameters 
#'   in a downscale object estimated from coarse grain sizes using 
#'   \code{\link{downscale}}. There is also a simple \code{\link{plot}}
#'   function.
#'
#' @usage predict(object, newdata, extent, tolerance = 1e-6, plot = FALSE)
#' 
#' @param object a fitted object of class \code{'downscale'}.
#' @param newdata vector of grain sizes (in squared units e.g. km^{2}) for which
#'   area of occupancy will be predicted.
#' @param extent total area in same units as newdata (only required for 
#'   \code{FNB} and \code{Thomas} models).
#' @param tolerance only applicable for the \code{Thomas} model. The tolerance 
#'   used during integration in the Thomas model during optimisation of 
#'   parameters. Lower numbers allow for greater accuracy but require longer 
#'   processing times (default = \code{1e-6}).
#' @param plot if \code{plot = TRUE} (default = \code{FALSE}) plots observed and
#'   predicted occupancies against grain size on a log-log plot.
#'   
#' @details The function takes the parameters for a downscaling model estimated 
#'   through \code{\link{downscale}} and uses the model to predict area of 
#'   occupancy at finer grain sizes. See \code{\link{downscale}} for details on
#'   the downscaling models and their parameterisation. Plotting can be called
#'   directly from \code{\link{predict.downscale}} or from
#'   \code{\link{plot.downscale}}.
#'   
#' @return Returns an object of class \code{'predict.downscale'} with three
#'   objects:
#'    \item{model}{Downscaling model used.} 
#'    \item{predicted}{Data frame containing two columns: 
#'    \tabular{lll}{
#'      \code{Cell.area} \tab  \tab Grain sizes for which occupancy have been
#'        estimated\cr 
#'      \code{Occupancy} \tab  \tab Predicted area of occupancy for each grain 
#'        size\cr
#'      }
#'    }
#'    \item{observed}{Data frame containing two columns: 
#'    \tabular{lll}{
#'      \code{Cell.area} \tab  \tab Grain sizes for which occupancy have been
#'        observed\cr 
#'      \code{Occupancy} \tab  \tab Observed area of occupancy for each grain 
#'        size\cr
#'      }
#'    }
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
#' @seealso See \code{\link{downscale}} for estimating parameters of a 
#'   downscaling function from observed occupancies at coarse grain sizes
#'   
#' @examples some examples
#' 
#' @export predict.downscale


################################################################################
# 
# predict.downscale.R
# Version 1.2
# 30/01/2015
#
# Updates:
#   03/02/2015: plot function added
#   03/02/2015: output defined as class 'downscale'
#   03/02/2015: observed data included in output
#   02/02/2015: Help file added to
#   30/01/2015: Thomas model added
#
# Predict area of occupancy at fine grain sizes using parameters in a downscale
# object estimated from coarse grain sizes using downscale
#
# Args:
#   object: model output of class 'downscale' (created from function downscale)
#   newdata: vector of grain sizes for model prediction
#   extent: total area (same units as newdata)- required only for FNB and Thomas
#           models
#   tolerance: tolerance for integration of Thomas model. Lower numbers allow
#              for greater accuracy but require longer processing times
#   plot: if TRUE plots observed and predicted occupancies against grain size on
#         a log-log plot
# 
# Returns:
#   list of three objects of class 'predict.downscale'
#     $model      Downscaling model used
#     $predicted  Dataframe of grain sizes and predicted occupancies
#     $observed   Dataframe of grain sizes and observed occupancies used for 
#                 modelling
#
################################################################################

predict.downscale <- function(object, 
                              newdata, 
                              extent, 
                              tolerance = 1e-6, 
                              plot = FALSE) {
  # error checking
  if (class(object) != "downscale"){
    stop("Input data not of class 'downscale'")
  }
  params <- as.list(object$pars)
  model <-object$model
  predict.function <- getFunction(paste("Predict", model, sep = ""))
  
  if ((model == "Nachman") | (model == "PL") | (model == "Logis") | 
        (model == "Poisson") | (model == "NB") | (model == "GNB") | 
        (model == "INB")){    
    AOO <- exp(predict.function(par = params, area = newdata))
  }
  
  if (model == "FNB") {
    AOO <- exp(predict.function(par = params, 
                                area = newdata, 
                                extent = extent))
  }
  
  if (model == "Thomas") {
    AOO <- exp(predict.function(par = params,
                                tolerance = tolerance,
                                area = newdata, 
                                extent = extent))
  }
  expected <- data.frame("Cell.area" = newdata, "Occupancy" = AOO)
  output <- list("model" = model,
                 "predicted" = expected,
                 "observed" = object$observed)
  class(output) <- "predict.downscale"
  
  if (plot == TRUE) {
    par.original <- par()
    par.original <- list(mfrow = par.original$mfrow, mar = par.original$mar)
    par(mfrow = c(1, 1), mar = c(5, 5, 3, 1))
    
    plot.downscale(output)
    
    par(mfrow = par.original$mfrow, mar = par.original$mar)
  }
  return(output)
}