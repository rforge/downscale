#' @title Ensemble modelling of multiple downscaling functions
#' @name ensemble.downscale
#' 
#' @description Predict area of occupancy at fine grain sizes for multiple 
#'   downscaling methods using \code{\link{downscale}} and 
#'   \code{\link{predict.downscale}}. The mean predicted occupancies of all 
#'   models is then calculated.
#'   
#' @param occupancy vector of observed area of occupancies in squared units 
#'   (e.g. km^{2}).
#' @param area vector of grain sizes (in same units as occupancy).
#' @param extent total area in same units as occupancy (only required for 
#'   \code{FNB} and \code{Thomas} models).
#' @param tolerance_mod only applicable for the \code{Thomas} model. The
#'   tolerance used during integration in the Thomas model during optimisation
#'   of parameters. Lower numbers allow for greater accuracy but require longer 
#'   processing times (default = \code{1e-6}).
#' @param tolerance_pred only applicable for the \code{Thomas} model. The
#'   tolerance used during the prediction stage.
#' @param newdata vector of grain sizes (in same units as occupancy) for which 
#'   area of occupancy will be predicted.
#' @param models vector of chosen downscaling models. Default \code{models = 
#'   "all"} runs all available models. See \code{\link{downscale}} for list of 
#'   available models.
#' @param plot if \code{TRUE} predictions of all models are plotted against
#'   grain size along with the mean of all models.
#' @param verbose if \code{TRUE} prints updates on modelling status.
#'   
#' @details The function is a simple ensemble technique that runs all available
#'   downscaling models for observed occupancies using \code{\link{downscale}},
#'   and uses the model outputs to predict occupancies at finer grain sizes
#'   using \code{\link{predict.downscale}}. It then calculates a simple mean of
#'   all models.
#'   
#' @return Returns a dataframe. The first column \code{cell.area} is the grain
#'   sizes used for predictions. The final column \code{Means} are the mean
#'   predictions of all models for each grain size. Intermediate columns are the
#'   predicted occupancies for the selected downscaling models.
#'   
#' @seealso See \code{\link{downscale}} and \code{\link{predict.downscale}}
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
#' @example R/Examples/Examples.R
#' 
#' @method ensemble.dowmscale
#' #@export ensemble.downscale


################################################################################
# 
# ensemble.downscale.R
# Version 1.1
# 04/02/2015
#
# Updates:
#   23/02/2015 different tolerances allowed for modelling and predicting
#   05/02/2015 plot = TRUE argument added
#   05/02/2015 improved warning control
#
# Ensemble modelling
#
# Args:
#   occupancy: vector of observed area of occupancies in squared units (e.g.
#              km^{2}).
#   area: vector of grain sizes (in same units as occupancy).
#   extent: total area in same units as occupancy tolerance_mod: The tolerance
#           used during integration in the Thomas model during optimisation of
#           parameters.
#   tolerance_pred: The tolerance used during the prediction stage. 
#   newdata: vector of grain sizes (in same units as occupancy) for which area 
#            of occupancy will be predicted.
#   models: vector of chosen downscaling models. Default models = "all" runs all
#           available models.
#   plot: if TRUE predictions of all models are plotted against grain size along
#         with the mean of all models.
#   verbose: if TRUE prints updates on modelling status.
# Returns:
#    a dataframe. The first column cell.area is the grain sizes used for 
#    predictions. The final column Means are the mean predictions of all models
#    for each grain size. Intermediate columns are the predicted occupancies for
#    the selected downscaling models.
#
################################################################################

ensemble.downscale <- function(occupancy,
                               area,
                               newdata,
                               extent,
                               tolerance_mod = 1e-6,
                               tolerance_pred = 1e-6,
                               models = "all",
                               plot = TRUE,
                               verbose = TRUE) {
  # error checking - model name is correct
  suppressWarnings(apply(as.data.frame(models), 1, function(x)
    if (x %in% c("Nachman", "PL", "Logis", "Poisson", "NB", "GNB", "INB",
                 "FNB", "Thomas","all") == FALSE) {
      stop("Model name invalid", call. = FALSE)
    }))
  
  if (length(models) == 1) {
    if (models == "all") {
      model.list <- c("Nachman","PL","Logis","Poisson","NB",
                      "GNB","INB","FNB","Thomas")
    }
    if (models != "all") {
      stop("Only one model selected: ensemble modelling not applicable", 
           call. = FALSE)
    }
  }
  if (length(models) > 1) {
    model.list <- models
  }
  
  all.predicted <- as.data.frame(matrix(NA, 
                                        ncol = (length(model.list) + 1),
                                        nrow = length(newdata)))
  colnames(all.predicted) <- c("Cell.area", model.list)
  all.predicted[, "Cell.area"] <- newdata
  
  # modelling
  for (i in 1:length(model.list)) {
    model.run<-model.list[i]
    
    if(verbose == TRUE){
      cat(paste(model.run, "model is running..."))
    }
    
    mod<-downscale(occupancy = occupancy,
                   area = area,
                   model = model.run,
                   extent = extent,
                   tolerance = tolerance_mod)
    est<-predict.downscale(object = mod,
                           newdata = newdata,
                           extent = extent,
                           tolerance = tolerance_pred,
                           plot = FALSE)
    all.predicted[, i + 1] <- est$predicted[, "Occupancy"]
    
    if(verbose == TRUE){
      cat(paste("  complete", "\n"))
    }
  }
  all.predicted$Means <- rowMeans(all.predicted[, -1], na.rm = TRUE)
  
  # plotting
  if (plot == TRUE) {
  par.original <- par()
  par.original <- list(mfrow = par.original$mfrow, mar = par.original$mar)
  par(mfrow = c(3, ceiling(length(model.list) / 3)), mar = c(5, 5, 3, 1))
  
  for (i in 1:length(model.list)) {
    plot(all.predicted[, 2] ~ all.predicted[, "Cell.area"],
         type = "n",
         log = "xy",
         xlim = c(min(c(all.predicted[, "Cell.area"], area)),
                  max(c(all.predicted[, "Cell.area"], area))),
         ylim = c(min(all.predicted[, -1]), 1),
         xlab = "Log cell area",
         ylab = "Log occupancy",
         main = paste(model.list[i], "model"))
    
    points(occupancy ~ area, type="b", lwd=2)
    points(all.predicted[, "Means"] ~ all.predicted[, "Cell.area"],
           type="b", lwd=2, col = "dark grey")
    points(all.predicted[, i + 1] ~ all.predicted[, "Cell.area"],
           type="b", lwd=2, col = "red")
  }
  par(mfrow = par.original$mfrow, mar = par.original$mar)
  }
  
  return(all.predicted)
}

