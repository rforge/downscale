#' @name predict.downscale
#' @method plot predict.downscale
#' @aliases plot
#' @param predict.object an object from \code{\link{predict.downscale}}.
#' @param \dots arguments, including graphical parameters for plot.downscale,
#'   passed to other methods.
#' @rdname predict.downscale



################################################################################
# 
# plot.downscale.R
# Version 1.0
# 03/02/2015
#
# Plot the observed and predicted area of occupancy against grain size on
# log-log axes.
#
# Args:
#   predict.object: an object of class 'predict.downscale' containing observed
#                   and predicted data 
#   ...: arguments, including graphical parameters, passed to other methods.
#
# Returns:
#   no object returned.
#
################################################################################

plot.downscale <- function(predict.object, ...) {
  # error checking
  if (class(predict.object) != "predict.downscale"){
    stop("Input data not of class 'predict.downscale'")
  }
  
  observed <- predict.object$observed
  predicted <- predict.object$predicted
  
  plot(observed[, "Occupancy"] ~ observed[, "Cell.area"],
       type = "n",
       log = "xy",
       xlim = c(min(c(observed[, "Cell.area"], 
                      predicted[, "Cell.area"])),
                max(c(observed[, "Cell.area"], 
                      predicted[, "Cell.area"]))),
       ylim = c(min(c(observed[, "Occupancy"], 
                      predicted[, "Occupancy"])), 1),
       xlab = "Log cell area",
       ylab = "Log occupancy",
       main = paste(model.run, "model"))
  
  points(observed[, "Occupancy"] ~ observed[, "Cell.area"],
         type="b",
         lwd=2)
  
  points(predicted[, "Occupancy"] ~ predicted[, "Cell.area"],
         type = "b",
         col = "red",
         lwd = 2)
}












