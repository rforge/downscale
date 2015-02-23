\name{Downscale-package}
\alias{downscale-package}
\docType{package}
\title{
Downscaling species occupancy
}
\description{
A set of functions that downscale species occupancy from coarse grain sizes to 
predict species occupancy at finer grain sizes.
}
\details{
\tabular{ll}{
Package: \tab downscale\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-23\cr
License: \tab GPL-2\cr
}
The package provides two main functions. \code{\link{downscale}} estimates the parameters for each downscaling function for log occupancy against grain size. \code{\link{predict.downscale}} then takes these parameters and estimates occupancy at finer grain sizes. A third function, \code{\link{ensemble.downscale}} will run \code{\link{downscale}}and \code{\link{predict.downscale}} for a number of selected downscaling functions and calculate the mean predicted occupancies across all models.

Occupancy should be measured as the proportion of occupied cells. Grain size should be measured as cell area and total extent should be in the same units.
}
\author{
Charles Marsh

Maintainer: Charles Marsh <charliem2003@gmail.com>
}
\references{
Azaele, S., Cornell, S.J., & Kunin, W.E. (2012). Downscaling species occupancy from coarse spatial scales. \emph{Ecological Applications} 22, 1004-1014.

Barwell, L.J., Azaele, S., Kunin, W.E., & Isaac, N.J.B. (2014). Can coarse-grain patterns in insect atlas data predict local occupancy? \emph{Diversity and Distributions} 20, 895-907.}

