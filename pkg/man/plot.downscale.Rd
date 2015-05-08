\name{plot.downscale}
\alias{plot.downscale}
\title{Plotting of downscaled occupancy at fine grain sizes}
\usage{
\method{plot.downscale}{predict.downscale}(predict.object, xlim = NULL, ylim = NULL,
               xlab = NULL, ylab = NULL, main = NULL, lwd.obs = NULL, lwd.pred = NULL,
               col.obs = NULL, col.pred = NULL, ...)
}
\arguments{
\item{predict.object}{Output object from \code{predict.downscale} and of class 
  \code{predict.downscale}}
\item{xlim, ylim}{limits of axes. Defaults to minimum and maximum of data.}
\item{xlab, ylab, main}{axis labels and title.}
\item{lwd.obs, lwd.pred}{line width of observed and predicted occupancies (\code{default} = 2).}
\item{col.obs, col.pred}{line and point colours of observed (\code{default} = black) and predicted (\code{default} = red) occupancies.}
\item{\dots}{arguments, including graphical parameters passed to other methods.}
}
\seealso{
See \code{\link{predict.downscale}} and \code{\link{hui.downscale}} for generating \code{predict.downscale} objects.
}
\description{
A simple plotting function of \code{predict.downscale} objects. Occupancy is plotted against grain size (cell area) in log-log space. Observed occupancy at large grain sizes are plotted in black, and occupancies predicted through \code{\link{predict.downscale}} plotted in red.
}
\examples{
## example species data
data.file <- system.file("extdata", "atlas_data.txt", package = "downscale")
atlas.data <- read.table(data.file, header = TRUE)

## upgrain data (using All Presences threshold)
occupancy <- upgrain(atlas.data,
                     cell.width = 10,
                     scales = 3,
                     method = "All_Presences")

## Logistic model
logis <- downscale(occupancies = occupancy,
                   model = "Logis")

## predict occupancies at fine scales
logis.pred <- predict.downscale(mod.fit = logis,
                                new.areas = c(1, 5, 25, 100, 400, 1600, 6400))

## plot predictions
plot.downscale(logis.pred)

## change some of the plotting arguments
plot.downscale(logis.pred,
               col.obs = "blue",
               pch = 16,
               ylim = c(0.07, 0.7))
}
\author{
Charles Marsh <\email{charliem2003@gmail.com}>.
}