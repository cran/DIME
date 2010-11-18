\name{gng.plot.fit}
\alias{gng.plot.fit}

\title{
Plot GNG Goodness of Fit
}
\description{
Plot the estimated GNG mixture model fitted using \code{\link{gng.fit}} along with 
it's estimated individual components, superimposed on the histogram of the 
observation data. This plot shows how good the fit of the estimated model to the
 data.
}
\usage{
gng.plot.fit(data, obj, resolution = 100, breaks = 100, legpos = NULL,
  xlim = NULL, main=NULL, ...)
}

\arguments{
\item{data}{
an \strong{R list} of vector of normalized intensities (counts). Each element can
 correspond to a particular chromosome. User can construct  their own list 
 containing only the chromosome(s) they want to analyze.
}
  \item{obj}{
a list object returned by \code{\link{gng.fit}} function.
}
  \item{resolution}{
optional bandwidth used to estimate the density function. Higher number smoother 
curve.
}
  \item{breaks}{
optional see \code{\link{hist}}, breaks parameters for histogram plot.
}
  \item{legpos}{
optional vector of (x,y) location for the legend position
}
  \item{xlim}{
optional x-axis limit (see \code{\link{par}}).
}
  \item{main}{
optional plot title (see \code{\link{par}}).
}
  \item{\dots}{
additional graphical arguments to be passed to methods (see \code{\link{par}}).
}
}
\details{
The components representing differential data are denoted by asterisk (*)
symbol on the plot legend.
}

\seealso{
\code{\link{gng.plot.comp}}, \code{\link{gng.plot.mix}}, \code{\link{hist}}
}
\examples{
library(DIME)
# generate simulated datasets with underlying exponential-normal components
N1 <- 1500; N2 <- 500; K <- 4; rmu <- c(-2.25,1.50); rsigma <- c(1,1); 
rpi <- c(.05,.45,.45,.05); rbeta <- c(12,10);
set.seed(1234);
chr1 <- c(-rgamma(ceiling(rpi[1]*N1),shape = 1,scale = rbeta[1]), 
  rnorm(ceiling(rpi[2]*N1),rmu[1],rsigma[1]), 
  rnorm(ceiling(rpi[3]*N1),rmu[2],rsigma[2]), 
  rgamma(ceiling(rpi[4]*N1),shape = 1,scale = rbeta[2]));
chr2 <- c(-rgamma(ceiling(rpi[1]*N2),shape = 1,scale = rbeta[1]), 
  rnorm(ceiling(rpi[2]*N2),rmu[1],rsigma[1]), 
  rnorm(ceiling(rpi[3]*N2),rmu[2],rsigma[2]), 
  rgamma(ceiling(rpi[4]*N2),shape = 1,scale = rbeta[2])); 
chr3 <- c(-rgamma(ceiling(rpi[1]*N2),shape = 1,scale = rbeta[1]), 
  rnorm(ceiling(rpi[2]*N2),rmu[1],rsigma[1]), 
  rnorm(ceiling(rpi[3]*N2),rmu[2],rsigma[2]), 
  rgamma(ceiling(rpi[4]*N2),shape = 1,scale = rbeta[2]));
# analyzing only chromosome 1 and chromosome 3
data <- list(chr1,chr3);

# Fitting a GNG model only with 2-normal components
bestGng <- gng.fit(data,K=2);

# Goodness of fit plot
gng.plot.fit(data,bestGng);
}
\keyword{ Graphics }
\keyword{ aplot }