\name{inudge.plot.fit}
\alias{inudge.plot.fit}
\title{
Plot iNUDGE Goodness of Fit
}
\description{
Plot the estimated iNUDGE mixture model fitted using \code{\link{inudge.fit}} 
along with it's estimated individual components, superimposed on the histogram 
of the observation data. This plot shows how good the fit of the estimated model
 to the data.
}
\usage{
inudge.plot.fit(data, obj, resolution = 100, breaks = 100, 
legpos = NULL, xlim = NULL, main = NULL,...)
}
\arguments{
  \item{data}{
an \strong{R list} of vector of normalized intensities (counts). Each element can
 correspond to particular chromosome. User can construct their own list 
 containing only the chromosome(s) they want to analyze.
}
  \item{obj}{
a list object returned by \code{\link{inudge.fit}} function.
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
\code{\link{gng.plot.comp}}, \code{\link{gng.plot.mix}}, 
\code{\link{hist}}
}
\examples{
library(DIME);
# generate simulated datasets with underlying uniform and 2-normal distributions
set.seed(1234);
N1 <- 1500; N2 <- 500; rmu <- c(-2.25,1.5); rsigma <- c(1,1); 
rpi <- c(.10,.45,.45); a <- (-6); b <- 6; 
chr4 <- list(c(-runif(ceiling(rpi[1]*N1),min = a,max =b),
  rnorm(ceiling(rpi[2]*N1),rmu[1],rsigma[1]), 
  rnorm(ceiling(rpi[3]*N1),rmu[2],rsigma[2])));
chr9 <- list(c(-runif(ceiling(rpi[1]*N2),min = a,max =b),
  rnorm(ceiling(rpi[2]*N2),rmu[1],rsigma[1]), 
  rnorm(ceiling(rpi[3]*N2),rmu[2],rsigma[2])));
# analyzing chromosome 4 and 9
data <- list(chr4,chr9);

# fit iNUDGE model with 2-normal components and maximum iterations = 20
set.seed(1234);
bestInudge <- inudge.fit(data, K = 2, max.iter=20);

# Goodness of fit plot
inudge.plot.fit(data,bestInudge,legpos=c(-6,0.3),ylim=c(0,0.3),breaks=40);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
\keyword{ dplot }% __ONLY ONE__ keyword per line
