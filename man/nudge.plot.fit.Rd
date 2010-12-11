\name{nudge.plot.fit}
\alias{nudge.plot.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot NUDGE Goodness of Fit
}
\description{
Plot the estimated NUDGE mixture model fitted using \code{\link{nudge.fit}} 
along with it's estimated individual components, superimposed on the histogram 
of the observation data. This plot shows how good the fit of the estimated model
 to the data.
}
\usage{
nudge.plot.fit(data, obj, resolution = 100, breaks = 100, 
xlim = NULL, legpos = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
an \strong{R list} of vector of normalized intensities (counts). Each element can
 correspond to a particular chromosome. User can construct their own list 
 containing only the chromosome(s) they want to analyze.
}
  \item{obj}{
a list object returned by \code{\link{nudge.fit}} function.
}
  \item{resolution}{
optional bandwidth used to estimate the density function. Higher number smoother 
curve.
}
  \item{breaks}{
optional see \code{\link{hist}}, breaks parameters for histogram plot.
}
  \item{xlim}{
optional limit for the x-axis.
}
  \item{legpos}{
optional vector of (x,y) location for the legend position
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
\code{\link{nudge.plot.comp}}, \code{\link{nudge.plot.mix}}, \code{\link{hist}}
}
\examples{
library(DIME);
# generate simulated datasets with underlying uniform and 1-normal components
set.seed(1234);
N1 <- 1500; N2 <- 500; rmu <- c(1.5); rsigma <- c(1); 
rpi <- c(.10,.90); a <- (-6); b <- 6; 
chr1 <- c(-runif(ceiling(rpi[1]*N1),min = a,max =b),
  rnorm(ceiling(rpi[2]*N1),rmu[1],rsigma[1]));
chr4 <- c(-runif(ceiling(rpi[1]*N2),min = a,max =b),
  rnorm(ceiling(rpi[2]*N2),rmu[1],rsigma[1]));  
# analyzing chromosome 1 and 4
data <- list(chr1,chr4);

# fit NUDGE model with maximum iterations = 20
set.seed(1234);
bestNudge <- nudge.fit(data, max.iter=20);

# goodness of fit plot
nudge.plot.fit(data,bestNudge,breaks=40);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
\keyword{ dplot }% __ONLY ONE__ keyword per line
