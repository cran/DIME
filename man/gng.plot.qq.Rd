\name{gng.plot.qq}
\alias{gng.plot.qq}
\title{
QQ-plot of GNG model vs. observed data
}
\description{
  Produces a QQ-plot for visual inspection of quality of fit with regards to 
  the exponential Gaussian (GNG) mixture model estimated using the function 
  \code{\link{gng.fit}}
}
\usage{
gng.plot.qq(data, obj, resolution=10, xlab=NULL, ylab=NULL, 
  main=NULL, pch=NULL,...)
}
\arguments{
\item{data}{
an \strong{R list} of vector of normalized intensities (counts). Each element can
 correspond to a particular chromosome. User can construct 
 their own list containing only the chromosome(s) they want to analyze.
}
  \item{obj}{
a list object returned by \code{\link{gng.fit}} function.
}
  \item{resolution}{
optional number of points used to sample the estimated density function.
}
\item{xlab}{
optional x-axis label (see \code{\link{par}}).
}
\item{ylab}{
optional y-axis label (see \code{\link{par}}).
}
\item{main}{
optional plot title (see \code{\link{par}}).
}
\item{pch}{
optional plotting symbol to use (see \code{\link{par}}).
}
\item{\dots}{
additional graphical arguments to be passed to methods (see \code{\link{par}}).
}
}
\author{
Cenny Taslim \email{taslim.2@osu.edu}, with contributions from Abbas Khalili
\email{khalili@stat.ubc.ca}, Dustin Potter \email{potterdp@gmail.com}, 
and Shili Lin \email{shili@stat.osu.edu}
}
\seealso{
\code{\link{gng.fit}}, \code{\link{gng.plot.fit}}
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

# Fitting a GNG model only
bestGng <- gng.fit(data,K=2);

# QQ-plot
gng.plot.qq(data,bestGng)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dplot }
\keyword{ aplot}