\name{gng.plot.mix}
\alias{gng.plot.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot GNG Mixture Component Function
}
\description{
Plot each estimated individual components of GNG mixture model fitted using 
\code{\link{gng.fit}}.
}
\usage{
gng.plot.mix(obj, amplify = 1, resolution = 100, new.plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a list object returned by \code{\link{gng.fit}} function.
}
  \item{amplify}{
optional scaling factor for visualization purposes.
}
  \item{resolution}{
optional bandwidth used to estimate the density function. Higher number makes a
 smoother curve.
}
  \item{new.plot}{
optional logical variable on whether to create new plot.
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
 \code{\link{gng.plot.mix}}, \code{\link{gng.plot.comp}}, \code{\link{gng.plot.fit}},
 \code{\link{gng.plot.qq}}, \code{\link{DIME.plot.fit}}, \code{\link{inudge.plot.fit}}. 
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

# Plot the estimated GNG model imposed on the histogram of the observed data
hist(unlist(data),freq=FALSE,breaks=100,xlim=c(-20,20))
gng.plot.mix(bestGng,resolution=1000,new.plot=FALSE,col="red");
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
\keyword{ dplot }% __ONLY ONE__ keyword per line
