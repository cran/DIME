\name{inudge.plot.mix}
\alias{inudge.plot.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot iNUDGE Mixture Component Function
}
\description{
Plot each estimated individual components of iNUDGE mixture model fitted using 
\code{\link{inudge.fit}}.
}
\usage{
inudge.plot.mix(obj, amplify = 1, resolution = 100, new.plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a list object returned by \code{\link{inudge.fit}} function.
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
\seealso{
 \code{\link{inudge.plot.comp}}, \code{\link{inudge.plot.fit}},
 \code{\link{inudge.plot.qq}}, \code{\link{DIME.plot.fit}}, \code{\link{gng.plot.fit}}. 
}
\examples{
library(DIME)

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

# fit iNUDGE model
set.seed(1234);
bestInudge <- inudge.fit(data, K = 2, max.iter=20);

# plot estimated iNUDGE model imposed on the histogram of observed data
hist(unlist(data),freq=FALSE,breaks=40);
inudge.plot.mix(bestInudge,resolution=1000,new.plot=FALSE,col="red");
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dplot }
\keyword{ aplot }% __ONLY ONE__ keyword per line
