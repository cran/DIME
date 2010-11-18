\name{nudge.plot.mix}
\alias{nudge.plot.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot NUDGE Mixture Component Function
}
\description{
Plot each estimated individual components of NUDGE mixture model fitted using 
\code{\link{nudge.fit}}.
}
\usage{
nudge.plot.mix(obj, amplify = 1, resolution = 100, new.plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a list object returned by \code{\link{nudge.fit}} function.
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
 \code{\link{nudge.plot.comp}}, \code{\link{nudge.plot.fit}},
 \code{\link{nudge.plot.qq}}, \code{\link{DIME.plot.fit}}, \code{\link{gng.plot.fit}}. 
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

# fit NUDGE model with maximum iterations = 20 only
bestNudge <- nudge.fit(data, max.iter=20);

# plot estimated iNUDGE model imposed on the histogram of observed data
hist(unlist(data),freq=FALSE,breaks=40);
nudge.plot.mix(bestNudge,resolution=1000,new.plot=FALSE,col="red");

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
\keyword{ dplot }% __ONLY ONE__ keyword per line
