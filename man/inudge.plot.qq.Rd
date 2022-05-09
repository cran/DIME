\name{inudge.plot.qq}
\alias{inudge.plot.qq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QQ-plot of GNG model vs. observed data
}
\description{
  Produces a QQ-plot for visual inspection of quality of fit with regards to 
  the uniform Gaussian (iNUDGE) mixture model estimated using the function 
  \code{\link{inudge.fit}}
}
\usage{
inudge.plot.qq(data, obj, resolution = 10, xlab = NULL, ylab = NULL,
  main = NULL, pch = NULL, ...)
}
\arguments{
  \item{data}{
an \strong{R list} of vector of normalized intensities (counts). Each element can
 correspond to particular chromosome. User can construct their own list 
 containing only the chromosome(s) they want to analyze.
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
\seealso{
\code{\link{inudge.fit}}, \code{\link{qqplot}}
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

# fit iNUDGE model with 2-normal components and maximum iteration =20 
set.seed(1234);
bestInudge <- inudge.fit(data, K=2, max.iter=20)

# QQ-plot
inudge.plot.qq(data,bestInudge);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dplot }
\keyword{ aplot }% __ONLY ONE__ keyword per line
