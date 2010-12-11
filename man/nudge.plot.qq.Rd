\name{nudge.plot.qq}
\alias{nudge.plot.qq}
\title{
QQ-plot of GNG model vs. observed data
}
\description{
  Produces a QQ-plot for visual inspection of quality of fit with regards to 
  the uniform Gaussian (NUDGE) mixture model estimated using the function 
  \code{\link{nudge.fit}}
}
\usage{
nudge.plot.qq(data, obj, resolution = 10, xlab = NULL, ylab = NULL, 
  main = NULL, pch = NULL, ...)
}
\arguments{
\item{data}{
an \strong{R list} of vector of normalized intensities (counts). Each element can
 correspond to a particular chromosome. User can construct their own list 
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
optional plotting character, i.e., symbol to use (see \code{\link{par}}).
}
  \item{\dots}{
additional graphical arguments to be passed to methods (see \code{\link{par}}).
}
}
\seealso{
\code{\link{nudge.fit}}, \code{\link{qqplot}}
}
\examples{
library(DIME)

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

# QQ-plot
nudge.plot.qq(data,bestNudge);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
\keyword{ dplot }% __ONLY ONE__ keyword per line
