\name{DIME.plot.fit}
\alias{DIME.plot.fit}
\alias{DIME.plot}
\title{
Plot Best Model Goodness of Fit
}
\description{
Plot the best mixture model fitted using \code{\link{DIME}} along with their 
estimated individual components, superimposed on the histogram of the 
observation data. This plot shows how good the fit of the estimated model to the
 data.
}
\usage{
DIME.plot.fit(data, obj, ...)
}
\arguments{
  \item{data}{
an \strong{R list} of vector of normalized intensities (counts). Each element can
 correspond to particular chromosome. User can construct their own list 
 containing only the chromosome(s) they want to analyze.
}
  \item{obj}{
a list object returned by \code{\link{DIME}} function.
}
  \item{\dots}{
additional graphical arguments to be passed to methods (see \code{\link{par}}).
}
}
\details{
The components representing differential data are denoted by asterisk (*)
symbol on the plot legend.
}
\author{
Cenny Taslim \email{taslim.2@osu.edu}, with contributions from Abbas Khalili
\email{khalili@stat.ubc.ca}, Dustin Potter \email{potterdp@gmail.com}, 
and Shili Lin \email{shili@stat.osu.edu}
}
\seealso{\code{\link{DIME}}, \code{\link{gng.plot.fit}},\code{\link{inudge.plot.fit}} }
\examples{
library(DIME);
# generate simulated datasets with underlying exponential-normal components
N1 <- 1500; N2 <- 500; K <- 4; rmu <- c(-2.25,1.50); rsigma <- c(1,1) 
rpi <- c(.05,.45,.45,.05); rbeta <- c(12,10)
set.seed(1234)
chr1 <- c(-rgamma(ceiling(rpi[1]*N1),shape = 1,scale = rbeta[1]), 
  rnorm(ceiling(rpi[2]*N1),rmu[1],rsigma[1]), 
  rnorm(ceiling(rpi[3]*N1),rmu[2],rsigma[2]), 
  rgamma(ceiling(rpi[4]*N1),shape = 1,scale = rbeta[2]));
chr3 <- c(-rgamma(ceiling(rpi[1]*N2),shape = 1,scale = rbeta[1]), 
  rnorm(ceiling(rpi[2]*N2),rmu[1],rsigma[1]), 
  rnorm(ceiling(rpi[3]*N2),rmu[2],rsigma[2]), 
  rgamma(ceiling(rpi[4]*N2),shape = 1,scale = rbeta[2]));
data <- list(chr1,chr3);

# run DIME with small maximum iterations and repetitions
set.seed(1234);
test <- DIME(data,gng.max.iter=10,gng.rep=1,inudge.max.iter=10,inudge.rep=1,
 nudge.max.iter=10,nudge.rep=1);
 
# plot best model
DIME.plot.fit(data,test); 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Graphics }
\keyword{ aplot }
