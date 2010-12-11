\name{inudge.plot.comp}
\alias{inudge.plot.comp}
\title{
Plot iNUDGE Individual Components
}
\description{
Plot each estimated individual components of iNUDGE model 
(mixture of uniform and \emph{k}-normal) fitted using \code{\link{inudge.fit}}.
}
\usage{
inudge.plot.comp(data, obj, new.plot = TRUE, legpos = NULL, xlim = NULL,
  ylim = NULL, xlab = NULL, ylab = NULL, main = NULL, lwd = NULL,...)
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
  \item{new.plot}{
optional logical variable on whether to create new plot.
}
  \item{legpos}{
optional vector of (x,y) location for the legend position
}
  \item{xlim}{
optional x-axis limit (see \code{\link{par}}).
}
  \item{ylim}{
optional y-axis limit (see \code{\link{par}}).
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
  \item{lwd}{
optional line width for lines in the plot (see \code{\link{par}}).
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
\seealso{
 \code{\link{inudge.plot.mix}}, \code{\link{inudge.plot.comp}},
 \code{\link{inudge.plot.fit}}, \code{\link{inudge.plot.qq}}, 
 \code{\link{DIME.plot.fit}}, \code{\link{gng.plot.fit}}.
}
\examples{
library(DIME);
# generate simulated datasets with underlying uniform and 2-normal distributions
set.seed(12);
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
set.seed(12);
bestInudge <- inudge.fit(data, K = 2, max.iter=20);

# plot individual components of iNUDGE
inudge.plot.comp(data,bestInudge);
# plot individual components of iNUDGE an it's mixture component on the same plot
inudge.plot.mix(bestInudge,resolution=1000,new.plot=FALSE);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dplot }
\keyword{ aplot }