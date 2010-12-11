\name{nudge.plot.comp}
\alias{nudge.plot.comp}
\title{
Plot NUDGE Individual Components
}
\description{
Plot each estimated individual components of NUDGE model 
(mixture of uniform and 1-normal) fitted using \code{\link{nudge.fit}}.
}
\usage{
nudge.plot.comp(data, obj, new.plot = TRUE, legpos = NULL, xlim = NULL,
  ylim = NULL, xlab = NULL, ylab = NULL, main = NULL, lwd = NULL, ...)
}
\arguments{
\item{data}{
an \strong{R list} of vector of normalized intensities (counts). Each element can
 correspond to a particular chromosome. User can construct 
 their own list containing only the chromosome(s) they want to analyze.
}
\item{obj}{
a list object returned by \code{\link{nudge.fit}} function.
}
  \item{new.plot}{
an \strong{R list} of vector of normalized intensities (counts). Each object can correspond to
particular chromosome that one want to fit.
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
\seealso{
 \code{\link{nudge.plot.mix}}, \code{\link{inudge.plot.comp}}, 
 \code{\link{nudge.plot.fit}}, \code{\link{nudge.plot.qq}}, 
 \code{\link{DIME.plot.fit}}, \code{\link{gng.plot.fit}}.
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

# plot individual components of NUDGE
nudge.plot.comp(data,bestNudge);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
\keyword{ dplot }% __ONLY ONE__ keyword per line
