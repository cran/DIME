\name{gng.plot.comp}
\alias{gng.plot.comp}
\title{
Plot GNG Individual Components
}
\description{
Plot each estimated individual components of GNG model 
(mixture of exponential and \emph{k}-normal) fitted using \code{\link{gng.fit}}.
}
\usage{
gng.plot.comp(data, obj, new.plot = TRUE, legpos = NULL, xlim=NULL,
  ylim=NULL, xlab=NULL, ylab=NULL, main=NULL,lwd=NULL,...)
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
 \code{\link{gng.plot.mix}}, \code{\link{gng.plot.comp}}, \code{\link{gng.plot.fit}},
 \code{\link{gng.plot.qq}}, \code{\link{DIME.plot.fit}}, \code{\link{inudge.plot.fit}}. 
}
\examples{
library(DIME);
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

# Fitting a GNG model with 2-normal component
bestGng <- gng.fit(data,K=2);

# plot individual components of GNG
gng.plot.comp(data,bestGng);
# plot mixture component on top of the individual components plot
gng.plot.mix(bestGng,resolution=1000,new.plot=FALSE);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dplot }
\keyword{ aplot }