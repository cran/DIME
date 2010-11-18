\name{gng.classify}
\alias{gng.classify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classification Based on GNG Model
}
\description{
Classifies observed data into differential and non-differential groups based on 
GNG model.
}
\usage{
gng.classify(data, obj, obj.cutoff = 0.1, obj.sigma.diff.cutoff = NULL, 
  obj.mu.diff.cutoff = NULL)
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
  \item{obj.cutoff}{
optional local \emph{fdr} cutoff for classifying data into differential
and non-differential groups based on GNG model.
}
  \item{obj.sigma.diff.cutoff}{
optional cut-off for standard deviation of the normal component in the best 
model to be declared as representing differential.
}
  \item{obj.mu.diff.cutoff}{
optional cut-off for standard deviation of the normal component in the best 
model to be declared as representing differential.
}
}
\value{
A list object passed as input with additional element $class containing vector 
of classifications for all the observations in data. A classification of 1
denotes that the data is classified as differential with fdr < obj.cutoff.
\item{mu.diff.cutoff}{normal component with mean > mu.diff.cutoff
  will be used to represent differential component.}
\item{sigma.diff.cutoff}{normal component with standard deviation > 
  sigma.diff.cutoff will be used to represent differential component.}
}
\author{
Cenny Taslim \email{taslim.2@osu.edu}, with contributions from Abbas Khalili
\email{khalili@stat.ubc.ca}, Dustin Potter \email{potterdp@gmail.com}, 
and Shili Lin \email{shili@stat.osu.edu}
}

\seealso{
\code{\link{gng.fit}}
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

# fit GNG model with 2 normal components
test <- gng.fit(data, K = 2);
# vector of classification. 1 represents differential, 0 denotes non-differential
gngClass <- test$class;
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
