\name{nudge.classify}
\alias{nudge.classify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classification Based on NUDGE Model
}
\description{
Classifies observed data into differential and non-differential groups based on
NUDGE model.
}
\usage{
nudge.classify(data, obj, obj.cutoff = 0.1, obj.sigma.diff.cutoff = NULL,
  obj.mu.diff.cutoff = NULL)
}
\arguments{
  \item{data}{
an \strong{R list} of vector of normalized intensities (counts). Each element can
 correspond to particular chromosome. User can construct their own list
 containing only the chromosome(s) they want to analyze.
}
  \item{obj}{
a list object returned by \code{\link{nudge.fit}} function.
}
  \item{obj.cutoff}{
optional local \emph{fdr} cutoff for classifying data into differential
and non-differential groups based on NUDGE model.
}
  \item{obj.sigma.diff.cutoff}{
optional cut-off for standard deviation of the normal component in NUDGE
model to be designated as representing differential.
}
  \item{obj.mu.diff.cutoff}{
optional cut-off for standard deviation of the normal component in NUDGE
model to be designated as representing differential.
}
}
\value{
A list object passed as input with additional element $class containing vector
of classifications for all the observations in data. A classification of 1
denotes that the data is classified as differential with fdr < obj.cutoff.
\item{mu.diff.cutoff}{normal component with mean > mu.diff.cutoff
  was used to represent differential component.}
\item{sigma.diff.cutoff}{normal component with standard deviation >
  sigma.diff.cutoff was used to represent differential component.}
}
\author{
Cenny Taslim \email{taslim.2@osu.edu}, with contributions from Abbas Khalili
\email{khalili@stat.ubc.ca}, Dustin Potter \email{potterdp@gmail.com},
and Shili Lin \email{shili@stat.osu.edu}
}
\seealso{
\code{\link{nudge.fit}}
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
set.seed(1234);
test <- nudge.fit(data, max.iter=20)
# vector of classification. 1 represents differential, 0 denotes non-differential
nudgeClass <- test$class;
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
