\name{nudge.fit}
\alias{nudge.fit}
\title{
Function for Fitting NUDGE model parameters
}
\description{
Function to estimate parameters for both NUDGE model, mixture of 
 uniform and 1-normal. Parameters are estimated using EM algorithm.
}
\usage{
nudge.fit(data, avg = NULL, weights = NULL, weights.cutoff = -1.345,
  pi = NULL, mu = NULL, sigma = NULL, tol = 1e-5, max.iter = 2000, z = NULL)
}
\arguments{
  \item{data}{
an \strong{R list} of vector of normalized intensities (counts). Each element can
 correspond to a particular chromosome. User can construct their own list 
 containing only the chromosome(s) they want to analyze.
}
  \item{avg}{
optional vector of mean data (or log intensities). Only required when any one
of huber weight (lower, upper or full) is selected.
}
  \item{weights}{
optional weights to be used for robust fitting. Can be a matrix the same 
  length as data, or a character description of the huber weight method to be 
  employed:
  "lower" - only value below weights.cutoff are weighted,\\ 
  "upper" - only value above weights.cutoff are weighted,\\
  "full"  - both values above and below weights.cutoff are weighted,\\ 
  If selected, mean of data (avg) is required. 
}
  \item{weights.cutoff}{
optional cutoff to be used with the Huber weighting scheme.
}
  \item{pi}{
optional vector containing initial estimates for proportion of the NUDGE mixture 
components. The first entry is for the uniform component, the middle \emph{k} entries
are for normal components.
}
  \item{mu}{
optional vector containing initial estimates of the Gaussian means in NUDGE model.
}
  \item{sigma}{
optional vector containing initial estimates of the Gaussian standard deviation
 in (i)NUDGE model. Must have K entries. 
}
  \item{tol}{
optional threshold for convergence for EM algorithm to estimate NUDGE parameters.
}
  \item{max.iter}{
optional maximum number of iterations for EM algorithm to estimate NUDGE parameters.
}
  \item{z}{
optional 2-column matrix with each row giving initial estimate of probability of
the region being non-differential and a starting estimate for the probability of
the region being differential. Each row must sum to 1. Number of row must be 
equal to data length.
}
}
\value{
A list of object:
\item{name}{the name of the model "NUDGE"}
\item{pi}{a vector of estimated proportion of each components in the model}  
\item{mu}{a vector of estimated Gaussian means for k-normal components.}
\item{sigma}{a vector of estimated Gaussian standard deviation for k-normal 
components.}
\item{loglike}{the log likelihood for the fitted mixture model.}
\item{iter}{the actual number of iterations run by the EM algorithm.}
\item{fdr}{the local false discover rate estimated based on NUDGE model.}
\item{phi}{a matrix of estimated NUDGE mixture component function.}
\item{AIC}{Akaike Information Criteria.}
\item{BIC}{Bayesian Information Criteria.}
}
\author{
Cenny Taslim \email{taslim.2@osu.edu}, with contributions from Abbas Khalili
\email{khalili@stat.ubc.ca}, Dustin Potter \email{potterdp@gmail.com}, 
and Shili Lin \email{shili@stat.osu.edu}
}
\seealso{
\code{\link{DIME}}, \code{\link{gng.fit}}, \code{\link{inudge.fit}}
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
bestNudge <- nudge.fit(data, max.iter=20);

# Getting the best fitted NUDGE model (parameters)
bestNudge$pi # estimated proportion of each component in NUDGE
bestNudge$mu # estimated mean of the normal component(s) in NUDGE
# estimated standard deviation of the normal component(s) in NUDGE
bestNudge$sigma 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
