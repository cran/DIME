\name{gng.fit}
\alias{gng.fit}
\title{
Function for Fitting GNG model parameters
}
\description{
Function to estimate parameters for GNG model, mixture of exponential and 
\emph{k}-normal. Parameters are estimated using EM algorithm.}
\usage{
gng.fit(data, avg = NULL, K = 2, weights = NULL, weights.cutoff = -1.345, 
  pi = NULL, mu = NULL, sigma = NULL, beta = NULL, tol = 1e-05, 
  max.iter = 2000, th = NULL)
}
\arguments{
  \item{data}{
an \strong{R list} of vector of normalized intensities (counts). Each element can
 correspond to particular chromosome. User can construct their own list 
 containing only the chromosome(s) they want to analyze.
}
\item{avg}{
optional vector of mean data (or log intensities). Only required when any one
of huber weight (lower, upper or full) is selected.
}
  \item{K}{
optional number of normal component that will be fitted in GNG model.
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
optional vector containing initial estimates for proportion of the GNG mixture 
components. The first and last entries are for the 
estimates of negative and positive exponentials, respectively. The middle k entries
are for normal components.
}
  \item{mu}{
optional vector containing initial estimates of the Gaussian means in GNG model.
}
  \item{sigma}{
optional vector containing initial estimates of the Gaussian standard deviation
 in GNG model. Must have K entries. 
}
  \item{beta}{
optional vector containing initial estimates for the shape parameter in 
exponential components in GNG model. Must have 2 entries, one for negative 
exponential the other for positive exponential components.
}
  \item{tol}{
optional threshold for convergence for EM algorithm to estimate GNG parameters.
}
  \item{max.iter}{
optional maximum number of iterations for EM algorithm to estimate GNG parameters.
}
  \item{th}{
optional location parameter used to fit the negative and positive exponential model.
}
}
\value{
A list of object:
\item{name}{the name of the model "GNG"}
\item{pi}{a vector of estimated proportion of each components in the model}  
\item{mu}{a vector of estimated Gaussian means for k-normal components.}
\item{sigma}{a vector of estimated Gaussian standard deviation for k-normal 
components.}
\item{beta}{a vector of estimated exponential shape values.}
\item{th1}{negative location parameter used to fit the negative exponential model.}
\item{th2}{positive location parameter used to fit the positive exponential model.}
\item{K}{the number of normal components in the corresponding mixture model.}
\item{loglike}{the log likelihood for the fitted mixture model.}
\item{iter}{the actual number of iterations run by the EM algorithm.}
\item{fdr}{the local false discover rate estimated based on GNG model.}
\item{phi}{a matrix of estimated GNG mixture component function.}
\item{AIC}{Akaike Information Criteria.}
\item{BIC}{Bayesian Information Criteria.}     
}
\author{
Cenny Taslim \email{taslim.2@osu.edu}, with contributions from Abbas Khalili
\email{khalili@stat.ubc.ca}, Dustin Potter \email{potterdp@gmail.com}, 
and Shili Lin \email{shili@stat.osu.edu}
}

\seealso{
\code{\link{DIME}}, \code{\link{inudge.fit}}, \code{\link{nudge.fit}}
}
\examples{
library(DIME)
# generate simulated datasets with underlying exponential-normal components
N1 <- 1500; N2 <- 500; K <- 4; rmu <- c(-2.25,1.50); rsigma <- c(1,1); 
rpi <- c(.05,.45,.45,.05); rbeta <- c(12,10);
set.seed(1234)
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

# Getting the best fitted GNG model (parameters)
test$pi # estimated proportion of each component in GNG
test$mu # estimated mean of the normal component(s) GNG
# estimated standard deviation of the normal component(s) in GNG
test$sigma 
# estimated shape parameter of the exponential components in best model  
test$beta
}
\keyword{ models }

