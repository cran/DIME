\name{DIME}
\alias{DIME}
\alias{DIME}
\title{
DIME (Differential Identification using Mixtures Ensemble) 
}
\description{
A robust differential identification method that considers ensemble of finite 
mixture models combined with a local false discovery rate \emph{(fdr)} for 
analyzing ChIP-seq data comparing two samples. \cr
This package can also be used to identify differential
in other high throughput data such as microarray, methylation etc.
}
\usage{
DIME(data, avg = NULL, gng.K = 2, gng.weights = NULL, gng.weights.cutoff= -1.345,
  gng.pi = NULL, gng.mu = NULL, 
  gng.sigma = NULL, gng.beta = NULL, gng.tol = 1e-05, gng.max.iter = 2000, 
  gng.th = NULL, gng.rep = 15, gng.fdr.cutoff = 0.1, 
  gng.sigma.diff.cutoff = NULL, gng.mu.diff.cutoff = NULL, 
  gng.var.thres = 1e2, gng.min.sd = NULL,
  inudge.K = 2, inudge.weights = NULL, inudge.weights.cutoff = -1.345, 
  inudge.pi = NULL, inudge.mu = NULL,
  inudge.sigma = NULL, inudge.tol = 1e-05, inudge.max.iter = 2000,
  inudge.z = NULL, inudge.rep = 15, inudge.fdr.cutoff = 0.1, 
  inudge.sigma.diff.cutoff = NULL, inudge.mu.diff.cutoff = NULL, 
  inudge.var.thres = 1e2, inudge.min.sd = NULL,
  nudge.z = NULL, nudge.tol = 1e-05, nudge.max.iter = 2000, 
  nudge.mu = NULL, nudge.sigma = NULL, nudge.rep = 15, 
  nudge.fdr.cutoff = 0.1, nudge.weights = NULL, nudge.weights.cutoff = -1.345, 
  nudge.pi = NULL)
}
\arguments{
  \item{data}{
an \strong{R list} of vector of normalized difference (log ratios). Each element 
can correspond to a particular chromosome. User can construct their own list 
containing only the chromosome(s) they want to analyze.
}
\item{avg}{
optional \strong{R list} of vector of mean data (or log intensities). 
Each element can correspond to a particular chromosome in data. 
Only required when any one of huber weight (lower, upper or full) is selected.
}
  \item{gng.K}{
optional maximum number of normal component that will be fitted in GNG model. 
For example: gng.K=2 will fit a model with 1 and 2 normal components and select 
the best k.
}
  \item{gng.weights}{
optional weights to be used for robust fitting. Can be a matrix the same 
  length as data with each row correspond to weights to be used in each 
  repetition or a character description of the huber-type method to be 
  employed:
  "lower" - only value below cutoff are weighted,\\ 
  "upper" - only value above cutoff are weighted,\\
  "full"  - both values above and below the cutoff are weighted,\\ 
  If selected, mean of data (avg) is required. 
}
  \item{gng.weights.cutoff}{
optional cutoff to be used with the Huber weighting scheme.
}
  \item{gng.pi}{
optional matrix containing initial estimates for proportion of the GNG mixture 
components. Each row is the initial pi to be used in each repetition. 
Each row must have gng.K+2 entries. The first and last entries are for the 
estimates of negative and positive exponentials, respectively. The middle k 
entries are for normal components. 
}
  \item{gng.mu}{
optional matrix containing initial estimates of the Gaussian means in GNG model. 
Each row is the initial means to be used in each repetition. Each row must 
have gng.K entries.
}
  \item{gng.sigma}{
optional maxtrix containing initial estimates of the Gaussian standard deviation
 in GNG model. Each row is the initial means to be used in each repetition. 
 Each row must have gng.K entries.
}
  \item{gng.beta}{
optional maxtrix containing initial estimates for the shape parameter in 
exponential components in GNG model. Each row is the initial beta's to be used 
in each repetition. Each row must have 2 entries, one for negative exponential 
followed by another for positive exponential.
}
  \item{gng.tol}{
optional threshold for convergence for EM algorithm to estimate GNG parameters.
}
  \item{gng.max.iter}{
optional maximum number of iterations for EM algorithm to estimate GNG parameters.
}
  \item{gng.th}{
optional 2-column matrix of threshold for the two location exponential components.
 First column is the initial estimates for negative exponential and the second 
 column is the initial estimates for positive exponential.
}
  \item{gng.rep}{
optional number of times to repeat the GNG parameter estimation using different 
starting estimates.
}
  \item{gng.fdr.cutoff}{
optional cut-off for local \emph{fdr} for classifying regions into differential
and non-differential using GNG mixture.
}
  \item{gng.sigma.diff.cutoff}{
optional cut-off for sigma of the normal component in GNG to be declared as 
representing differential. For example: gng.sigma.diff.cutoff = 2 then if a 
normal component has sigma > 2 then this component is considered as differential
 component. Default = (1.5*iqr(data)-gng$mu)/2.
Where gng$mu is mean of non-differential normal components in iNUDGE.
}
  \item{gng.mu.diff.cutoff}{
optional cut-off for mu of the normal component in GNG to be declared as 
representing differential. For example: gng.mu.diff.cutoff = 2 then if a normal 
component has mean > 2 then this component is considered as differential 
component.
}
\item{gng.var.thres}{
optional threshold to detect huge imbalance in variance. 
max(gng.variance)/min(gng.variance) <= gng.var.thres.
}
\item{gng.min.sd}{
optional threshold to detect very small sigma. 
all normal components in GNG model has to have sigma > gng.min.sd.
Default = 0.1 * sd(data)
}
  \item{inudge.K}{
optional maximum number of normal component that will be fitted in iNUDGE model. 
For example: inudge.K=2 will fit a model with 1 and 2 normal components and select 
the best k.
}
  \item{inudge.weights}{
optional weights to be used for robust fitting. Can be a matrix the same 
  length as data with each row correspond to weights to be used in each 
  repetition or a character description of the huber-type method to be 
  employed:
  "lower" - only value below cutoff are weighted,\\ 
  "upper" - only value above cutoff are weighted,\\
  "full"  - both values above and below the cutoff\\ 
            are weighted, 
  any other character - "lower" are used (default). \\
  If selected, mean of data (avg) is required.
}
  \item{inudge.weights.cutoff}{
optional cutoff to be used with the Huber weighting scheme.
}
  \item{inudge.pi}{
optional matrix of initial estimates for proportion of the iNUDGE mixture components.
Each row correspond to the intial proportion to be used in each repetition. Each 
row must have inudge.K+1 entries corresponding to proportion of negative 
exponential, proportion of k-normal and proportion of exponential, respectively.
}
  \item{inudge.mu}{
optional maxtrix of initial estimates of the Gaussian means in iNUDGE model. 
Each row correspond to the intial means to be used in each repetition. Each row 
must have inudge.K entries.
}
  \item{inudge.sigma}{
optional matrix of initial estimates for Gaussian standard deviation in iNUDGE model. 
Each row correspond to the intial means to be used in each repetition. Each row 
must have inudge.K entries.
}
  \item{inudge.tol}{
optional threshold for convergence for EM algorithm to estimate iNUDGE parameters.
}
  \item{inudge.max.iter}{
optional maximum number of iterations for EM algorithm to estimate iNUDGE parameters.
}
  \item{inudge.z}{
optional 2-column matrix with each row giving initial estimate of probability of
the region being non-differential and a starting estimate for the probability of
the region being differential. Each row must sum to 1. Number of row must be 
equal to data length.
}
  \item{inudge.rep}{
optional number of times to repeat the iNUDGE parameter estimation using different 
starting estimates.
}
  \item{inudge.fdr.cutoff}{
optional cut-off for local \emph{fdr} for classifying regions into differential
and non-differential based on iNUDGE mixture.
}
  \item{inudge.sigma.diff.cutoff}{
optional cut-off for sigma of the normal component in GNG to be declared as 
representing differential. For example: gng.sigma.diff.cutoff = 2 then if a 
normal component has sigma > 2 then this component is considered as 
differential component. Default = (1.5*iqr(data)-inudge$mu^)/2.
Where inudge$mu^ is mean of non-differential normal components in iNUDGE.
}
  \item{inudge.mu.diff.cutoff}{
optional cut-off for mu of the normal component in GNG to be declared as
representing differential. For example: gng.mu.diff.cutoff = 2 then if a 
normal component has mean > 2 then this component is considered as 
differential component.
}
\item{inudge.var.thres}{
optional threshold to detect huge imbalance in variance. 
max(inudge.variance)/min(inudge.variance) <= inudge.var.thres.
}
\item{inudge.min.sd}{
optional threshold to detect very small sigma. 
all normal components in iNUDGE model has to have sigma > inudge.min.sd.
Default = 0.1 * sd(data)
}
  \item{nudge.z}{
optional 2-column matrix with each row giving initial estimate of probability of
the region being non-differential and a starting estimate for the probability of
the region being differential. Each row must sum to 1. Number of row must be 
equal to data length.
}
  \item{nudge.tol}{
optional threshold for convergence for EM algorithm to estimate NUDGE parameters.
}
  \item{nudge.max.iter}{
optional maximum number of iterations for EM algorithm to estimate iNUDGE 
parameters.
}
  \item{nudge.mu}{
optional maxtrix of initial estimates of the Gaussian means in NUDGE model. 
Each row correspond to the intial means to be used in each repetition. Each row 
must have 1 entry.
}
  \item{nudge.sigma}{
optional initial estimates of the Gaussian standard deviation in NUDGE model. 
Each row correspond to the intial standard deviation to be used in each 
repetition. Each row must have 1 entry.
}
  \item{nudge.rep}{
optional number of times to repeat the NUDGE parameter estimation using 
different starting estimates.
}
  \item{nudge.fdr.cutoff}{
optional cut-off for local \emph{fdr} for classifying regions into differential
and non-differential based on NUDGE mixture.
}
  \item{nudge.weights}{
optional weights to be used for robust fitting. Can be a matrix the same 
  length as data with each row correspond to weights to be used in each 
  repetition or a character description of the huber-type method to be 
  employed:
  "lower" - only value below cutoff are weighted,\\ 
  "upper" - only value above cutoff are weighted,\\
  "full"  - both values above and below the cutoff\\ 
            are weighted, 
  any other character - "lower" are used (default). \\
  If selected, mean of data (avg) is required. 
}
  \item{nudge.weights.cutoff}{
optional cutoff to be used with the Huber weighting scheme.
}
  \item{nudge.pi}{
optional initial estimates for proportion of the NUDGE mixture components.
Each row is the initial pi to be used in each repetition. Each row must have 2 
entries: proportion of uniform and proportion of normal components, respectively . 
}
}
\details{
After normalization, a Gamma-Normal(k)-Gamma (GNG), a Uniform-Normal(k) (iNUDGE)
and a Uniform-Normal (NUDGE) mixture are fitted to the data. Two-phase selection
method is used to choose the best model. 
The (k)-normal component can represent either differential
regions or non-differential regions depending on their location and shape, 
making the model more robust to different underlying distributions. The
exponential or uniform represents differential sites. Local \emph{(fdr)} is
computed from the best fitted model. Parameters estimation is performed 
using EM algorithm.
}
\value{
A list with 4 components (i.e. best, gng, inudge and nudge) which in itself is 
another list containing the estimated parameters of each model fitted 
correspondingly. "best" lists the model chosen as the best overall model, 
i.e. if the best model is inudge then best$name = "iNUDGE" and its content 
is the same as inudge. Thus, depending on the model, the components are:
\item{name}{the name of the model {"GNG", "iNUDGE","NUDGE"} where
  GNG: normal(k)-exponential (a special case of gamma), iNUDGE: normal(k)-uniform, 
  or NUDGE: normal-uniform models}
\item{pi}{a vector of estimated proportion of each components in the model}  
\item{mu}{a vector of estimated Gaussian means for k-normal components.}
\item{sigma}{a vector of estimated Gaussian standard deviation for k-normal 
components.}
\item{beta}{a vector of estimated exponential shape values. Only available in gng. }
\item{th1}{negative location parameter used to fit the negative exponential model.
  Only available in gng.}
\item{th2}{positive location parameter used to fit the positive exponential model.
  Only available in gng.}
\item{a}{the minimum value of the normalized data. Only available in (i)nudge.}
\item{b}{the maximum value of the normalized data. Only available in (i)nudge.}
\item{K}{the number of normal components in the corresponding mixture model. For
 inudge, K=1.}
\item{loglike}{the log likelihood for the fitted mixture model.}
\item{iter}{the number of iterations run by the EM algorithm until either 
  convergence or iteration limit was reached.}
\item{fdr}{the local false discover rate estimated based on the corresponding 
  model.}
\item{class}{a vector of classifications for the observations in data. 
  A classification of 0 denotes that the regions could not be classified as 
  differential with fdr < <model>.fdr.cutoff, 1 denotes differential.}
\item{diffPiIdx}{a vector of index of the normal components that are defined as
  capturing differential regions based on their shape and locations.}
\item{phi}{a vector of estimated mixture function}
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
\code{\link{gng.fit}}, \code{\link{inudge.fit}}, \code{\link{nudge.fit}}
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

# run DIME
set.seed(1234)
test <- DIME(data,gng.max.iter=10,gng.rep=1,inudge.max.iter=10,inudge.rep=1,
 nudge.max.iter=10,nudge.rep=1)

# Getting the best fitted model (parameters)
test$best$name # name of the best fitted model
test$best$pi # estimated proportion of each component in the best model
test$best$mu # estimated mean of the normal component(s) the best model
# estimated standard deviation of the normal component(s) in best model
test$best$sigma 
# estimated shape parameter of the exponential components in best model  
test$best$beta
# class indicator inferred using best model chosen. 1 means differential, 0 o.w
bestClass = test$best$class 

# plot best model
DIME.plot.fit(data,test)

# Eg. getting Gaussian mean from iNUDGE model
test$inudge$mu

# Eg. getting Gaussian mean from NUDGE model
test$nudge$mu

# Eg. getting parameters from GNG model
test$gng$mu

# provide initial estimates means of Gaussian in GNG model
test <- DIME(data,gng.max.iter=5,gng.rep=1,inudge.max.iter=5,inudge.rep=1,
 nudge.max.iter=5,nudge.rep=1,gng.K=2,gng.mu=rbind(c(1,2)))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }