\name{DIME-package}
\alias{DIME_package}
\alias{DIME.package}
\docType{package}
\title{
DIME (Differential Identification using Mixtures Ensemble) 
}
\description{
A robust differential identification method that considers an ensemble of finite 
mixture models combined with a local false discovery rate \emph{(fdr)} for 
analyzing ChIP-seq data comparing two samples. \cr
This package can also be used to identify differential
in other high throughput data such as microarray, methylation etc.\cr
After normalization, an Exponential-Normal(k) or a Uniform-Normal(k) mixture is
fitted to the data. The (k)-normal component can represent either differential
regions or non-differential regions depending on their location and spread. The
exponential or uniform represent differentially sites. local \emph{(fdr)} are
computed from the fitted model.
Unique features of the package: 
\enumerate{\item{Accurate modeling of data that comes from 
any distribution by the use of multiple normal components (any distribution can be
accurately represented by mixture of normal).}
\item{Using ensemble of mixture models
allowing data to be accurately & efficiently represented. Then two-phase selection ensure the 
selection of the best overall model.} 
\item{This method can be used as a general
program to fit a mixture of uniform-normal or uniform-k-normal or exponential-k-normal} 
}
}
\details{
\tabular{ll}{
Package: \tab DIME\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-11-19\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Cenny Taslim \email{taslim.2@osu.edu}, with contributions 
from Abbas Khalili \email{khalili@stat.ubc.ca}, 
Dustin Potter \email{potterdp@gmail.com}, and 
Shili Lin \email{shili@stat.osu.edu}\cr
Maintainer: Cenny Taslim \email{taslim.2@osu.edu} or 
Shili Lin  \email{shili@stat.osu.edu}
}
\references{
\itemize{
\item{Khalili, A., Huang, T., and Lin, S. (2009). A robust unified approach to 
analyzing methylation and gene expression data. \emph{Computational Statistics 
and Data Analysis}, 53(5), 1701-1710.}
\item{Dean, N. and Raftery, A. E. (2005). Normal uniform mixture differential 
gene expression detection for cDNA microarrays. \emph{BMC Bioinformatics}, 6, 173.}
\item{Taslim, C., Wu, J., Yan, P., Singer, G., Parvin, J., Huang, T., Lin, S., 
and Huang, K. (2009). Comparative study on chip-seq data: normalization and 
binding pattern characterization. \emph{Bioinformatics}, 25(18), 2334-2340.}
}
}
\keyword{package}
