\name{huber}
\alias{huber}
\title{ Huber's weight function }
\description{
 A weight functions used to downweigh outliers.
}
\usage{
huber(input, co, shape = c("full", "lower", "upper"))
}
\arguments{
  \item{input}{ an \strong{R list} of vector of normalized mean 
    (log intensities). }
  \item{co}{ cutoff used in determining weights.  }
  \item{shape}{ parameter determining which outliers are weighted: 
    "full"  - both values above and below -threshold are downweighted;\\
    "lower" - only values below threshold are downweighted; \\
    "upper" - only values above threshold are downweighted.}
}
\value{
 a vector of weights.	
}
\references{ Huber, P. J. (1981) Robust Statistics. John Wiley & Sons }
\author{ Dustin Potter }