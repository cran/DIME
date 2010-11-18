\name{gng.qq.plot.internal}
\alias{gng.qq.plot.internal}
\title{
Internal function for QQ-plot of GNG model.
}
\description{
  Internal function needed to generate a QQ-plot for GNG model. 
  \code{\link{gng.plot.qq}}
}
\usage{
gng.qq.plot.internal(i, obj)
}
\arguments{
\item{i}{
a number which is used to generate random sample of GNG model.
}
  \item{obj}{
a list object returned by \code{\link{gng.fit}} function.
}
}
\author{
Cenny Taslim \email{taslim.2@osu.edu}, with contributions from Abbas Khalili
\email{khalili@stat.ubc.ca}, Dustin Potter \email{potterdp@gmail.com}, 
and Shili Lin \email{shili@stat.osu.edu}
}

\seealso{
\code{\link{gng.fit}}, \code{\link{gng.plot.qq}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}