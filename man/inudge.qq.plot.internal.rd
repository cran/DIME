\name{inudge.qq.plot.internal}
\alias{inudge.qq.plot.internal}
\title{
Internal function for QQ-plot of iNUDGE model
}
\description{
  Internal function needed to generate a QQ-plot for iNUDGE model. 
  \code{\link{inudge.plot.qq}}
}
\usage{
inudge.qq.plot.internal(i, obj)
}
\arguments{
\item{i}{
a number which is used to generate random sample of inudge model.
}
  \item{obj}{
a list object returned by \code{\link{inudge.fit}} function.
}
}
\author{
Cenny Taslim \email{taslim.2@osu.edu}, with contributions from Abbas Khalili
\email{khalili@stat.ubc.ca}, Dustin Potter \email{potterdp@gmail.com}, 
and Shili Lin \email{shili@stat.osu.edu}
}

\seealso{
\code{\link{inudge.fit}}, \code{\link{inudge.plot.qq}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}