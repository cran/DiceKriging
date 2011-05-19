\name{covparam2vect}
\alias{covparam2vect}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Auxiliary function }
\description{  Gather the covariance parameters in a single vector. This is useful in the estimation step. Not for direct use.}
\usage{
covparam2vect(covStruct)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{covStruct}{ an object of class covTensorProduct. See \code{\link{covTensorProduct-class}}}
}

\value{ A vector containing the covariance parameters. }

\author{ O. Roustant, D. Ginsbourger}

\seealso{ \code{\link{vect2covparam}} }

\keyword{models}
\keyword{internal}
