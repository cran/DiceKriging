\name{vect2covparam}
\alias{vect2covparam}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Auxiliary function }
\description{  Extract the covariance parameters from a single vector. Not for direct use.}

\usage{
vect2covparam(param, covStruct)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{param}{ a vector containing the covariance parameters.}
   \item{covStruct}{ an object of class \code{\link{covTensorProduct-class}}. }
}

\value{ An object of class covTensorProduct }

\author{ O. Roustant, D. Ginsbourger}

\seealso{ \code{\link{covparam2vect}} }

\keyword{models}
\keyword{internal}
