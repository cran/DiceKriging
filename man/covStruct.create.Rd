\name{covStruct.create}
\alias{covStruct.create}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Spatial covariance - Class constructor }

\description{  Creates a covariance structure. }

\usage{
covStruct.create(covtype, d, coef.cov=NULL, coef.var=NULL, nugget = NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covtype}{ a character string specifying the covariance structure. See \code{\link{covTensorProduct-class}}.}
  \item{d}{ an integer containing the spatial dimension.}
  \item{coef.cov}{ an optional vector containing the values for covariance parameters.  }
  \item{coef.var}{ an optional number containing the variance value. }
  \item{nugget}{ an optional variance value standing for the homogenous nugget effect. Default is NULL.}
}

\value{ A formal S4 class of type "covTensorProduct".}

\author{ O. Roustant, D. Ginsbourger}

\seealso{ \code{\link{covTensorProduct-class}} }

\keyword{models}
\keyword{internal}
