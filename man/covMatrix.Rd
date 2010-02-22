\name{covMatrix} 
\alias{covMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Covariance matrix}
\description{
  Computes the covariance matrix at a set of locations for a spatial random process with a given covariance structure.
}

\usage{
covMatrix(X, covStruct, noise.var = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{X}{ a matrix whose columns represent locations. }
 \item{covStruct}{an object of class \code{\link{covTensorProduct-class}}, specifying the covariance structure.}
 \item{noise.var}{for noisy observations : an optional vector containing the noise variance at each observation}
}

\value{
  a list with the following items :
  \item{C0}{ a matrix representing the covariance matrix for the locations specified in the X argument, when there is no nugget effect nor observation noise.}
  \item{C}{ a matrix equal to the covariance matrix (including a possible nugget or observation noise).}
  }


\author{ Olivier Roustant, David Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{  \code{\link{covMatrixDerivative}} }

\keyword{models}