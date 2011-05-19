\name{covMatrixDerivative}
\alias{covMatrixDerivative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Covariance matrix derivatives }
\description{
  Computes a partial derivative of the covariance matrix \code{C0} in function \code{\link{covMatrix}}.
}

\usage{
covMatrixDerivative(object, X, C, k, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object specifying the covariance structure. }
  \item{X}{a matrix whose columns represent locations. }
  \item{C}{a matrix corresponding to the covariance matrix for the locations specified in the X argument, when there is no nugget effet nor observation noise. }
  \item{k}{an integer representing the partial derivative index. }
  \item{\dots}{additional parameters, typically an environment used for storage (when the \code{object} type is \code{covAffineScaling})}
}

\value{
  a list with the following items :
  \item{dC}{ a matrix representing the partial derivative of \code{C}.}
  \item{dlogR}{ a matrix representing the partial derivative of logarithm of R = C/sigma2, where sigma2 is the variance of the stationary part of the stochastic process.}
}


\author{ Olivier Roustant, David Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{  \code{\link{covMatrix}} }

\keyword{models}
\keyword{htest}
\keyword{internal}