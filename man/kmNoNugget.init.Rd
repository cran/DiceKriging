\name{kmNoNugget.init}
\alias{kmNoNugget.init}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Fitting Kriging Models}
\description{
  \code{kmNoNugget.init} is used to give initial values to fit kriging models when there is no nugget effect nor noisy observations.
}

\usage{
kmNoNugget.init(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{  an object of class \code{km}. }
 }

\value{
  An object of class \code{km}.
}

\details{
  The procedure can be summarized in 2 stages:
  \tabular{ll}{ 
   1) \tab If no initial value is provided by the user for the covariance parameters, simulate them uniformly inside the domain delimited by \code{model@lower} and \code{model@upper}. The number of simulations is the one given in \code{model@control$pop.size}. \cr
   2) \tab Compute the likelihood for each parameters set, and select the one that gives the highest value.
   	}
}


\author{ O. Roustant, David Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{ \code{\link{km}}, \code{\link{kmNoNugget}} }

\keyword{models}
\keyword{htest}
\keyword{internal}




