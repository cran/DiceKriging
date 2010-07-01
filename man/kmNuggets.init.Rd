\name{kmNuggets.init}
\alias{kmNuggets.init}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Fitting Kriging Models}
\description{
  \code{kmNuggets.init} is used to give initial values to fit kriging models, in presence of noisy observations.
}
\usage{
kmNuggets.init(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ an object of class \code{km}. }
}

\value{
  An object of class \code{km}.
}

\details{
  The procedure can be summarized in 4 stages:
  \tabular{ll}{ 
  1) \tab Compute the variogram and give a first estimation of the process variance, as well as lower and upper bounds. \cr
  2) \tab Simulate several values for the process variance, around the estimation obtained at stage 1). The number of simulations is the one given in \code{model@control$pop.size}. \cr
  3) \tab If no initial value is provided for the other covariance parameters, simulate them uniformly inside the domain delimited by \code{model@lower} and \code{model@upper}. The number of simulations is the same as in stage 2). \cr
  4) \tab Compute the likelihood at each simulated "point" (variance + other covariance parameters), and take the best one. This point gives the first initial values.
  }
}


\author{ O. Roustant, David Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{ \code{\link{km}}, \code{\link{kmNuggets}} }

\keyword{models}
\keyword{htest}
\keyword{internal}



