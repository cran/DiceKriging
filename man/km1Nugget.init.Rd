\name{km1Nugget.init}
\alias{km1Nugget.init}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Fitting Kriging Models}
\description{
  \code{km1Nugget.init} is used to give good initial values to fit kriging models when there is an unknown nugget effect to be estimated.
}
\usage{
km1Nugget.init(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ an object of class \code{km}. }
}

\value{
  An object of class \code{km}.
}

\details{
  The procedure can be summarized in 4 stages : 
  \tabular{ll}{
  1) \tab Compute the variogram and deduce a first estimation of the total variance. If an initial value is provided for \code{nugget}, check its compatibility with the estimated variance. If not, use again the variogram to give a first estimation of the nugget effect. \cr
  2) \tab Simulate several values for the nugget effect and the process variance, around the estimations obtained at stage 1). The number of simulations is the one given in \code{model@control$pop.size}. \cr
  3) \tab If no initial value is provided for the other covariance parameters, simulate them uniformly inside the domain delimited by \code{model@lower} and \code{model@upper}. The number of simulations is the same as in stage 2). \cr
  4) \tab Compute the likelihood at each simulated "point" (variance + nugget effect + other covariance parameters), and take the best one. This point gives the first initial values.
   }
}


\author{ O. Roustant, David Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{ \code{\link{km}}, \code{\link{kmEstimate}} }

\keyword{models}
\keyword{htest}
\keyword{internal}




