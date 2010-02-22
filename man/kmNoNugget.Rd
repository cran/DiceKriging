\name{kmNoNugget}
\alias{kmNoNugget}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Fitting Kriging Models }
\description{
  \code{kmNoNugget} is used to fit kriging models when there is no nugget effect. This function should not be called directly, due to the environments defined in \code{km} to avoid computing twice \code{nxn} matrices. Call \code{km} instead.
}
\usage{
kmNoNugget(model, envir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{  an object of class \code{km}. }
  \item{envir}{  an environment specifying where to assign intermediate values for future gradient calculations. }
}

\value{
  An object of class \code{km}.
}


\references{ Park J-S, Baek J. (2001), Efficient computation of maximum likelihood estimators
                      in a spatial linear model with power exponential covariogram, \emph{Computer Geosciences},
                      \bold{27}, 1-7. }
\author{ O. Roustant, D. Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{ \code{\link{kmNoNugget.init}}, \code{\link{km}}, \code{\link{kmNuggets}}, \code{\link{km1Nugget}}  }

\keyword{models}
\keyword{htest}
\keyword{internal}



