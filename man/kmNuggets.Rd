\name{kmNuggets}
\alias{kmNuggets}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Fitting Kriging Models }
\description{
  \code{kmNuggets} is used to fit kriging models for noisy observations. This function should not be called directly, due to the environments defined in \code{km} to avoid computing twice \code{nxn} matrices. Call \code{km} instead and specify the item \code{noise.var}.
}
\usage{
kmNuggets(model, envir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{  an object of class \code{km}. }
  \item{envir}{  an environment specifying where to assign intermediate values for future gradient calculations. }
 }

\value{
  An object of class \code{km}.
}


\author{ O. Roustant, David Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{ \code{\link{kmNuggets.init}}, \code{\link{km}}, \code{\link{kmNoNugget}}, \code{\link{km1Nugget}}  }

\keyword{models}
\keyword{htest}
\keyword{internal}




