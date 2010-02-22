\name{km1Nugget}
\alias{km1Nugget}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Fitting Kriging Models }
\description{
  \code{km1Nugget} is used to fit kriging models with an unknown nugget effect. This function should not be called directly, due to the environments defined in \code{km} to avoid computing twice \code{nxn} matrices. Call \code{km} instead and specify the item \code{nugget}.
}
\usage{
km1Nugget(model, envir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{  an object of class \code{km}. }
  \item{envir}{  an environment specifying where to assign intermediate values for future gradient calculations. }
 }

\note{
	When a nugget effect is to be estimated, optimization if performed with respect to the correlation parameters and \code{alpha = sigma^2/(sigma^2 +  delta^2)}, with \code{sigma^2} the variance of unnoisy part, and \code{delta^2} the noise variance. Thus, \code{alpha} is the proportion of variance explained by the unnoisy part of the process. It can easily be shown that all other parameters ML estimators are function of the correlation parameters and \code{alpha} estimators.
	}

\value{
  An object of class \code{km}.
}


\author{ O. Roustant, David Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{ \code{\link{km1Nugget.init}}, \code{\link{km}}, \code{\link{kmNoNugget}}, \code{\link{kmNuggets}}  }

\keyword{models}
\keyword{htest}
\keyword{internal}




