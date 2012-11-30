\name{leaveOneOutGrad}
\alias{leaveOneOutGrad}

\title{
Leave-one-out least square criterion - Analytical gradient
}

\description{
Returns the analytical gradient of \code{\link{leaveOneOutFun}}.
}

\usage{
leaveOneOutGrad(param, model, envir)
}
\arguments{
  \item{param}{ a vector containing the optimization variables.}
  \item{model}{ an object of class \code{km}.}
  \item{envir}{ an environment specifying where to get intermediate values calculated in \code{leaveOneOutFun}.}
}

\value{
the gradient of \code{leaveOneOutFun} at \code{param}.
}

\references{

F. Bachoc (+2012), Cross validation and maximum likelihood estimations of hyper-parameters of Gaussian processes with model misspecification. Submitted.

O. Dubrule (1983), Cross validation of Kriging in a unique neighborhood. \emph{Mathematical Geology}, \bold{15}, 687-699.

}

\author{
O. Roustant, Ecole des Mines de St-Etienne
}

\seealso{
\code{\link{leaveOneOutFun}}
}

\keyword{models}
\keyword{htest}
