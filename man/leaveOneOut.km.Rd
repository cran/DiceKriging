\name{leaveOneOut.km}
\alias{leaveOneOut.km}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Leave-one-out for a km object }
\description{
  Cross validation by leave-one-out for a \code{km} object.
}
\usage{
leaveOneOut.km(model, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ an object of class "km".  }
  \item{type}{ a character string corresponding to the kriging family, to be chosen between simple kriging ("SK"), or universal kriging ("UK"). At this stage, only "UK" is available.}
}

\details{
  Leave-one-out consists of computing the prediction at a design point when the corresponding observation is removed from the learning set (and this, for all design points).}
  
\section{Warning}{
 Kriging parameters are not re-estimated when removing one observation. With few points, the re-estimated values can be far from those obtained with the entire learning set.
}

\value{ A list composed of:
  \item{mean }{ a vector of length \emph{n}. The ith coordinate is equal to the kriging mean (including the trend) at the ith observation number when removing it from the learning set, }
  \item{sd }{ a vector of length \emph{n}. The ith coordinate is equal to the kriging standard deviation at the ith observation number when removing it from the learning set,}
where \emph{n} is the total number of observations.
}

\references{ 

N.A.C. Cressie (1993), \emph{Statistics for spatial data}, Wiley series in probability and mathematical statistics.

J.D. Martin and T.W. Simpson (2005), Use of kriging models to approximate deterministic computer models, \emph{AIAA Journal}, \bold{43} no. 4, 853-863.

M. Schonlau (1997), \emph{Computer experiments and global optimization}, Ph.D. thesis, University of Waterloo.	
	
}

\author{ O. Roustant, D. Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{ \code{\link{predict,km-method}},  \code{\link{plot,km-method}} }

\keyword{models}

