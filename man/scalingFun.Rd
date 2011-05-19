\name{scalingFun}
\alias{scalingFun}

\title{
Scaling function 
}
\description{
Parametric transformation of the input space variables. The transformation is obtained coordinatewise by integrating piecewise affine marginal "densities" parametrized by a vector of knots and a matrix of density values at the knots. See references for more detail.
}
\usage{
scalingFun(X, knots, eta)
}

\arguments{
  \item{X}{ an n*d matrix standing for a design of n experiments in d-dimensional space }
  \item{knots}{ a (K+1) vector of knots parametrizing the transformation. The knots are here the same in all dimensions. }
  \item{eta}{ a d*(K+1) matrix of coefficients parametrizing the d marginal transformations. Each line stands for a set of (K+1) marginal density values at the knots defined above.}
}

\value{
The image of X by a scaling transformation of parameters knots and eta
}

\references{
%N.A.C. Cressie (1993), \emph{Statistics for spatial data}, Wiley %series in probability and mathematical statistics.

Y. Xiong, W. Chen, D. Apley, and X. Ding (2007), \emph{Int. J. Numer. Meth. Engng}, A non-stationary covariance-based Kriging method for metamodelling in engineering design. 
}

\examples{

###############################
# A 1D Transform a la Xiong et al...
###############################

knots <- c(0,0.5,1)
eta <- matrix(c(2,0.4,1.2), nrow=1)
t <- matrix(seq(0,1,,400), ncol=1)
f <- scalingFun(t,knots,eta)

plot(t,f, type="l", ylim=c(0,1))

###############################
# ...and the corresponding density!
###############################

ScalingDensity1d <- function(t, knots_j, eta_j){
ind_right <- which.max(knots_j-t>=0)
if(ind_right==1){f <- eta_j[ind_right]}
else{
f <- eta_j[ind_right-1] + ((eta_j[ind_right]-eta_j[ind_right-1])/(knots_j[ind_right]-knots_j[ind_right-1]))*(t-knots_j[ind_right-1]) 
}
return(f)
}

knots <- c(0,0.5,1)
eta <- c(2,0.4,1.2)
t <- seq(0,1,,400)
f <- apply(as.matrix(t),1,ScalingDensity1d, knots_j=knots, eta_j=eta)

plot(t,f, type="l", ylim=c(0,2.2))
text(knots,0, labels=c(expression(zeta[1]),expression(zeta[2]),expression(zeta[3])),  cex = 1.6 )
text(knots,eta, labels=c(expression(eta[1]),expression(eta[2]),expression(eta[3])),  cex = 1.6 )

}
\keyword{models}
