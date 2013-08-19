\name{scalingFun}
\alias{scalingFun}

\title{
Scaling function 
}
\description{
Parametric transformation of the input space variables. The transformation is obtained coordinatewise by integrating piecewise affine marginal "densities" parametrized by a vector of knots and a matrix of density values at the knots. See references for more detail.
}
\usage{
scalingFun(X, knots, eta, plot=FALSE)
}

\arguments{
  \item{X}{ an n*d matrix standing for a design of n experiments in d-dimensional space }
  \item{knots}{ a list of knots parametrizing the transformation. }
  \item{eta}{ a list of coefficients parametrizing the d marginal transformations. Each element stands for a set of marginal density values at the knots defined above.}
  \item{plot}{ if TRUE plots the image of the columns of X according to the corresponding marginal transformations.}
}

\value{
The image of X by a scaling transformation of parameters knots and eta
}

\references{
Y. Xiong, W. Chen, D. Apley, and X. Ding (2007), \emph{Int. J. Numer. Meth. Engng}, A non-stationary covariance-based Kriging method for metamodelling in engineering design. 
}
\examples{
###############################
# 1D Transform of Xiong et al...
###############################

knots <- list(c(0,0.3,1))
eta <- list(matrix(c(2,0.4,1.2), nrow = 1))
eps <- 0 #0.0001
t <- matrix(seq(0 + eps, 1 - eps, , 399), ncol = 1)
f <- scalingFun(t, knots, eta)

plot(t,f, type="l", main = expression(f[j]), xlab = expression(x[j]), ylab = expression(f[j]))
text(0.7, 0.3, expression(f[j](x[j]) == integral(g[j](t)*dt, 0, x[j] ) ), cex = 2)

###############################
# ...and the corresponding density!
###############################

ScalingDensity1d <- function(t, knots_j, eta_j){
  ind_right <- which.max(knots_j - t >= 0)
  if (ind_right == 1){
    f <- eta_j[ind_right]
  } else {
    f <- eta_j[ind_right - 1] + ((eta_j[ind_right] - eta_j[ind_right - 1])
    /(knots_j[ind_right]-knots_j[ind_right - 1]))*(t-knots_j[ind_right - 1]) 
  }
  return(f)
}

knots <- c(0,0.3,1)
eta <- c(2,0.4,1.2)
t <- seq(0,1,,400)
f <- apply(as.matrix(t), 1, ScalingDensity1d, knots_j = knots, eta_j = eta)

plot(t, f, type = "l", ylim = c(0,2.2), xlab = expression(x[j]), ylab = expression(g[j]),
           main = expression(g[j])) 

eps <- 0.15
text(knots, eta + eps, labels = c(expression(eta[1]^j), expression(eta[2]^j),
                                  expression(eta[3]^j)), cex = 1.8)

mtext(expression(zeta[1]^j), side=1, cex = 1.8, at=0, padj=1.7)
mtext(expression(zeta[2]^j), side=1, cex = 1.8, at=0.3, padj=1.7)
mtext(expression(zeta[3]^j), side=1, cex = 1.8, at=1, padj=1.7)

segments(x0=knots[1], y0=0, x1=knots[1], y1=eta[1], lty=2)
segments(x0=knots[2], y0=0, x1=knots[2], y1=eta[2], lty=2)
segments(x0=knots[3], y0=0, x1=knots[3], y1=eta[3], lty=2)

polygon(x=c(t,rep(1,400),rev(t),rep(0,400)), 
y = c(f, seq(eta[3],0,,400), rep(0,400), seq(0,eta[1],,400)),
density = 15, angle = 45, col="orange", border = NA)
}
\keyword{models}
