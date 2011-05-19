\name{covAffineScaling-class}
\Rdversion{1.1}
\docType{class}
\alias{covAffineScaling-class}
\alias{coefficients,covAffineScaling-method}
\alias{covMat1Mat2,covAffineScaling-method}
\alias{covMatrix,covAffineScaling-method}
\alias{covMatrixDerivative,covAffineScaling-method}
\alias{covParametersBounds,covAffineScaling-method}
\alias{inputnames,covAffineScaling-method}
%\alias{inputnames<-,covAffineScaling,character-method}
\alias{kernelname,covAffineScaling-method}
\alias{ninput,covAffineScaling-method}
\alias{nuggetflag,covAffineScaling-method}
\alias{nuggetvalue,covAffineScaling-method}
\alias{nuggetvalue<-,covAffineScaling,numeric-method}
\alias{show,covAffineScaling-method}
\alias{slice,covAffineScaling,numeric-method}
\alias{summary,covAffineScaling-method}

\title{Class "covAffineScaling" }
\description{ Composition of isotropic kernels with coordinatewise 
non-linear scaling obtained by integrating affine functions }

\section{Objects from the Class}{
In 1-dimension, the covariance kernels are parameterized as in (Rasmussen, Williams, 2006). Denote by \code{theta} the range parameter, \code{p} the exponent parameter (for power-exponential covariance), \code{s} the standard deviation, and \code{h=|x-y|}. Then we have \code{C(x,y) = s^2 * k(x,y)}, with:

\tabular{ll}{
Gauss \tab \code{ k(x,y) = exp(-1/2*(h/theta)^2) } \cr
Exponential \tab \code{ k(x,y) = exp(-h/theta) } \cr
Matern(3/2) \tab \code{ k(x,y) = (1+sqrt(3)*h/theta)*exp(-sqrt(3)*h/theta) } \cr
Matern(5/2) \tab \code{ k(x,y) = (1+sqrt(5)*h/theta+(1/3)*5*(h/theta)^2)} \cr 
\tab \code{          *exp(-sqrt(5)*h/theta)} \cr
}

Here, in every dimension, the corresponding one-dimensional stationary kernel \code{k(x,y)} is replaced by \code{k(f(x),f(y))}, where \code{f} is a continuous monotonic function indexed by 2 parameters (see the references for more detail). 

}

\section{Slots}{
  \describe{
	 \item{\code{d}:}{Object of class \code{"integer"}. The spatial dimension. }
    \item{\code{knots}:}{Object of class \code{"numeric"}. A vector specifying the position of the two knots, common to all dimensions.}
    \item{\code{eta}:}{Object of class \code{"matrix"}. A d*2 matrix of scaling coefficients, parametrizing the coordinatewise transformations in the d dimensions. }
	 \item{\code{name}:}{Object of class \code{"character"}. The covariance function name. To be chosen between \code{"gauss", "matern5_2", "matern3_2", "exp"}, and \code{"powexp"} }
	 \item{\code{paramset.n}:}{Object of class \code{"integer"}. 1 for covariance depending only on the ranges parameters, 2 for "powexp" which also depends on exponent parameters. }
	 \item{\code{var.names}:}{Object of class \code{"character"}. The variable names. }
	 \item{\code{sd2}:}{Object of class \code{"numeric"}. The variance of the stationary part of the process. }
	 \item{\code{known.covparam}:}{Object of class \code{"character"}. Internal use. One of: "None", "All". }
    \item{\code{nugget.flag}:}{Object of class \code{"logical"}. Is there a nugget effect? }
    \item{\code{nugget.estim}:}{Object of class \code{"logical"}. Is the nugget effect estimated or known? }
    \item{\code{nugget}:}{Object of class \code{"numeric"}. If there is a nugget effect, its value (homogeneous to a variance). }
    \item{\code{param.n}:}{Object of class \code{"integer"}. The total number of parameters. }
  }
}

\section{Extends}{
Class \code{"\linkS4class{covKernel}"}, directly.
}
\section{Methods}{
  \describe{
    \item{coefficients}{\code{signature(object = "covAffineScaling")}: ... }
    \item{covMat1Mat2}{\code{signature(object = "covAffineScaling")}: ... }
    \item{covMatrix}{\code{signature(object = "covAffineScaling")}: ... }
    \item{covMatrixDerivative}{\code{signature(object = "covAffineScaling")}: ... }
    \item{covParametersBounds}{\code{signature(object = "covAffineScaling")}: ... }
    \item{inputnames}{\code{signature(x = "covAffineScaling")}: ... }
    \item{kernelname}{\code{signature(x = "covAffineScaling")}: ... }
    \item{ninput}{\code{signature(x = "covAffineScaling")}: ... }
    \item{nuggetflag}{\code{signature(x = "covAffineScaling")}: ... }
    \item{nuggetvalue}{\code{signature(x = "covAffineScaling")}: ... }
    \item{show}{\code{signature(object = "covAffineScaling")}: ... }
    \item{slice}{\code{signature(object = "covAffineScaling", index = "numeric")}: ... }
    \item{summary}{\code{signature(object = "covAffineScaling")}: ... }
	 }
}

\author{ O. Roustant, D. Ginsbourger}

\references{
	
Y. Xiong, W. Chen, D. Apley, and X. Ding (2007), \emph{Int. J. Numer. Meth. Engng}, A non-stationary covariance-based Kriging method for metamodelling in engineering design. 

}


\seealso{
	\code{\linkS4class{km}}
	\code{\linkS4class{covTensorProduct}}
	\code{\linkS4class{covIso}}
	\code{\linkS4class{covKernel}}
}


\examples{
showClass("covAffineScaling")

###############################
# 1D example with automatic scaling 
###############################

f11_xiong <- function (x) {
return( sin(30*(x - 0.9)^4)*cos(2*(x - 0.9)) + (x - 0.9)/2)}

t <- seq(0,1,,300)
f <- f11_xiong(t)

par(mfrow=c(1,1))
plot(t,f,type="l", ylim=c(-1,0.6))

doe <- data.frame(seq(0,1,,20))
resp <- f11_xiong(doe)

m <- km(design=doe, response=resp, scaling=TRUE, gr=TRUE, covtype="matern5_2", coef.var=1, coef.trend=0)

p <- predict(m, t, "UK")$mean
s <- predict(m, t, "UK")$sd

lines(t,p,col="blue")
lines(t,p + 2*s,col="blue")
lines(t,p - 2*s,col="blue")

###############################
# 2D example with automatic scaling 
###############################

# a 16-points factorial design, and the corresponding response
d <- 2; n <- 16
design.fact <- expand.grid(seq(0,1,length=4), seq(0,1,length=4))
design.fact <- data.frame(design.fact); 
names(design.fact)<-c("x1", "x2")
y <- branin(design.fact) 

# graphics 
n.grid <- 50
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
response.grid <- apply(design.grid, 1, branin)

par(mfrow=c(2,1))
knots <- list(seq(from=0, to=1, length=4), seq(from=0, to=1, length=4))
m <- km(design=design.fact, response=y, scaling=TRUE, gr=TRUE)
predicted.values.model <- predict(m, design.grid, "UK")$mean
sd.model <- predict(m, design.grid, "UK")$sd
contour(x.grid, y.grid, matrix(predicted.values.model, n.grid, n.grid), 50, main="Kriging Mean")
contour(x.grid, y.grid, matrix(sd.model, n.grid, n.grid), 50, 
        main="Kriging standard deviation")

###
}
\keyword{classes}
