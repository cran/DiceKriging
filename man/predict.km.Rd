\name{predict.km}
\alias{predict.km}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Predict values and confidence intervals at newdata for a km object }
\description{
Predicted values and (marginal of joint) conditional variances based on a \code{km} model.  95 \% confidence intervals are given, based on strong assumptions: Gaussian process assumption, specific prior distribution on the trend parameters, known covariance parameters. This might be abusive in particular in the case where estimated covariance parameters are plugged in.
}
\usage{
predict.km(object, newdata, type, se.compute=TRUE, 
                                  cov.compute=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of class \code{km}. }
  \item{newdata}{ a vector, matrix or data frame containing the points where to perform predictions. }
  \item{type}{ a character string corresponding to the kriging family, to be chosen between simple kriging ("SK"), or universal kriging ("UK").}
  \item{se.compute}{ an optional boolean. If \code{FALSE}, only the kriging mean is computed. If \code{TRUE}, the kriging variance (actually, the corresponding standard deviation) and confidence intervals are computed too.}
  \item{cov.compute}{ an optional boolean. If \code{TRUE} the conditional covariance matrix is computed.}
  \item{...}{ no other argument for this method.}
}

\value{
  \item{mean }{ kriging mean (including the trend) computed at \code{newdata}. }
  \item{sd }{ kriging standard deviation computed at \code{newdata}. Not computed if 
              \code{se.compute=FALSE}. }
  \item{cov }{ kriging conditional covariance matrix. Not computed if \code{cov.compute=FALSE} (default).}
  \item{lower95, }{ }
  \item{upper95 }{ bounds of the 95 \% confidence interval computed at \code{newdata} (to be interpreted with special care when parameters are estimated, see description above). Not computed if \code{se.compute=FALSE}.}
  \item{c }{ an auxiliary matrix, containing all the covariances between newdata and the initial design points. }
  \item{Tinv.c }{ an auxiliary vector, equal to \code{T^(-1)*c}.}
}

\details{ 
When \code{type = "UK"}, the estimated variance and covariance are multiplied by \code{n/(n-p)}, where \code{n} and \code{p} are respectively the number of rows and the number of columns of the design matrix F. This would lead to an unbiased estimate if the correlation parameters were known (but this is not the case here...).}

\section{Warning}{
Beware that the only consistency check between \code{newdata} and the experimental design is to test whether they have same number of columns. In that case, the columns of \code{newdata} are interpreted in the same order as the initial design.
}

\references{

N.A.C. Cressie (1993), \emph{Statistics for spatial data}, Wiley series in probability and mathematical statistics.

A.G. Journel and C.J. Huijbregts (1978), \emph{Mining Geostatistics}, Academic Press, London.

D.G. Krige (1951), A statistical approach to some basic mine valuation problems on the witwatersrand, \emph{J. of the Chem., Metal. and Mining Soc. of South Africa}, \bold{52} no. 6, 119-139.

J.D. Martin and T.W. Simpson (2005), Use of kriging models to approximate deterministic computer models, \emph{AIAA Journal}, \bold{43} no. 4, 853-863.

G. Matheron (1963), Principles of geostatistics, \emph{Economic Geology}, \bold{58},
1246-1266.

G. Matheron (1969), Le krigeage universel, \emph{Les Cahiers du Centre de Morphologie Mathematique de Fontainebleau}, \bold{1}.

J.-S. Park and J. Baek (2001), Efficient computation of maximum likelihood estimators in a spatial linear model with power exponential covariogram, \emph{Computer Geosciences}, \bold{27} no. 1, 1-7.

C.E. Rasmussen and C.K.I. Williams (2006), \emph{Gaussian Processes for Machine Learning}, the MIT Press, \url{http://www.GaussianProcess.org/gpml}

J. Sacks, W.J. Welch, T.J. Mitchell, and H.P. Wynn (1989), Design and analysis of computer experiments, \emph{Statistical Science}, \bold{4}, 409-435.

Y. Xiong, W. Chen, D. Apley, and X. Ding (2007), \emph{Int. J. Numer. Meth. Engng}, A non-stationary covariance-based Kriging method for metamodelling in engineering design. 

}

\author{ O. Roustant, D. Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{  \code{\link{km}},  \code{\link{plot.km}} }

\keyword{models}
\keyword{methods}

\examples{
# ------------
# a 1D example
# ------------

x <- c(0, 0.4, 0.6, 0.8, 1);
y <- c(-0.3, 0, 0, 0.5, 0.9)

formula <- y~x   # try also   y~1  and  y~x+I(x^2)

model <- km(formula=formula, design=data.frame(x=x), response=data.frame(y=y), 
            covtype="matern5_2")

tmin <- -0.5; tmax <- 2.5
t <- seq(from=tmin, to=tmax, by=0.005)
color <- list(SK="black", UK="blue")

# Results with Universal Kriging formulae (mean and and 95% intervals)
p.UK <- predict(model, newdata=t, type="UK")
plot(t, p.UK$mean, type="l", ylim=c(min(p.UK$lower95),max(p.UK$upper95)), 
                   xlab="x", ylab="y")
lines(t, p.UK$lower95, col=color$UK, lty=2)
lines(t, p.UK$upper95, col=color$UK, lty=2)
points(x, y, col="red", pch=19)
abline(h=0)

# Results with Simple Kriging (SK) formula. The difference between the width of
# SK and UK intervals are due to the estimation error of the trend parameters 
# (but not to the range parameters, not taken into account in the UK formulae).
p.SK <- predict(model, newdata=t, type="SK")
lines(t, p.SK$mean, type="l", ylim=c(-7,7), xlab="x", ylab="y")
lines(t, p.SK$lower95, col=color$SK, lty=2)
lines(t, p.SK$upper95, col=color$SK, lty=2)
points(x, y, col="red", pch=19)
abline(h=0)

legend.text <- c("Universal Kriging (UK)", "Simple Kriging (SK)")
legend(x=tmin, y=max(p.UK$upper), legend=legend.text, 
       text.col=c(color$UK, color$SK), col=c(color$UK, color$SK), 
       lty=3, bg="white")


# ---------------------------------------------------------------------------------
# a 1D example (following)- COMPARISON with the PREDICTION INTERVALS for REGRESSION
# ---------------------------------------------------------------------------------
# There are two interesting cases: 
# *  When the range parameter is near 0 ; Then the intervals should be nearly 
#    the same for universal kriging as for regression. This is because the 
#    uncertainty around the range parameter is not taken into account in the
#    Universal Kriging formula.
# *  Where the predicted sites are "far" (relatively to the spatial correlation) 
#    from the design points ; in this case, the kriging intervals are not equal 
#    but nearly proportional to the regression ones, since the variance estimate 
#    for regression is not the same than for kriging (that depends on the 
#    range estimate)

x <- c(0, 0.4, 0.6, 0.8, 1);
y <- c(-0.3, 0, 0, 0.5, 0.9)

formula <- y~x   # try also   y~1  and  y~x+I(x^2)
upper <- 0.05 	 # this is to get something near to the regression case. 
                 # Try also upper=1 (or larger) to get usual results.

model <- km(formula=formula, design=data.frame(x=x), response=data.frame(y=y), 
               covtype="matern5_2", upper=upper)

tmin <- -0.5; tmax <- 2.5
t <- seq(from=tmin, to=tmax, by=0.005)
color <- list(SK="black", UK="blue", REG="red")

# Results with Universal Kriging formulae (mean and and 95% intervals)
p.UK <- predict(model, newdata=t, type="UK")
plot(t, p.UK$mean, type="l", ylim=c(min(p.UK$lower95),max(p.UK$upper95)), 
                   xlab="x", ylab="y")
lines(t, p.UK$lower95, col=color$UK, lty=2)
lines(t, p.UK$upper95, col=color$UK, lty=2)
points(x, y, col="red", pch=19)
abline(h=0)

# Results with Simple Kriging (SK) formula. The difference between the width of
# SK and UK intervals are due to the estimation error of the trend parameters 
# (but not to the range parameters, not taken into account in the UK formulae).
p.SK <- predict(model, newdata=t, type="SK")
lines(t, p.SK$mean, type="l", ylim=c(-7,7), xlab="x", ylab="y")
lines(t, p.SK$lower95, col=color$SK, lty=2)
lines(t, p.SK$upper95, col=color$SK, lty=2)
points(x, y, col="red", pch=19)
abline(h=0)

# results with regression given by lm (package stats)
m.REG <- lm(formula)
p.REG <- predict(m.REG, data.frame(x=t), interval="prediction")
lines(t, p.REG[,1], col=color$REG)
lines(t, p.REG[,2], col=color$REG, lty=2)
lines(t, p.REG[,3], col=color$REG, lty=2)

legend.text <- c("Universal Kriging (UK)", "Simple Kriging (SK)", "Regression")
legend(x=tmin, y=max(p.UK$upper), legend=legend.text, 
       text.col=c(color$UK, color$SK, color$REG), 
       col=c(color$UK, color$SK, color$REG), lty=3, bg="white")

# ----------------------------------
# A 2D example - Branin-Hoo function
# ----------------------------------

# a 16-points factorial design, and the corresponding response
d <- 2; n <- 16
fact.design <- expand.grid(seq(0,1,length=4), seq(0,1,length=4))
fact.design <- data.frame(fact.design); names(fact.design)<-c("x1", "x2")
branin.resp <- data.frame(branin(fact.design)); names(branin.resp) <- "y" 

# kriging model 1 : gaussian covariance structure, no trend, 
#                   no nugget effect
m1 <- km(~1, design=fact.design, response=branin.resp, covtype="gauss")

# predicting at testdata points
testdata <- expand.grid(s <- seq(0,1, length=15), s)
predicted.values.model1 <- predict(m1, testdata, "UK")

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