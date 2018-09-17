require("DiceKriging")
require("testthat")




context("Checking scalingFun.dx analytical calculation equals numerical result")

f = function(x) scalingFun1d(x,c(0,.5,1),c(10.5,.5,1.75))
dfdx_a = function(x) scalingFun1d.dx(x,c(0,.5,1),c(10.5,.5,1.75))
dfdx_num = function(x) numDeriv::grad(function(xx) scalingFun1d(matrix(xx,nrow=1), c(0,.5,1),c(10.5,.5,1.75)),x)

for (x in runif(10)) {
  dfdx_ax = dfdx_a(x)
  dfdx_numx = dfdx_num(x)
  test_that(desc=paste0("numerical derivative equals analytical at x=",x," : ",dfdx_ax,"!=", dfdx_numx),expect_true(abs(dfdx_ax - dfdx_numx) < 1E-3))
}




context("Checking km covariance derivative 'covVector.dx'")

f = function(x) {1-1/2*(sin(12*x)/(1+x)+2*cos(7*x)*x^5+0.7)}

X <- matrix(c(0,.25,.5,.75,1.0),ncol=1)
y <- f(X)

set.seed(123)
m <- km(formula=~1, design=X, response=y,scaling=F,control=list(trace=FALSE))

# plot covariance function (of x)
x = seq(0,1,,101) #101 because we need to have X values in x to check derivative is null
cx = matrix(NaN,nrow=nrow(X),ncol=length(x)) # covariance function (in 5d space)
cdx = matrix(NaN,nrow=nrow(X),ncol=length(x)) # derivative of covariance function
for (i in 1:length(x)) {
  mx = predict(m,x[i],type="SK",se.compute=T)
  cx[,i] = covMatrix(m@covariance,rbind(X=m@X,x[i]))$C[1:length(X),1+length(X)]
  cdx[,i] = covVector.dx(m@covariance,x=x[i],X=m@X,mx$c)
}

par(mfrow=c(2,3))
for (j in 1:length(X)) {
  plot(x, cx[j,], type='l',main=paste0('Covariance to X_',j))
  abline(v=X[j,],col='blue')
  for (i in 1:(length(x)/10))
    arrows(x[10*i], cx[j,10*i], x[10*i]+0.1, cx[j,10*i] + cdx[j,10*i]/10,length=0.05,col='red')

  ij = which(x==X[j,1])
  test_that(desc="zero derivative at X_i",expect_true(cdx[j,ij] == 0))
}
par(mfrow=c(1,1))








context("Checking km _scaling_ (1 node) covariance derivative 'covVector.dx'")

library(numDeriv)

f = function(x) {1-1/2*(sin(12*x)/(1+x)+2*cos(7*x)*x^5+0.7)}

X <- matrix(c(0,.25,.5,.75,1.0),ncol=1)
y <- f(X)

set.seed(123)
m_scaling0 <- km(formula=~1, design=X, response=y,scaling=T,knots=list(design=c(.5)),control=list(trace=FALSE))

# plot covariance function (of x)
x = seq(0,1,,101) #101 because we need to have X values in x to check derivative is null
cx = matrix(NaN,nrow=nrow(X),ncol=length(x)) # covariance function (in 5d space)
cdx = matrix(NaN,nrow=nrow(X),ncol=length(x)) # derivative of covariance function
for (i in 1:length(x)) {
  mx = predict(m_scaling0,x[i],type="SK",se.compute=T)
  cx[,i] = covMatrix(m_scaling0@covariance,rbind(X=m_scaling0@X,x[i]))$C[1:length(X),1+length(X)]
  cdx[,i] = covVector.dx(m_scaling0@covariance,x=x[i],X=m_scaling0@X,mx$c)
}

par(mfrow=c(2,3))
for (j in 1:length(X)) {
  plot(x, cx[j,], type='l',main=paste0('Covariance to X_',j))
  abline(v=X[j,],col='blue')
  for (i in 1:(length(x)/10))
    arrows(x[10*i], cx[j,10*i], x[10*i]+0.1, cx[j,10*i] + cdx[j,10*i]/10,length=0.05,col='red')

  ij = which(x==X[j,1])
  test_that(desc="zero derivative at X_i",expect_true(cdx[j,ij] == 0))
}
par(mfrow=c(1,1))









context("Checking km _scaling_ covariance derivative 'covVector.dx'")

f = function(x) {x <- x^2; 1-1/2*(sin(12*x)/(1+x)+2*cos(7*x)*x^5+0.7)}

X <- matrix(c(0,0.125,.25,0.375,.5,0.625,.75,0.875,1.0),ncol=1)
y <- f(X)

set.seed(123)
m_scaling <- km(formula=~1, design=X, response=y,scaling=T,knots=list(design=c(0,.5,1)),control=list(trace=FALSE))

# plot covariance function (of x)
x = seq(0,1,,201) #101 because we need to have X values in x to check derivative is null
cx = matrix(NaN,nrow=nrow(X),ncol=length(x)) # covariance function (in 5d space)
cdx = matrix(NaN,nrow=nrow(X),ncol=length(x)) # derivative of covariance function
for (i in 1:length(x)) {
  mx = predict(m_scaling,x[i],type="SK",se.compute=T)
  cx[,i] = covMatrix(m_scaling@covariance,rbind(X=m_scaling@X,x[i]))$C[1:length(X),1+length(X)]
  cdx[,i] = covVector.dx(m_scaling@covariance,x=x[i],X=m_scaling@X,mx$c)
}

par(mfrow=c(3,3))
for (j in 1:length(X)) {
  plot(x, cx[j,], type='l',main=paste0('Covariance to X_',j))
  abline(v=X[j,],col='blue')
  for (i in 1:(length(x)/10))
    arrows(x[10*i], cx[j,10*i], x[10*i]+0.1, cx[j,10*i] + cdx[j,10*i]/10,length=0.05,col='red')

  ij = which(x==X[j,1])
  test_that(desc="zero derivative at X_i",expect_true(cdx[j,ij] == 0))
}
par(mfrow=c(1,1))





context("Checking km _affine scaling_ (no node given) covariance derivative 'covVector.dx'")

library(numDeriv)

f = function(x) {x <- x^2; 1-1/2*(sin(12*x)/(1+x)+2*cos(7*x)*x^5+0.7)}

X <- matrix(c(0,0.125,.25,0.375,.5,0.625,.75,0.875,1.0),ncol=1)
y <- f(X)

set.seed(123)
m_ascaling <- km(formula=~1, design=X, response=y,scaling=T,knots=NULL,control=list(trace=FALSE))

# plot covariance function (of x)
x = seq(0,1,,201) #101 because we need to have X values in x to check derivative is null
cx = matrix(NaN,nrow=nrow(X),ncol=length(x)) # covariance function (in 5d space)
cdx = matrix(NaN,nrow=nrow(X),ncol=length(x)) # derivative of covariance function
for (i in 1:length(x)) {
  mx = predict(m_ascaling,x[i],type="SK",se.compute=T)
  cx[,i] = covMatrix(m_ascaling@covariance,rbind(X=m_ascaling@X,x[i]))$C[1:length(X),1+length(X)]
  cdx[,i] = covVector.dx(m_ascaling@covariance,x=x[i],X=m_ascaling@X,mx$c)
}

par(mfrow=c(3,3))
for (j in 1:length(X)) {
  plot(x, cx[j,], type='l',main=paste0('Covariance to X_',j))
  abline(v=X[j,],col='blue')
  for (i in 1:(length(x)/10))
    arrows(x[10*i], cx[j,10*i], x[10*i]+0.1, cx[j,10*i] + cdx[j,10*i]/10,length=0.05,col='red')

  #ij = which(x==X[j,1])
  #test_that(desc="zero derivative at X_i",expect_true(cdx[j,ij] == 0))
}
par(mfrow=c(1,1))






