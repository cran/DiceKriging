# ----------------
# CLASS definition
# ----------------

# km Class

setClass("km", 		
	representation( 
		d = "integer",						   ## spatial dimension
		n = "integer",						   ## observations number
			## data
		X = "matrix",				      	   ## the design of experiments, size nxd
		y = "matrix", 					      ## the observations, size nx1
			## trend information
		p = "integer",						   ## 1+number of trend basis functions
		F = "matrix",						   ## the experimental matrix, size nxp
		trend.formula = "formula",		   ## trend form
		trend.coef = "numeric",			   ## trend coefficients, size px1
			## covariance
		covariance = "covKernel",  ## covariance structure (new S4 class, see covStruct.R)
			## noisy observations
		noise.flag = "logical",  	   	   ## Are observations noisy ? 
		noise.var = "numeric",		      ## vector of length n
			## model information
		known.param = "character",		   ## known parameters: "None", "All" or "Trend"
		case = "character",				   ## "NoNugget" : deterministic observations, no nugget effect
												   ## "1Nugget"  : homogenous nugget effect, to be estimated
												   ## "Nuggets"  : known nugget or (exclusive) noisy observations 
			## optimisation 
		param.estim = "logical",			   ## Are parameter estimated ??
		method = "character", 			   ## Statistical criterion : "MLE" or "PMLE"
		penalty = "list",		   			   ## fun ("SCAD"), fun.derivative, value
		optim.method = "character",		   ## Optimisation algorithm : "BFGS" or "gen"
		lower = "numeric",					   ##
		upper = "numeric",					   ## boudaries for parameter estimation. Length covariance@param.n
		control = "list",					   ## pop.size, wait.generations, max.generations, BFGSburnin 
		gr = "logical",						   ## is analytical gradient to be used ?
		call = "language", 				   ## user call
		parinit = "numeric", 				   ## initial values used (given by user or computed)
		logLik = "numeric",				   ## objective function value (concentrated log-Likelihood)
			## auxilliary variables 
		T = "matrix",                     ## Upper triang. factor of the Choleski dec. 
												   ## of the cov. matrix : t(T)%*%T = C. Size nxn
		z = "numeric",						   ## t(T)^(-1)*(y - F%*%beta). Size nx1
		M = "matrix"							   ## t(T)^(-1)%*%F. Size nxp
	), 
	validity = function(object) {
		if (object@n <= object@d) stop("the number of experiments must be larger than the spatial dimension")
		
		if (ncol(object@y) != 1) stop("the response must be 1-dimensional")
		
		if (!identical(nrow(object@X), nrow(object@y))) stop("the number of observations is not equal to the number of experiments")
		
		if (object@covariance@nugget.flag & object@noise.flag) stop("'nugget' and 'noise' cannot be specified together")
		
		if (object@noise.flag) {	
			if (!identical(length(object@noise.var), object@n)) stop("the length of the vector 'noise.var' must be equal to the number of experiments")
		}
		
		if (!identical(object@trend.coef, numeric(0))) {
			if (!identical(length(object@trend.coef), object@p)) stop("the number of trend coefficients is not compatible with the trend formula")
		}
	}
)



if(!isGeneric("show")) {
  setGeneric(name    = "show",
     def = function(object) standardGeneric("show")
  )
}

setMethod("show", "km", 
	function(object){
		show.km(object)		
	}
)

# **********************************************
# 						P L O T  METHOD
# **********************************************


plot.km <- function(x, kriging.type="UK", ...) {
		
	model <- x
	pred <- leaveOneOut.km(model, kriging.type)
	y <- as.matrix(model@y)
	yhat <- pred$mean
	sigma <- pred$sd
	
	resid <- (y-yhat)/sigma
		#par(ask=TRUE)
	xmin <- min(min(yhat), min(y))
	xmax <- max(max(yhat), max(y))
	
	par(mfrow=c(3,1))
	plot(y[,1], yhat, xlim=c(xmin, xmax), ylim=c(xmin, xmax),
		 xlab="Exact values", ylab="Fitted values", main="Leave-one-out", ...)
	lines(c(xmin,xmax), c(xmin, xmax))
	plot(resid, xlab="Index", ylab="Standardized residuals", main="Standardized residuals", ...)
	qqnorm(resid, main="Normal QQ-plot of standardized residuals") 
	qqline(resid)
  par(mfrow=c(1,1))
	
	invisible(pred)
}

if(!isGeneric("plot")) {
	setGeneric(name = "plot",
		def = function(x, y, ...) standardGeneric("plot")
	)
}

setMethod("plot", 
          signature(x="km"), 
	        function(x, y=NULL, kriging.type="UK", ...){
		        plot.km(x=x, kriging.type=kriging.type, ...)
	        }
)



# **********************************************
# 						P R E D I C T  METHOD
# **********************************************


predict.km <- function(object, newdata, type, se.compute=TRUE, cov.compute=FALSE, checkNames=TRUE, ...) {
	# newdata : n x d
	
	X <- object@X
	y <- object@y
	
	if (checkNames) {
		newdata <- checkNames(X1=X, X2=newdata, X1.name="the design", X2.name="newdata")
	} else {
    newdata <- as.matrix(newdata)
  	d.newdata <- ncol(newdata)
	  if (!identical(d.newdata, object@d)) stop("newdata must have the same numbers of columns than the experimental design")
    if (!identical(colnames(newdata), colnames(X))) {
    #  warning("column names mismatch between 'newdata' and the experimental design - the columns of 'newdata' are interpreted in the same order as the experimental design names")
    colnames(newdata) <- colnames(X)
    }
	}
	
	T <- object@T
	z <- object@z
	M <- object@M
		
	beta <- object@trend.coef
		
#	F <- model.matrix(object@trend.formula, data=data.frame(X))
#	aux <- data.frame(newdata,runif(m))
#	names(aux) <- names(data.frame(X,y))
#	F.newdata <- model.matrix(object@trend.formula, data=aux)

	F.newdata <- model.matrix(object@trend.formula, data=data.frame(newdata))
	y.predict.trend <- F.newdata%*%beta
	
	c.newdata <- covMat1Mat2(object@covariance, X1=X, X2=newdata, nugget.flag=object@covariance@nugget.flag)    # compute c(x) for x = newdata ; remark that for prediction (or filtering), cov(Yi, Yj)=0  even if Yi and Yj are the outputs related to the equal points xi and xj.
	
	Tinv.c.newdata <- backsolve(t(T), c.newdata, upper.tri=FALSE)
	y.predict.complement <- t(Tinv.c.newdata)%*%z
	y.predict <- y.predict.trend + y.predict.complement
	y.predict <- as.numeric(y.predict)
	
	output.list <- list(mean = y.predict, c=c.newdata, Tinv.c=Tinv.c.newdata)
	
	if (se.compute==TRUE) {		
		
		s2.predict.1 <- apply(Tinv.c.newdata, 2, crossprod)         # compute c(x)'*C^(-1)*c(x)   for x = newdata
		
		# A FAIRE : 
		# REMPLACER total.sd2 par cov(Z(x),Z(x)) ou x = newdata
		# partout dans les formules ci-dessous
		# c'est utile dans le cas non stationnaire
		
		
		if (object@covariance@nugget.flag) {
			total.sd2 <- object@covariance@sd2 + object@covariance@nugget
		} else total.sd2 <- object@covariance@sd2

		
		if (type=="SK") {
			s2.predict <- pmax(total.sd2 - s2.predict.1, 0)
			s2.predict <- as.numeric(s2.predict)
			q95 <- qnorm(0.975)
		}
		else if (type=="UK") {
			T.M <- chol(t(M)%*%M)   # equivalently : qrR <- qr.R(qr(M))
			s2.predict.mat <- backsolve(t(T.M), t(F.newdata - t(Tinv.c.newdata)%*%M) , upper.tri=FALSE)

			s2.predict.2 <- apply(s2.predict.mat, 2, crossprod)
			s2.predict <- pmax(total.sd2 - s2.predict.1 + s2.predict.2, 0)
			s2.predict <- as.numeric(s2.predict)
#			n <- nrow(X)
#			p <- ncol(F)
			s2.predict <- s2.predict * object@n/(object@n - object@p)
			q95 <- qt(0.975, object@n - object@p)
		}
	
		lower95 <- y.predict - q95*sqrt(s2.predict)
		upper95 <- y.predict + q95*sqrt(s2.predict)
		
#		output.list <- c(output.list, sd = sqrt(s2.predict), lower95 = lower95, upper95=upper95)
		output.list$sd <- sqrt(s2.predict)
		output.list$lower95 <- lower95
		output.list$upper95 <- upper95
	}
	
	if (cov.compute==TRUE) {		
			
		if (object@covariance@nugget.flag) {
			total.sd2 <- object@covariance@sd2 + object@covariance@nugget
		} else total.sd2 <- object@covariance@sd2
		
		C.newdata <- covMatrix(object@covariance, newdata)[[1]]
		cond.cov <- C.newdata - crossprod(Tinv.c.newdata)
		
		if (type=="UK") {	
			T.M <- chol(t(M)%*%M)   # equivalently : qrR <- qr.R(qr(M))
			s2.predict.mat <- backsolve(t(T.M), t(F.newdata - t(Tinv.c.newdata)%*%M), upper.tri=FALSE)
			cond.cov <- cond.cov + crossprod(s2.predict.mat)
			cond.cov <- cond.cov * object@n/(object@n - object@p)
		}
	
		output.list$cov <- cond.cov
	
	}
	
	return(output.list)
			
}



if(!isGeneric("predict")) {
	setGeneric(name = "predict",
		def = function(object, ...) standardGeneric("predict")
	)
}

setMethod("predict", "km", 
	function(object, newdata, type, se.compute=TRUE, cov.compute=FALSE, checkNames=TRUE, ...) {
		predict.km(object=object, newdata=newdata, type=type, se.compute=se.compute, 
               cov.compute=cov.compute, checkNames=checkNames, ...)
	}
)


# **********************************************
# 			     S I M U L A T E  METHOD
# **********************************************

simulate.km <- function(object, nsim=1, seed=NULL, newdata=NULL, cond=FALSE, nugget.sim=0, checkNames=TRUE, ...) {
	
	if (!is.numeric(nugget.sim)) stop("'nugget.sim' must be a number")
	if (nugget.sim<0) stop("nugget.sim (homogenous to a variance) must not be negative")
	if (!is.logical(cond)) stop("'cond' must be TRUE/FALSE")
	if ((!is.null(newdata)) & (checkNames)) newdata <- checkNames(X1=object@X, X2=newdata, X1.name="the design", X2.name="newdata")
    
	#if (object@noise.flag) {
	#	if (!is.null(newdata)) warning("for heteroscedastic noise, simulation is possible only at design points: 'newdata' is ignored")
	#	newdata <- object@X		
	#	F.newdata <- object@F
	#	T.newdata <- object@T
	#} else {
	if (is.null(newdata)) {
		newdata <- object@X
		F.newdata <- object@F
		T.newdata <- object@T
	} else {
		newdata <- as.matrix(newdata)
		m <- nrow(newdata)
		if (!identical(ncol(newdata), object@d)) 
   		stop("newdata must have the same numbers of columns than the experimental design")
   		if (!identical(colnames(newdata), colnames(object@X))) {
         colnames(newdata) <- colnames(object@X)
   		}
		  F.newdata <- model.matrix(object@trend.formula, data = data.frame(newdata))
		  Sigma <- covMatrix(object@covariance, newdata)[[1]]
		  T.newdata <- chol(Sigma + diag(nugget.sim, m, m))
	  }
  #}
		
	y.trend <- F.newdata %*% object@trend.coef
	
	m <- nrow(newdata)
	
	if (!cond) {			# non conditional simulations
		white.noise <- matrix(rnorm(m*nsim), m, nsim)
		y.rand <- t(T.newdata) %*% white.noise
		y <- matrix(y.trend, m, nsim) + y.rand
	} else {				# simulations conditional to the observations
		#if (object@noise.flag) {
		#	stop("conditional simulations not available for heterogeneous observations")
		#} else {
		Sigma21 <- covMat1Mat2(object@covariance, X1=object@X, X2=newdata, nugget.flag=FALSE)           # size n x m
		Tinv.Sigma21 <- backsolve(t(object@T), Sigma21, upper.tri = FALSE)     # t(T22)^(-1) * Sigma21,  size  n x m
		y.trend.cond <- y.trend + t(Tinv.Sigma21) %*% object@z                 # size m x 1
			
		if (!is.null(newdata)) {
			Sigma11 <- Sigma
		} else Sigma11 <- t(object@T) %*% object@T	
			
		Sigma.cond <- Sigma11 - t(Tinv.Sigma21) %*% Tinv.Sigma21          # size m x m
		T.cond <- chol(Sigma.cond + diag(nugget.sim, m, m))			
		white.noise <- matrix(rnorm(m*nsim), m, nsim)
		y.rand.cond <- t(T.cond) %*% white.noise
		y <- matrix(y.trend.cond, m, nsim) + y.rand.cond	
  }
	
	return(t(y))

}


if(!isGeneric("simulate")) {
	setGeneric(name = "simulate",
		def = function(object, nsim=1, seed=NULL, ...) standardGeneric("simulate")
	)
}

setMethod("simulate", "km", 
	function(object, nsim=1, seed=NULL, newdata=NULL, cond=FALSE, nugget.sim=0, checkNames=TRUE, ...) {
		simulate.km(object=object, nsim=nsim, newdata=newdata, cond=cond, 
                nugget.sim=nugget.sim, checkNames=checkNames, ...)
	}
)


