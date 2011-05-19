`km` <-
function(formula=~1, design, response, covtype="matern5_2",  
			 coef.trend=NULL, coef.cov=NULL, coef.var=NULL,
         nugget=NULL, nugget.estim=FALSE, noise.var=NULL, penalty=NULL, 
         optim.method="BFGS", lower=NULL, upper=NULL, parinit=NULL, control=NULL, gr=TRUE, iso=FALSE, scaling=FALSE) {
	
	model <- new("km")
	
	model@call <- match.call()
	
	#if (!is.data.frame(design)) stop("Argument design must be a data.frame. This is to avoid variables names misspecification.")

	# formula : remove automatically the response from it
	data = data.frame(design)
	model@trend.formula <- formula <- drop.response(formula, data=data)
	F <- model.matrix(formula, data=data)
	
	X <- as.matrix(design)
	y <- as.matrix(response)
	model@X <- X
	model@y <- y
	model@d <- ncol(X)
	model@n <- nrow(X)
	model@F <- F
	model@p <- ncol(F)
	model@noise.flag <- (length(noise.var)!=0)
	model@noise.var <- as.numeric(noise.var)

	known.param <- ((length(coef.trend)!=0) & (length(coef.cov)!=0) & (length(coef.var)!=0))
	if (known.param) {
		nugget.estim <- FALSE
		known.covparam <- "All"
	} else {
		coef.var <- coef.cov <- NULL
		known.covparam <- "None"
	}
	
	model@covariance <- covStruct.create(covtype=covtype, d=model@d, known.covparam=known.covparam, coef.cov=coef.cov, coef.var=coef.var, nugget=nugget, nugget.estim=nugget.estim, nugget.flag=((length(nugget)!=0) | nugget.estim),  iso=iso, scaling=scaling, var.names=colnames(X))

	# Now, at least some parameters are unknown
	
	if (known.param) {
		model@trend.coef <- as.numeric(coef.trend)
		model@param.estim <- FALSE
		model@known.param <- "All"	
#		validObject(model, complete=TRUE)
		model <- computeAuxVariables(model)
		return(model)
	}		
	 
	if ((length(coef.trend)!=0) & (length(coef.cov)==0) & (length(coef.var)==0)) {
		model@trend.coef <- as.numeric(coef.trend)
		model@known.param <- "Trend"
	} else {
		model@known.param <- "None"
	}
	
	if (length(penalty)==0) {
		model@method <- method <- "MLE"}
	else {
		if (covtype!="gauss") stop("At this stage, Penalized Maximum Likelihood is coded only for Gaussian covariance")
		penalty.set<- c("SCAD")
		if (!is.element(penalty$fun, penalty.set)) stop("At this stage, the penalty #function has to be one of : SCAD")
		if (length(penalty$value)==0) penalty$value <- sqrt(2*log(model@n)/model@n)*seq(from=1, by=0.5, length=15)
		penalty$fun.derivative <- paste(penalty$fun, ".derivative", sep="")
		model@penalty <- penalty
		model@method <- method <- "PMLE"
	}
	
	model@param.estim <- TRUE
	model@optim.method <- as.character(optim.method)
	
	if ((length(lower)==0) | (length(upper)==0)) {
		bounds <- covParametersBounds(model@covariance, design)
		if (length(lower)==0) lower <- bounds$lower
		if (length(upper)==0) upper <- bounds$upper
	}
 	
 	model@lower <- as.numeric(lower)
	model@upper <- as.numeric(upper)
	model@parinit <- as.numeric(parinit)
	
	if (optim.method=="BFGS") {
		if (length(control$pop.size)==0) control$pop.size <- 20
		if (identical(control$trace, FALSE)) {
			control$trace <- 0}
		else control$trace <- 3
	}
	if (optim.method=="gen") {
		d <- ncol(design)
		if (length(control$pop.size)==0) control$pop.size <- min(20, floor(4+3*log(d)))
		if (length(control$max.generations)==0) control$max.generations <- 5
		if (length(control$wait.generations)==0) control$wait.generations <- 2
		if (length(control$BFGSburnin)==0) control$BFGSburnin <- 0
		if (identical(control$trace, FALSE)) {
			control$trace <- 0}
		else control$trace <- 2
	}
	
	upper.alpha <- control$upper.alpha
	if (length(upper.alpha)==0) {
		control$upper.alpha <- 1 - 1e-8
	} else if ((upper.alpha<0) | (upper.alpha>1)) {
		control$upper.alpha <- 1 - 1e-8
	}
			 
	model@control <- control

	model@gr <- as.logical(gr)
 	
 	envir.logLik <- new.env()
	environment(logLikFun) <- environment(logLikGrad) <- envir.logLik
	environment(kmNuggets) <- environment(kmNuggets.init) <- envir.logLik
	environment(kmNoNugget) <- environment(kmNoNugget.init) <- envir.logLik
	environment(km1Nugget) <- environment(km1Nugget.init) <- envir.logLik
	
	validObject(model, complete=TRUE)
	
 	if ((length(noise.var)!=0) | ((length(nugget)!=0) & (nugget.estim==FALSE))) {
		model@case <- "Nuggets"
		f <- kmNuggets
 	}
 	
 	if ((length(nugget)==0) & (nugget.estim == FALSE) & (length(noise.var)==0)) {
		model@case <- "NoNugget"
		f <- kmNoNugget
	} 
	
	if ((length(noise.var)==0) & (nugget.estim == TRUE)) {
		model@case <- "1Nugget"
		f <- km1Nugget
	}  

	if (identical(model@method, "PMLE")) {
		cv <- function(lambda, object, f) {
			object@penalty$value <- lambda
			object@control$trace <- 0
			object <- f(object, envir=envir.logLik)
			criterion <- sum((object@y-leaveOneOut.km(object, type="UK")$mean)^2)
			return(criterion)
		}
		lambda.val <- model@penalty$value
		nval <- length(lambda.val)
		u <- rep(0, nval)
		for (i in 1:nval) {
			u[i] <- cv(lambda.val[i], object=model, f)
		}
		plot(lambda.val, u)
		lambda <- lambda.val[which.min(u)]
		model@penalty$value <- lambda
		model <- f(model, envir=envir.logLik)
	} else {
		model <- f(model, envir=envir.logLik)
	}

	return(model)
}

