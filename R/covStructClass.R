
# -----------------
# CLASS definitions
# -----------------

# covTensorProduct : separable (or tensor product) covariances, depending on 1 set of parameters
#     Examples - Gaussian, Exponential, Matern with fixed nu=p+1/2, Power-Exponential
#

setClass("covTensorProduct", 		
	representation(
		d = "integer",            	## (spatial) dimension
		name = "character",			## "powexp"
		paramset.n = "integer", 		## number of parameters sets 
											##   gauss, exp : 1;  powexp : 2
		var.names = "character",  	## e.g.  c("Lat", "Long") length d
			## s.d. of dor the non-nugget part of error
		sd2 = "numeric",       		## variance (stationarity)
			## nugget part
		known.covparam = "character",  ## known covariance parameters (except nugget): "All" or "Known"
		nugget.flag = "logical",  	## logical : is there a nugget effect ?
  		nugget.estim = "logical", 	## logical : is it estimated (TRUE) or known ?
  		nugget = "numeric",    		## nugget (variance)
  			## total number of parameters (except sigma and nugget)
  		param.n = "integer",			## range.n + shape.n
  			## range part 
  		range.n = "integer",      	## number of distinct range parms
  		range.names = "character",	## their name (usually "theta")
  		range.val = "numeric",    	## their values
  			## shape part, if any 
  		shape.n = "integer",      	## number of distinct shape parms
  		shape.names = "character",	## their name ("p", "nu", "alpha", etc.)
  		shape.val = "numeric" 		## their values
	),
	validity = function(object) {
 	
 		covset <- c("gauss", "exp", "matern3_2", "matern5_2", "powexp")
		if (!is.element(object@name, covset)) {
			cat("The list of available covariance functions is:\n", covset, "\n")
			stop("invalid character string for 'covtype' argument")
		}
	
		if (!identical(object@sd2, numeric(0))) {
			if (object@sd2<0) stop("The model variance should be non negative")
		}
		
		if (length(object@nugget)>1) stop("Nugget must be a single non-negative number. For heteroskedastic noisy observations, use noise.var instead.")
		
		if (!identical(object@nugget, numeric(0))) {
			if (object@nugget<0) stop("The nugget effect should be non negative")
		}
	
		if (!identical(object@range.val, numeric(0))) {
			if (min(object@range.val)<0) stop("The range parameters must have positive values")
			if (length(object@range.val)!=object@d) stop("Incorrect number of range parameters")
		}
		
		if (!identical(object@shape.val, numeric(0))) {
			if (min(object@shape.val)<0) stop("The shape parameters must have positive values")
			if (length(object@shape.val)!=object@d) stop("Incorrect number of shape parameters")
			if (identical(object@name, "powexp") & (max(object@shape.val)>2)) stop("The exponents must be <= 2 for a Power-Exponential covariance")
		}
 	}
)

setClass("covIso", 	
	representation(
		d = "integer",            	## (spatial) dimension
		name = "character",				## "gauss"
		paramset.n = "integer", 		## number of parameters sets 
											##   gauss, exp : 1;  powexp : 2
		var.names = "character",  	## e.g.  c("Lat", "Long") length d
			## s.d. of the non-nugget part of error
		sd2 = "numeric",       		## variance (stationarity)
			## nugget part
		known.covparam = "character",  ## known covariance parameters (except nugget): "All" or "Known"
		nugget.flag = "logical",  	## logical : is there a nugget effect ?
  		nugget.estim = "logical", 	## logical : is it estimated (TRUE) or known ?
  		nugget = "numeric",    		## nugget (variance)
  			## total number of parameters (except sigma and nugget)
  		param.n = "integer",			## 1
  			## range part 
  		range.names = "character",	## their name (usually "theta")
  		range.val = "numeric"     	## their values
	),
	validity = function(object) {
 	
 		covset <- c("gauss", "exp", "matern3_2", "matern5_2")
		if (!is.element(object@name, covset)) {
			cat("The list of available covariance functions is:\n", covset, "\n")
			stop("invalid character string for 'covtype' argument")
		}
	
		if (!identical(object@sd2, numeric(0))) {
			if (object@sd2<0) stop("The model variance should be non negative")
		}
		
		if (length(object@nugget)>1) stop("Nugget must be a single non-negative number. For heteroskedastic noisy observations, use noise.var instead.")
		
		if (!identical(object@nugget, numeric(0))) {
			if (object@nugget<0) stop("The nugget effect should be non negative")
		}
	
		if (!identical(object@range.val, numeric(0))) {
			if (length(object@range.val)>1) stop("Only one positive value for an isotropic kernel")
			if (min(object@range.val)<0) stop("The range parameter must have positive values")
		}
		
 	}
)


setClass("covAffineScaling", 	
	representation(
		d = "integer",				## (spatial) dimension
		knots = "numeric",			## 2 nodes 
		eta = "matrix",				## value at nodes, matrix d x 2
		name = "character",			## "gauss"
		paramset.n = "integer", 	## number of parameters sets 
			##   gauss, exp : 1;  powexp : 2
		var.names = "character",  	## e.g.  c("Lat", "Long") length d
			## s.d. of the non-nugget part of error
		sd2 = "numeric",       		## variance (stationarity)
			## nugget part
		known.covparam = "character",  ## known covariance parameters (except nugget): "All" or "Known"
		nugget.flag = "logical",  	## logical : is there a nugget effect ?
  		nugget.estim = "logical", 	## logical : is it estimated (TRUE) or known ?
  		nugget = "numeric",    		## nugget (variance)
  			## total number of parameters (except sigma and nugget)
  		param.n = "integer"			## 2*d
	),
	validity = function(object) {
 	
 		covset <- c("gauss", "exp", "matern3_2", "matern5_2")
		if (!is.element(object@name, covset)) {
			cat("The list of available covariance functions is:\n", covset, "\n")
			stop("invalid character string for 'covtype' argument")
		}
	
		if (!identical(object@sd2, numeric(0))) {
			if (object@sd2<0) stop("The model variance should be non negative")
		}
		
		if (length(object@nugget)>1) stop("Nugget must be a single non-negative number. For heteroskedastic noisy observations, use noise.var instead.")
		
		if (!identical(object@nugget, numeric(0))) {
			if (object@nugget<0) stop("The nugget effect should be non negative")
		}
		
 	}
)


setClass("covScaling",   
	representation(
		d = "integer",				## (spatial) dimension
		knots = "list",			## a named list containing the knots
		eta = "list",				## a list containing the values at knots
		#knots.n = "integer",    ## number of knots for each dimension
    name = "character",			## "gauss"
		paramset.n = "integer", 	## number of parameters sets 
			##   gauss, exp : 1;  powexp : 2
		var.names = "character",  	## e.g.  c("Lat", "Long") length d
			## s.d. of the non-nugget part of error
		sd2 = "numeric",       		## variance (stationarity)
			## nugget part
		known.covparam = "character",  ## known covariance parameters (except nugget): "All" or "Known"
		nugget.flag = "logical",  	## logical : is there a nugget effect ?
  		nugget.estim = "logical", 	## logical : is it estimated (TRUE) or known ?
  		nugget = "numeric",    		## nugget (variance)
  			## total number of parameters (except sigma and nugget)
  		param.n = "integer"			## length of knots
	),
	validity = function(object) {
 	
 	 covset <- c("gauss", "exp", "matern3_2", "matern5_2")
		if (!is.element(object@name, covset)) {
			cat("The list of available covariance functions is:\n", covset, "\n")
			stop("invalid character string for 'covtype' argument")
		}
    
    names.knots <- names(object@knots)
    n.knots <- length(object@knots)
    
    if (n.knots>0) {
      if (length(names.knots)==0) {
        stop("the list containing the knots must be named")
      } else if (!all(sapply(object@knots, is.numeric))) {
        stop("the knots must be numeric")
      } else if ((!all(names.knots%in%object@var.names)) | (!all(object@var.names%in%names.knots))) {
        stop("mismatch between the names of knots and input variables")
      }
      for (i in 1:n.knots) {
        knots.i <- object@knots[[i]]
        if (any(is.na(knots.i))) stop("missing values not allowed in knots")
        if (is.unsorted(knots.i)) stop("the knots must be sorted")
        if (any(diff(knots.i) <= 0.0)) stop("dupplicated values in knots")
        if (length(knots.i) < 2L) stop("knots must be of length >=2") 
      }
      
      n.eta <- length(object@eta)
      if (n.eta>0) {
        if (n.eta != n.knots) {
          stop("the number of values at knots is different from the number of knots")
        } else if (!all(sapply(object@eta, is.numeric))) {
          stop("the values at knots must be numeric")
        } 
        for (i in 1:n.eta) {
          eta.i <- object@eta[[i]]
          if (any(eta.i <= 0.0)) stop("the values at knots must be positive")
        }
      }
    }

		if (!identical(object@sd2, numeric(0))) {
			if (object@sd2<0) stop("The model variance should be non negative")
		}
		
		if (length(object@nugget)>1) stop("Nugget must be a single non-negative number. For heteroskedastic noisy observations, use noise.var instead.")
		
		if (!identical(object@nugget, numeric(0))) {
			if (object@nugget<0) stop("The nugget effect should be non negative")
		}
		
 	}
)


setClassUnion("covKernel", c("covTensorProduct", "covIso", "covAffineScaling", "covScaling"))  

