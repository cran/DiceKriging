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
		}
		
		if (!identical(object@shape.val, numeric(0))) {
			if (min(object@shape.val)<0) stop("The shape parameters must have positive values")
			if (identical(object@name, "powexp") & (max(object@shape.val)>2)) stop("The exponents must be <= 2 for a Power-Exponential covariance")
		}
 	}
)


## ******************
##   METHODS - coef 
## ******************

if(!isGeneric("coef")) {
  setGeneric(name    = "coef",
             def     = function(object, ...) standardGeneric("coef")
           ##  ,package = "DiceKriging"
             )
}

setMethod("coef", "covTensorProduct", 
          function(object){
          	if (object@paramset.n==1) {
            		res <- c(object@range.val)
            		names(res) <- c(object@range.names)
            	} else {
            		res <- c(object@range.val, object@shape.val)
            		names(res) <- c(object@range.names, object@shape.names)
            	}
            	res
          }
)




## ******************
##   METHODS - show 
## ******************

if(!isGeneric("show")) {
  setGeneric(name    = "show",
             def     = function(object) standardGeneric("show")
             )
}

setMethod("show", "covTensorProduct", 
          function(object){
          	
          	range.names <- paste(object@range.names, "(", object@var.names, ")", sep = "")
            	range.names <- formatC(range.names, width = 12)
            	val.mat <- matrix(object@range.val, object@d, 1)
            	tab <- t(formatC(val.mat, width = 10, digits = 4, format = "f", flag = " "))
            	if (identical(object@known.covparam, "All")) {
            		dimnames(tab) <- list("", range.names)
            	} else {
            		dimnames(tab) <- list("  Estimate", range.names)
            	}
          	
          	if (object@paramset.n==2) {
          		shape.names <- paste(object@shape.names, "(", object@var.names, ")", sep = "")
            		shape.names <- formatC(shape.names, width=12)
            		tab <- rbind(tab, shape.names, deparse.level = 0)
            		tab <- rbind(tab, formatC(object@shape.val, width = 10, digits = 4, format = "f", flag = " "))
            		if (identical(object@known.covparam, "All")) {
            			row.names(tab)[3] <- "          "
            		} else {
            			row.names(tab)[3] <- "  Estimate"
            		}
            	}
	      	    
	      	  	cat("\n")
            	cat("Covar. type  :", object@name, "\n")
            
            	cat("Covar. coeff.:\n")
            	print(t(tab), quote=FALSE)
            	cat("\n")
            
            	if (identical(object@known.covparam, "All")) {
            		cat("Variance:", object@sd2)
            	} else {	         
            		cat("Variance estimate:", object@sd2)
            	}
	         	cat("\n")
	         
	         	if (object@nugget.flag) {
	         		if (object@nugget.estim) {
	         			cat("\nNugget effect estimate:", object@nugget)
	         		} else cat("\nNugget effect :", object@nugget)
	         		cat("\n\n")  
	    		}
	    }
)