## *******************************
##   METHODS - covParametersBounds 
## *******************************

if(!isGeneric("covParametersBounds")) {
  setGeneric(name    = "covParametersBounds",
             def     = function(object, X) standardGeneric("covParametersBounds")
           ##  ,package = "DiceKriging"
             )
}

setMethod("covParametersBounds", "covTensorProduct", 
	function(object, X){
		if (object@paramset.n==1) {
			k <- object@range.n
			lower <- rep(1e-10, k)
			upper <- 2 * diff(apply(X, 2, range))
			upper <- as.vector(upper)
		} else if (identical(object@name, "powexp")) {              
		# coef. order : theta_1, ..., theta_d, p_1, ..., p_d
			lower <- rep(1e-10, object@param.n)
			upper <- 2 * diff(apply(X, 2, range))
			k <- object@shape.n
			upper <- as.vector(c(upper, rep(2, k)))
		} else stop("No default values for covariance parameters bounds, the inputs 'lower' and 'upper' are required")
		return(list(lower=lower, upper=upper))
	}
)

setMethod("covParametersBounds", "covIso", 
	function(object, X){
		if (object@paramset.n==1) {
			lower <- 1e-10
			upper <- 2 * diff(apply(X, 2, range))
			upper <- max(upper)
		} else stop("No default values for covariance parameters bounds, the inputs 'lower' and 'upper' are required")
		return(list(lower=lower, upper=upper))
	}
)

setMethod("covParametersBounds", "covAffineScaling", 
	function(object, X){
    object <- as(extract.covIso(object), "covTensorProduct")
		bounds <- covParametersBounds(object=object, X=X)
    bounds <- list(lower=rep(1/bounds$upper, each=2), upper=rep(1/bounds$lower, each=2))
    return(bounds)
	}
)

setMethod("covParametersBounds", "covScaling", 
  function(object, X){
    knots.n <- sapply(object@knots, length)
    object.tp <- as(extract.covIso(object), "covTensorProduct")
		bounds <- covParametersBounds(object=object.tp, X=X)
    bounds <- list(lower=rep(1/bounds$upper, times=knots.n), upper=rep(1/bounds$lower, times=knots.n))
    return(bounds)
	}
)

