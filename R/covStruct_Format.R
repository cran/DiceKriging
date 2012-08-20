## ****************************************
##   METHODS - covparam2vect, vect2covparam 
## ****************************************

if(!isGeneric("covparam2vect")) {
  setGeneric(name    = "covparam2vect",
             def     = function(object) standardGeneric("covparam2vect")
           ##  ,package = "DiceKriging"
             )
}

covparam2vect.fun <- function(object){
	if (object@paramset.n==1) {
		param <- object@range.val
	} else param <- c(object@range.val, object@shape.val)
	return(as.numeric(param))
}


setMethod("covparam2vect", "covTensorProduct", 
	covparam2vect.fun
)

setMethod("covparam2vect", "covIso", 
	covparam2vect.fun
)


setMethod("covparam2vect", "covAffineScaling", 
	function(object){
		param <- matrix(t(object@eta), 2*object@d, 1)
		return(as.numeric(param))
	}
)

setMethod("covparam2vect", "covScaling", 
	function(object){
		param <- unlist(object@eta)
		return(as.numeric(param))
	}
)



####

if(!isGeneric("vect2covparam")) {
	setGeneric(name    = "vect2covparam",
			   def     = function(object, param) standardGeneric("vect2covparam")
##  ,package = "DiceKriging"
			   )
}

vect2covparam.fun <- function(object, param){
	if (length(param)>0) {
		if (object@paramset.n==1) {
			object@range.val <- param
		} else {	
			range.n <- object@range.n
			object@range.val <- param[1:range.n]
			object@shape.val <- param[(range.n+1):length(param)]
		}
	}
	return(object)
}


setMethod("vect2covparam", "covTensorProduct", 
	vect2covparam.fun
)

setMethod("vect2covparam", "covIso", 
	vect2covparam.fun
)

setMethod("vect2covparam", "covAffineScaling", 
	function(object, param){
		if (length(param)>0) {
			object@eta <- t(matrix(param, 2, object@d))
		}
		return(object)
	}
)

setMethod("vect2covparam", "covScaling", 
function(object, param){
	if (length(param)>0) {
		knots.n <- sapply(object@knots, length)
		ind <- rep(names(knots.n), times=knots.n)
		df <- data.frame(values=param, ind=ind)
		object@eta <- unstack(df)
	}
	return(object)
}
)
