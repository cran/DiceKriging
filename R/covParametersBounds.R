`covParametersBounds` <-
function(X, covStruct) {
	
	if (covStruct@paramset.n==1) {
		k <- covStruct@range.n
		lower <- rep(1e-10, k)
		upper <- 2 * diff(apply(X, 2, range))
		upper <- as.vector(upper[1:k])
		
	} else if (identical(covStruct@name, "powexp")) {              
		# coef. order : theta_1, ..., theta_d, p_1, ..., p_d
		
		lower <- rep(1e-10, covStruct@param.n)
		upper <- 2 * diff(apply(X, 2, range))
		k <- covStruct@shape.n
		upper <- as.vector(c(upper, rep(2, k)))
			
	} else stop("No default values for covariance parameters bounds, the inputs 'lower' and 'upper' are required")
	
	return(list(lower=lower, upper=upper))
}

