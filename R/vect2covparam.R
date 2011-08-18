vect2covparam <- function(param, covStruct) {
	
		if (length(param)>0) {
	
		if ((!is(covStruct, "covAffineScaling")) & (!is(covStruct, "covScaling"))) {
			
			if (covStruct@paramset.n==1) {
				covStruct@range.val <- param
			} else {	
				range.n <- covStruct@range.n
				covStruct@range.val <- param[1:range.n]
				covStruct@shape.val <- param[(range.n+1):length(param)]
			}
	
		} else if (is(covStruct, "covAffineScaling")) {
			covStruct@eta <- t(matrix(param, 2, covStruct@d))
			#covStruct@eta <- matrix(param, covStruct@d, 2)
		} else {
      knots.n <- sapply(covStruct@knots, length)
      ind <- rep(names(knots.n), times=knots.n)
      df <- data.frame(values=param, ind=ind)
      covStruct@eta <- unstack(df)
		}
	}
	
	return(covStruct)
  
}