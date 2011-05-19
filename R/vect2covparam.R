vect2covparam <- function(param, covStruct) {
	
	param <- as.numeric(param)

	if (!identical(param, numeric(0))) {
	
		if (!is(covStruct, "covAffineScaling")) {
			
			if (covStruct@paramset.n==1) {
				covStruct@range.val <- param
			} else {	
				range.n <- covStruct@range.n
				covStruct@range.val <- param[1:range.n]
				covStruct@shape.val <- param[(range.n+1):length(param)]
			}
	
		} else {
			covStruct@eta <- t(matrix(param, 2, covStruct@d))
			#covStruct@eta <- matrix(param, covStruct@d, 2)
		}
	}
		
	return(covStruct)

}