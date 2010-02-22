vect2covparam <- function(param, covStruct) {
	
	param <- as.numeric(param)

	if (!identical(param, numeric(0))) {
	
		if (covStruct@paramset.n==1) {
			covStruct@range.val <- param
		} else {	
			range.n <- covStruct@range.n
			covStruct@range.val <- param[1:range.n]
			covStruct@shape.val <- param[(range.n+1):length(param)]
		}
	
	}
		
	return(covStruct)

}