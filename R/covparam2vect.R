covparam2vect <- function(covStruct) {
	
	if (!is(covStruct, "covAffineScaling")) {
		if (covStruct@paramset.n==1) {
			param <- covStruct@range.val
		} else param <- c(covStruct@range.val, covStruct@shape.val)
	} else param <- as.numeric(matrix(t(covStruct@eta), 2*covStruct@d, 1))
  #param <- as.numeric(matrix(covStruct@eta, 2*covStruct@d, 1))
	return(param)
}