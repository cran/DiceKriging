covparam2vect <- function(covStruct) {
	
	if ((!is(covStruct, "covAffineScaling")) & (!is(covStruct, "covScaling"))) {
		if (covStruct@paramset.n==1) {
			param <- covStruct@range.val
		} else param <- c(covStruct@range.val, covStruct@shape.val)
	} else if (is(covStruct, "covAffineScaling")) {
    param <- matrix(t(covStruct@eta), 2*covStruct@d, 1)    
	} else {
    param <- unlist(covStruct@eta)
  }
  
  return(as.numeric(param))
}