covparam2vect <- function(covStruct) {
	
	if (covStruct@paramset.n==1) {
		param <- covStruct@range.val
	} else param <- c(covStruct@range.val, covStruct@shape.val)
		
	return(param)
}