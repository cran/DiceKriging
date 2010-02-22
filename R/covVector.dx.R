covVector.dx <- function(x, X, covStruct, c) {

	n <- nrow(X);
	d <- ncol(X);

	param <- covparam2vect(covStruct)

	out <- .C("C_covVector_dx", 
				as.double(x), as.double(X), 
				as.integer(n), as.integer(d),
				as.double(param), as.character(covStruct@name), 
				as.double(c), 
				ans = double(n * d))
   	
   	return(matrix(out$ans, n, d))
	
}
