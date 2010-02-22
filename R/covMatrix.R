`covMatrix` <-
function(X, covStruct, noise.var=NULL) {
		
	d <- ncol(X)
	n <- nrow(X)
	
	param <- covparam2vect(covStruct)
		
	out <- .C("C_covMatrix", 
					as.double(X), as.integer(n), as.integer(d), 
					as.double(param), as.double(covStruct@sd2), as.character(covStruct@name), 
					ans = double(n * n),
					PACKAGE="DiceKriging")
   	
	C <- matrix(out$ans, n, n)
	
	C0 <- C			# covariance matrix when there is no nugget effect
	
	if (covStruct@nugget.flag) {
		C <- C0 + diag(covStruct@nugget, n, n)
	} else if (length(noise.var)>0) {
		C <- C0 + diag(noise.var)
	} else C <- C0
		
	return(list(C,C0))	
}

