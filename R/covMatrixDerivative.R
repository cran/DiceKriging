`covMatrixDerivative` <-
function(X, covStruct, C, k) {
	# X : n x d
	n <- nrow(X)
	d <- ncol(X)
				
	# beware that the index k  starts at 0 in C language 
	
	param <- covparam2vect(covStruct)
	
	out <- .C("C_covMatrixDerivative", 
					as.double(X), as.integer(n), as.integer(d), 
					as.double(param), as.character(covStruct@name),
					as.integer(k), as.double(C),
					ans = double(n * n),
					PACKAGE="DiceKriging")

	return(matrix(out$ans, n, n))
	
}

