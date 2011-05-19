computeAuxVariables <- function(model) {
	  
	aux <- covMatrix(model@covariance, X=model@X, noise.var=model@noise.var)
	C <- aux[[1]]
	#if (model@covariance@nugget.flag) {
	#	model@C0 <- aux[[2]] 
	#}
		
   T <- chol(C)
	
   x <- backsolve(t(T), model@y, upper.tri = FALSE)
   M <- backsolve(t(T), model@F, upper.tri = FALSE)
	z <- backsolve(t(T), model@y-model@F%*%as.matrix(model@trend.coef), upper.tri=FALSE) 	
	#Q <- qr.Q(qr(M))
   	#H <- Q %*% t(Q)
   	#z <- x - H %*% x
	
	#model@C <- C
	model@T <- T
	model@z <- as.numeric(z)
	model@M <- M
	return(model)
} 
