`leaveOneOut.km` <-
function(model, type) {
	
	X <- model@X
	y <- model@y
	T <- model@T
	C <- t(T)%*%T            # should be improved in future versions  
	n <- nrow(as.matrix(X))
	yhat <- sigma2 <- matrix(NA,n,1) 

	beta <- model@trend.coef
	F <- model.matrix(model@trend.formula, data=data.frame(X,y))
	norm2 <- function(u) t(u)%*%u					

	if (type=="UK")	{	# universal kriging
	
		for (i in 1:n) {
			F.i <- F[i,] 
			y.predict.trend <- F.i %*% beta
			
			y.but.i <- y[-i,]
			F.but.i <- F[-i,]
			C.but.i <- C[-i,-i]
			c.but.i <- C[-i,i]
			T.but.i <- chol(C.but.i)					
			x <- backsolve(t(T.but.i), y.but.i, upper.tri=FALSE)
			M <- backsolve(t(T.but.i), F.but.i, upper.tri=FALSE)
			M <- as.matrix(M)
			z <- x - M%*%beta
			Tinv.c <- backsolve(t(T.but.i), c.but.i, upper.tri=FALSE)      # only a vector in this case
			y.predict.complement <- t(Tinv.c) %*% z
			
			y.predict <- y.predict.trend + y.predict.complement
			yhat[i] <- as.numeric(y.predict)	
			
			sigma2.1 <- norm2(Tinv.c)
			#qrR <- qr.R(qr(M))
			T.M <- chol(t(M)%*%M)
			sigma2.mat <- backsolve(t(T.M), t(F.i - t(Tinv.c)%*%M) , upper.tri=FALSE)
			sigma2.2 <- apply(sigma2.mat, 2, norm2)
		
			sigma2[i] <- pmax(model@covariance@sd2 - sigma2.1 + sigma2.2, 0)
			sigma2[i] <- as.numeric(sigma2[i])
			
		}		
	}	
	return(list(mean=yhat, sd=sqrt(sigma2)))
}
