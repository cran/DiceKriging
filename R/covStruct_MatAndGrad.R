## ***********************
##   METHODS - covMatrix
## ***********************

covMatrix.covTensorProduct <- function(object, X, noise.var=NULL) {
	
	d <- ncol(X)
	n <- nrow(X)
	
	param <- covparam2vect(object)
		
	out <- .C("C_covMatrix", 
			  as.double(X), as.integer(n), as.integer(d), 
			  as.double(param), as.double(object@sd2), as.character(object@name), 
			  ans = double(n * n),
			  PACKAGE="DiceKriging")
   	
	C <- matrix(out$ans, n, n)
	
	C0 <- C			# covariance matrix when there is no nugget effect
	
	if (object@nugget.flag) {
		C <- C0 + diag(object@nugget, n, n)
	} else if (length(noise.var)>0) {
		C <- C0 + diag(noise.var)
	} else C <- C0
	
	return(list(C,C0))	
}


if(!isGeneric("covMatrix")) {
	setGeneric(name    = "covMatrix",
			   def     = function(object, X, noise.var=NULL) standardGeneric("covMatrix")
			   )
}

setMethod("covMatrix", "covTensorProduct", 
function(object, X, noise.var=NULL) {
	covMatrix.covTensorProduct(object=object, X=X, noise.var=noise.var)
}
)

setMethod("covMatrix", "covIso", 
function(object, X, noise.var=NULL) {
	covMatrix.covTensorProduct(object=as(object, "covTensorProduct"), X=X, noise.var=noise.var)
}
)

setMethod("covMatrix", "covAffineScaling", 
function(object, X, noise.var=NULL) {
	covMatrix(extract.covIso(object), X=affineScalingFun(X, knots=object@knots, eta=object@eta), noise.var=noise.var)
}
)

setMethod("covMatrix", "covScaling", 
function(object, X, noise.var=NULL) {
  covMatrix(extract.covIso(object), X=scalingFun(X, knots=object@knots, eta=object@eta), noise.var=noise.var)
}
)


## *******************************
##   METHODS - covMatrixDerivative
## *******************************

covMatrixDerivative.covTensorProduct <- function(object, X, C, k) {
	# X : n x d
	n <- nrow(X)
	d <- ncol(X)
				
	# beware that the index k  starts at 0 in C language 
	
	param <- covparam2vect(object)
	
	out <- .C("C_covMatrixDerivative", 
					as.double(X), as.integer(n), as.integer(d), 
					as.double(param), as.character(object@name),
					as.integer(k), as.double(C),
					ans = double(n * n),
					PACKAGE="DiceKriging")

	return(matrix(out$ans, n, n))
}

covMatrixDerivative.dx.covTensorProduct <- function(object, X, C, k) {
  # X : n x d
	n <- nrow(X)
	d <- ncol(X)
				
	# beware that the index k  starts at 0 in C language 
	
	param <- covparam2vect(object)
	
	out <- .C("C_covMatrixDerivative_dx", 
					as.double(X), as.integer(n), as.integer(d), 
					as.double(param), as.character(object@name),
					as.integer(k), as.double(C),
					ans = double(n * n),
					PACKAGE="DiceKriging")

	return(matrix(out$ans, n, n))
}


if(!isGeneric("covMatrixDerivative")) {
	setGeneric(name = "covMatrixDerivative",
			   def = function(object, X, C, k, ...) standardGeneric("covMatrixDerivative")
			   )
}

setMethod("covMatrixDerivative", "covTensorProduct", 
function(object, X, C, k) {
	covMatrixDerivative.covTensorProduct(object=object, X=X, C=C, k=k)
}
)

setMethod("covMatrixDerivative", "covIso", 
function(object, X, C, k) {
	dC <- matrix(0, nrow(C), ncol(C))
	object <- as(object, "covTensorProduct")
	for (j in 1:object@d) {
		dC <- dC + covMatrixDerivative.covTensorProduct(object=object, X=X, C=C, k=j)
	}
	return(dC)
}
)

envir.covScaling <- envir.covAffineScaling <- new.env()

setMethod("covMatrixDerivative", "covAffineScaling", 
function(object, X, C, k, envir=envir.covAffineScaling) {
  # NOTE : this function MUST be used in a loop over the index k, from 1 to k.max
  i <- k
  k <- floor((i-1)/2)+1
  l <- (i-1)%%2+1
  if (l==1) {
   object.covTensorProduct <- as(extract.covIso(object), "covTensorProduct")
   fX <- affineScalingFun(X, knots=object@knots, eta=object@eta)
   Dk <- covMatrixDerivative.dx.covTensorProduct(object=object.covTensorProduct, 
                       X=fX, C=C, k=k)
   assign("Dk", Dk, envir=envir)
  } else {
    Dk <- get("Dk", envir=envir) 
  }
  df.dkl <- affineScalingGrad(X=X, knots=object@knots, k=k)[,l]
  A <- outer(df.dkl, df.dkl, "-")
  return(Dk*A)
}
)

setMethod("covMatrixDerivative", "covScaling", 
function(object, X, C, k, envir=envir.covScaling) {
  # NOTE : this function MUST be used in a loop over the index k, from 1 to k.max
  i <- k
  if (i==1) {
    knots.n <- as.numeric(sapply(object@knots, length))
    k.vec <- rep(1:object@d, times=knots.n)
    l.vec <- sequence(knots.n)
    assign("k.vec", k.vec, envir=envir)
    assign("l.vec", l.vec, envir=envir)
  } else {
    k.vec <- get("k.vec", envir=envir)
    l.vec <- get("l.vec", envir=envir)
  }
  k <- k.vec[i]
  l <- l.vec[i]
  
  if (l==1) {
    object.covTensorProduct <- as(extract.covIso(object), "covTensorProduct")
    fX <- scalingFun(X, knots=object@knots, eta=object@eta)
    Dk <- covMatrixDerivative.dx.covTensorProduct(object=object.covTensorProduct, 
                       X=fX, C=C, k=k)   
    assign("Dk", Dk, envir=envir)
  } else {
    Dk <- get("Dk", envir=envir) 
  }
  df.dkl <- scalingGrad(X=X, knots=object@knots, k)[,l]
  A <- outer(df.dkl, df.dkl, "-")
  return(Dk*A)
}
)




## ************************
##   METHODS - covVector.dx
## ************************

covVector.dx.covTensorProduct <- function(object, x, X, c) {

	n <- nrow(X);
	d <- ncol(X);

	param <- covparam2vect(object)

	out <- .C("C_covVector_dx", 
				as.double(x), as.double(X), 
				as.integer(n), as.integer(d),
				as.double(param), as.character(object@name), 
				as.double(c), 
				ans = double(n * d))
   	
   	return(matrix(out$ans, n, d))
	
}

if(!isGeneric("covVector.dx")) {
	setGeneric(name = "covVector.dx",
			   def = function(object, x, X, c) standardGeneric("covVector.dx")
			   )
}

setMethod("covVector.dx", "covTensorProduct", 
function(object, x, X, c) {
	covVector.dx.covTensorProduct(object=object, x=x, X=X, c=c)
}
)

setMethod("covVector.dx", "covIso", 
function(object, x, X, c) {
	covVector.dx.covTensorProduct(object=as(object, "covTensorProduct"), x=x, X=X, c=c)
}
)


## ***********************
##   METHODS - covMat1Mat2
## ***********************

covMat1Mat2.covTensorProduct <- function(object, X1, X2, nugget.flag=FALSE) {

	# X1 : matrix n1 x d - containing training points
	# X2 : matrix n2 x d - containing test points
	
	# X1 <- as.matrix(X1)
	# X2 <- checkNames(X1, X2)
	
	n1 <- nrow(X1)
	n2 <- nrow(X2)
	d <- ncol(X1)
	
	param <- covparam2vect(object)
	
	out <- .C("C_covMat1Mat2", 
				as.double(X1), as.integer(n1),
				as.double(X2), as.integer(n2), 
				as.integer(d),
				as.double(param), as.double(object@sd2), as.character(object@name),
				ans = double(n1 * n2), PACKAGE="DiceKriging")
				
	M <- matrix(out$ans, n1, n2)
	
	if ((!nugget.flag) | (!object@nugget.flag)) {
		return(M)
	} else {
		out <- .C("C_covMat1Mat2", 
				as.double(X1), as.integer(n1),
				as.double(X2), as.integer(n2), 
				as.integer(d),
				as.double(param), as.double(object@nugget), "whitenoise",
				ans = double(n1 * n2), PACKAGE="DiceKriging")
		N <- matrix(out$ans, n1, n2)
		return(M+N)
	}
		
}

if(!isGeneric("covMat1Mat2")) {
	setGeneric(name    = "covMat1Mat2",
			   def     = function(object, X1, X2, nugget.flag=FALSE) standardGeneric("covMat1Mat2")
			   )
}

setMethod("covMat1Mat2", "covTensorProduct", 
function(object, X1, X2, nugget.flag=FALSE) {
	covMat1Mat2.covTensorProduct(object=object, X1=X1, X2=X2, nugget.flag=nugget.flag)
}
)

setMethod("covMat1Mat2", "covIso", 
function(object, X1, X2, nugget.flag=FALSE) {
	covMat1Mat2.covTensorProduct(object=as(object, "covTensorProduct"), X1=X1, X2=X2, nugget.flag=nugget.flag)
}
)

setMethod("covMat1Mat2", "covAffineScaling", 
function(object, X1, X2, nugget.flag=FALSE) {
	covMat1Mat2(extract.covIso(object), X1=affineScalingFun(X1, knots=object@knots, eta=object@eta), X2=affineScalingFun(X2, knots=object@knots, eta=object@eta), nugget.flag=nugget.flag)
}
)

setMethod("covMat1Mat2", "covScaling", 
function(object, X1, X2, nugget.flag=FALSE) {
  covMat1Mat2(extract.covIso(object), X1=scalingFun(X1, knots=object@knots, eta=object@eta), X2=scalingFun(X2, knots=object@knots, eta=object@eta), nugget.flag=nugget.flag)
}
)
