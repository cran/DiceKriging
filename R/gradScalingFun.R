dfj_1pt_2param <- function(t, knots){ 
    dfj <- vector(length=length(knots))
    zeta0 <- knots[1]
    zeta1 <- knots[2]
    aux <- (t-zeta0)^2/(2*(zeta1-zeta0))
    dfj[1] <- (t-zeta0) - aux 
    dfj[2] <- aux 
    return(dfj)
}

scalingFunGrad <- function(X, knots, k) {
  df <- apply(X[, k, drop=FALSE], 1, dfj_1pt_2param, knots)
  return(t(df))
}


