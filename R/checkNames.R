checkNames <- function(X1, X2, X1.name="X1", X2.name="X2") {
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  d1 <- ncol(X1)
  d2 <- ncol(X2)
  
  if (d1!=d2) stop(paste(X1.name, "and", X2.name, "must have the same numbers of columns"))
  d <- d1
  
	
  # check names
  nm1 <- unique(colnames(X1))
  if (sum(nchar(nm1))==0) stop(paste(X1.name, "does not contain any names"))
  if (length(nm1)<d) stop(paste("not enough names (ties?) found in", X1.name))
  
  nm2 <- unique(colnames(X2))
  if (sum(nchar(nm2))==0) {
    #warning(paste(X2.name, "is not named. I consider that the columns of", X1.name, "and", X2.name, "correspond to the same variables and in the same order."))
    colnames(X2) <- colnames(X1)
  } else {
    if (length(nm2)<d) stop(paste("not enough names (ties?) found in", X2.name))
    if (!all(nm2 %in% nm1)) stop(paste("one name in", X2.name, "is not in", X1.name))
    if (!all(nm1 %in% nm2)) stop(paste("one name in", X1.name, "is not in", X2.name))
    ind <- 1:length(nm1)
    names(ind) <- nm1
    X2 <- X2[, ind[nm2], drop=FALSE]
  }
  return(X2)
}