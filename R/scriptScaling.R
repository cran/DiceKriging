sf1d <- function(x, knots, eta){
# Maybe split it into 2 or more functions 
# on the longer term

############################################
################# TESTS #################### 
# equality test between the elements of knots
if(prod(diff(knots))==0){
stop("There are (at least) two overlapping knots. 
Knots must not overlap!")}

# ordering test for the elements of knots
if(is.unsorted(knots)){
ordered_knots <- sort(knots,index.return = TRUE)
knots <- ordered_knots$x
eta <- eta[ordered_knots$ix]}
# stop("the vector of knots knots must be ordered!")}

# testing the fact that x lies within knots's range
if( prod(x >= min(knots))*prod(x <= max(knots))==0 ){
stop("x must lie within knots's range!")}
# Test has to be made for all potential inputs!
############################################

# knots is a set of knots, say K + 1 of them. 
#"Length" to be replaced in higher dims. 
K <- length(knots)-1

# eta is a set of transformation parameters. 
# eta should also be K+1 dimensional
Restricted_eta <- eta[seq(2,K+1)]
Shifted_eta <- eta[seq(1,K)]

# Same shift/restriction for knots:
Restricted_knots <- knots[seq(2,K+1)]
Shifted_knots <- knots[seq(1,K)]

# Computation of the vectors of coefficients 
# (defining the piecewise linear density)
a <- (Restricted_knots*Shifted_eta
-Shifted_knots*Restricted_eta)/
(Restricted_knots-Shifted_knots)

b <- (Restricted_eta-Shifted_eta)/
(Restricted_knots-Shifted_knots)

# Trick to compute the integral

tricky <- matrix(0,nrow=K, ncol=length(x))

for(i in seq(1,K)){
tricky[i,] <- pmax(knots[i], pmin(x,knots[i+1]))
}

transformed_x <-  as.vector(a%*%(tricky-Shifted_knots)+ 
0.5*b%*%(tricky^2-Shifted_knots^2))

return(transformed_x)
}





scalingFun <- function(X, knots, eta) {
	# Here X is meant to be a n*d matrix,
# knots is a (K+1) vector (in this special case),
# and eta is a d*(K+1) matrix of coefficients.
d <- dim(X)[2]

transformed_X <- matrix(NA, nrow=nrow(X), ncol=ncol(X))

for(i in seq(1,d)){
transformed_X[,i] <- sf1d(X[,i], knots, eta[i,])
}

return(transformed_X)

}