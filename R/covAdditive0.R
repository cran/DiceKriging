## -----------------
## CLASS definitions
## -----------------

## covAdditive0 (Matern 5/2)

# setClass("covAdditive0", 		
#          representation(
#                         d = "integer",            	## (spatial) dimension
#                         name = "character",         ## "matern5_2add0"
#                         var.names = "character",    ## e.g.  c("Lat", "Long") length d
#                         sd2 = "numeric",         	## variance (stationarity)
#                         known.covparam = "character",   ## known covariance parameters (except nugget): "All" or "Known"
#                         param.n = "integer",       ## 2d+1
#                         weight = "numeric",           ## variances for each summand kernel
#                         weight.names = "character",
#                         range.val = "numeric",          ## their values                  
#                         range.names = "character",
#                         nugget = "numeric",        ## nugget (variance)
#                         nugget.flag = "logical",
#                         nugget.estim = "logical"
#                         )
# )



#                        range.n = "integer",            ## number of distinct range parms
#                        range.names = "character",  ## their name (usually "theta")
#                        name = "character",             ## "matern5_2"
#                        paramset.n = "integer",         ## number of parameters sets 
                        ##   gauss, exp : 1;  powexp : 2
			## s.d. of dor the non-nugget part of error
			## nugget part
#                        nugget.flag = "logical",  	## logical : is there a nugget effect ?
#                        nugget.estim = "logical", 	## logical : is it estimated (TRUE) or known ?
  			## total number of parameters (except sigma and nugget)
  			## range part 
  			## shape part, if any 
#                        shape.n = "integer",            ## number of distinct shape parms
#                        shape.names = "character",      ## their name ("p", "nu", "alpha", etc.)
#                        shape.val = "numeric"           ## their values
#                         ),
#          validity = function(object) {
#            
#            covset <- c("gauss", "exp", "matern3_2", "matern5_2", "powexp")
#            if (!is.element(object@name, covset)) {
#              cat("The list of available covariance functions is:\n", covset, "\n")
#              return("invalid character string for 'covtype' argument")
#            }
#            
#            if (!identical(object@sd2, numeric(0))) {
#              if (object@sd2 < 0) {
#                return("The model variance should be non negative")
#              }
#            }
#            
#            if (length(object@nugget) > 1L) {
#              return("'nugget' must be a single non-negative number. For heteroskedastic noisy observations, use 'noise.var' instead.")
#            }
#            
#            if (!identical(object@nugget, numeric(0))) {
#              if (object@nugget < 0) {
#                return("The nugget effect should be non negative")
#              }
#            }
#            
#            if (!identical(object@range.val, numeric(0))) {
#              if (min(object@range.val) < 0) {
#                return("The range parameters must have positive values")
#              }
#              if (length(object@range.val) != object@d) {
#                return("Incorrect number of range parameters")
#              }
#            }
#            
#            if (!identical(object@shape.val, numeric(0))) {
#              if (min(object@shape.val) < 0) {
#                return("The shape parameters must have positive values")
#              }
#              if (length(object@shape.val) != object@d) {
#                return("Incorrect number of shape parameters")
#              }
#              if (identical(object@name, "powexp") && (max(object@shape.val) > 2)) {
#                return("The exponents must be <= 2 for a Power-Exponential covariance")
#              }
#            }
#            TRUE
#          }
#          )


# setMethod("covMatrix", 
#           signature = "covAdditive0", 
#           definition = function(object, X, noise.var=NULL) {
#             n <- nrow(X)
#             d <- object@d
#             C <- matrix(0, nrow=n, ncol=n)
#             for (i in 1:d) {
#               out <- .C("C_covMatrix", 
#                         as.double(X[,i]), as.integer(n), as.integer(1), 
#                         as.double(object@range.val[i]), as.double(object@weight[i]),
#                         as.character("matern5_2"),
#                         ans = double(n * n),
#                         PACKAGE="DiceKriging")
#               Ci <- matrix(out$ans, n, n)          
#               C <- C + Ci
#             }
#             vn <- rep(object@nugget, n)
#             C <- C + diag(vn)  
#             return(list(C=C, vn=vn))
#           } # end definition
# )
# 
# 
# setMethod("covMat1Mat2", 
#           signature = "covAdditive0", 
#           definition = function(object, X1, X2, nugget.flag=FALSE) {
#             n1 <- nrow(X1)
#             n2 <- nrow(X2)
#             d <- object@d
#             M <- matrix(0, nrow=n1, ncol=n2)
#             for (i in 1:d) {
#               out <- .C("C_covMat1Mat2", 
#                         as.double(X1[,i]), as.integer(n1),
#                         as.double(X2[,i]), as.integer(n2), 
#                         as.integer(1),
#                         as.double(object@range.val[i]), as.double(object@weight[i]), 
#                         as.character("matern5_2"),
#                         ans = double(n1 * n2), PACKAGE="DiceKriging")
#               
#               Mi <- matrix(out$ans, n1, n2)          
#               M <- M + Mi
#             }
#             out <- .C("C_covMat1Mat2", 
#                       as.double(X1), as.integer(n1),
#                       as.double(X2), as.integer(n2), 
#                       as.integer(d),
#                       as.double(0), as.double(object@nugget), "whitenoise",
#                       ans = double(n1 * n2), PACKAGE="DiceKriging")
#             N <- matrix(out$ans, n1, n2)
#             return(M+N)
#           } # definition
# )
# 
# 
# setMethod("covparam2vect", 
#           signature ="covAdditive0", 
#           definition = function(object){
#             v <- c(object@range.val, object@weight, object@nugget)
#             names(v) <- c(object@range.names, object@weight.names, "nugget")
#             return(v)
#           }
# )
# 
# setMethod("vect2covparam", 
#           signature = "covAdditive0", 
#           definition = function(object, param){
#             d <- object@d
#             object@range.val <- param[1:d]
#             object@weight <- param[(d+1):(2*d)]
#             object@nugget <- param[2*d+1]
#             object@sd2 <- sum(object@weight)   
#             return(object)
#           }
# )

# setMethod("coef", 
#           signature = signature(object = "covAdditive0"), 
#           definition = function(object, type = "all", as.list = TRUE){
#             val <-
#               switch(type,
#                      "all" = list(range = object@range.val, weight = object@weight, 
#                                   nugget = object@nugget),
#                      NULL
#               )
#             if (as.list) return(val)
#             else return(unlist(val, use.names = FALSE))
#           }
# )
# 
# setMethod("coef<-",
#           signature = signature(object = "covAdditive0"),
#           definition = function(object, type = "all", checkValidity = TRUE,
#                                 ..., value) {
#           }
# )
            

# setMethod("covParametersBounds", 
#           signature = "covAdditive0", 
#           definition = function(object, X){
#             d <- object@d
#             lower <- rep(1e-10, 2*d+1)
#             upper.range <- 2 * diff(apply(X, 2, range))
#             upper <- c(upper.range, rep(1e10, d), 1e10)
#             names(lower) <- names(upper) <- c(object@range.names, object@weight.names, "nugget")            
#             return(list(lower=lower, upper=upper))
#           }
# )
# 
# setMethod("paramSample",
#           signature = "covAdditive0", 
#           definition = function(object, n, lower, upper, y){
#             d <- object@d
#             thetaSample <- matrix(runif(n*d), nrow=d, ncol=n)
#             thetaSample <- lower + thetaSample*(upper[1:d] - lower[1:d])
#             weightSample <- matrix(rchisq(n*(d+1), df=1)*var(y), nrow=d+1, ncol=n)
#             return(rbind(thetaSample, weightSample))
#           }
# )
