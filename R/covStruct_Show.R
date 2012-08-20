## ******************
##   METHODS - show 
## ******************

if(!isGeneric("show")) {
  setGeneric(name    = "show",
             def     = function(object) standardGeneric("show")
             )
}

setMethod("show", "covTensorProduct", 
          function(object){
          	
          	range.names <- paste(object@range.names, "(", object@var.names, ")", sep = "")
            	range.names <- formatC(range.names, width = 12)
            	val.mat <- matrix(object@range.val, object@d, 1)
            	tab <- t(formatC(val.mat, width = 10, digits = 4, format = "f", flag = " "))
            	if (identical(object@known.covparam, "All")) {
            		dimnames(tab) <- list("", range.names)
            	} else {
            		dimnames(tab) <- list("  Estimate", range.names)
            	}
          	
          	if (object@paramset.n==2) {
          		shape.names <- paste(object@shape.names, "(", object@var.names, ")", sep = "")
            		shape.names <- formatC(shape.names, width=12)
            		tab <- rbind(tab, shape.names, deparse.level = 0)
            		tab <- rbind(tab, formatC(object@shape.val, width = 10, digits = 4, format = "f", flag = " "))
            		if (identical(object@known.covparam, "All")) {
            			row.names(tab)[3] <- "          "
            		} else {
            			row.names(tab)[3] <- "  Estimate"
            		}
            	}
	      	    
	      	  	cat("\n")
            	cat("Covar. type  :", object@name, "\n")
            
            	cat("Covar. coeff.:\n")
            	print(t(tab), quote=FALSE)
            	cat("\n")
            
            	if (identical(object@known.covparam, "All")) {
            		cat("Variance:", object@sd2)
            	} else {	         
            		cat("Variance estimate:", object@sd2)
            	}
	         	cat("\n")
	         
	         	if (object@nugget.flag) {
	         		if (object@nugget.estim) {
	         			cat("\nNugget effect estimate:", object@nugget)
	         		} else cat("\nNugget effect :", object@nugget)
	         		cat("\n\n")  
	    		}
	    }
)



setMethod("show", "covIso", 
          function(object){
          	
          		range.names <- object@range.names
            	range.names <- formatC(range.names, width = 12)
            	val.mat <- matrix(object@range.val, 1, 1)
            	tab <- t(formatC(val.mat, width = 10, digits = 4, format = "f", flag = " "))
            	if (identical(object@known.covparam, "All")) {
            		dimnames(tab) <- list("", range.names)
            	} else {
            		dimnames(tab) <- list("  Estimate", range.names)
            	}
	      	    
	      	  	cat("\n")
            	cat("Covar. type  :", object@name)
           		cat(", isotropic")
            	cat("\n")
            
            	cat("Covar. coeff.:\n")
            	print(t(tab), quote=FALSE)
            	cat("\n")
            
            	if (identical(object@known.covparam, "All")) {
            		cat("Variance:", object@sd2)
            	} else {	         
            		cat("Variance estimate:", object@sd2)
            	}
	         	cat("\n")
	         
	         	if (object@nugget.flag) {
	         		if (object@nugget.estim) {
	         			cat("\nNugget effect estimate:", object@nugget)
	         		} else cat("\nNugget effect :", object@nugget)
	         		cat("\n\n")  
	    		}
	    }
)


setMethod("show", "covAffineScaling", 
          function(object){
          	
          		range.names <- paste("eta", "(", object@var.names, ")", sep = "")
            	range.names <- formatC(range.names, width = 12)
            	tab <- t(formatC(object@eta, width = 10, digits = 4, format = "f", flag = " "))
            	if (identical(object@known.covparam, "All")) {
            		dimnames(tab) <- list(c("",""), range.names)
            	} else {
            		dimnames(tab) <- list(c("  Estimate", "  Estimate"), range.names)
            	}
          	    
	      	  	cat("\n")
            	cat("Covar. type  :", object@name, ", with affine scaling \n")
            
            	cat("Covar. coeff.:\n")
            	print(t(tab), quote=FALSE)
            	cat("\n")
            
            	if (identical(object@known.covparam, "All")) {
            		cat("Variance:", object@sd2)
            	} else {	         
            		cat("Variance estimate:", object@sd2)
            	}
	         	cat("\n")
	         
	         	if (object@nugget.flag) {
	         		if (object@nugget.estim) {
	         			cat("\nNugget effect estimate:", object@nugget)
	         		} else cat("\nNugget effect :", object@nugget)
	         		cat("\n\n")  
	    		}
	    }
)

setMethod("show", "covScaling", 
          function(object){
 
              cat("\n")
              cat("Covar. type  :", object@name, ", with scaling \n")
            
              cat("Covar. coeff.")
              if (!identical(object@known.covparam, "All")) cat(", with estimated values for eta")
              cat(":\n")
 
          		for (i in 1:object@d) {
                knots.names <- paste("knots", "(", object@var.names[i], ")", sep = "")
                eta.names <- paste("eta", "(", object@var.names[i], ")", sep = "")
                param.names <- c(eta.names, knots.names)
                param.names <- formatC(param.names, width = 12)
                tab <- t(formatC(cbind(object@eta[[i]], object@knots[[i]]), width = 10, digits = 4, format = "f", flag = " "))
                n.i <- length(object@knots[[i]])
              	dimnames(tab) <- list(param.names, rep("", n.i))
            	 print(tab, quote=FALSE)
           	}
               
          	    
	      	  	cat("\n")
            
            	if (identical(object@known.covparam, "All")) {
            		cat("Variance:", object@sd2)
            	} else {	         
            		cat("Variance estimate:", object@sd2)
            	}
	         	cat("\n")
	         
	         	if (object@nugget.flag) {
	         		if (object@nugget.estim) {
	         			cat("\nNugget effect estimate:", object@nugget)
	         		} else cat("\nNugget effect :", object@nugget)
	         		cat("\n\n")  
	    		}
	    }
)

setMethod("show", "covUser", 
          function(object){
            cat("\n")
            cat("Covar. type  : user type \n")                    
            print(object@kernel)
            
            if (object@nugget.flag) {
              cat("\nNugget effect :", object@nugget)
              cat("\n\n")
            }
          }
)




