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


