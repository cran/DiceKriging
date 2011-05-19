`covStruct.create` <- 
function(covtype, d, known.covparam, coef.cov=NULL, coef.var=NULL, nugget=NULL, nugget.estim=FALSE, nugget.flag=FALSE, iso=FALSE, scaling=FALSE, var.names=NULL) {
	
  if (scaling & iso) {
    iso <- FALSE
    warning("At this stage no isotropic version is available, regular scaling is applied.")
  }
  
	covsetI <- c("gauss", "exp", "matern3_2", "matern5_2")
	covsetII <- c("powexp")

	classeType <- "covTensorProduct"
	if (iso) classeType <- "covIso"
	if (scaling) classeType <- "covAffineScaling"
			 
	covStruct <- new(classeType, d=as.integer(d), name=as.character(covtype), 
		sd2 = as.numeric(coef.var), var.names=as.character(var.names), 
		nugget = as.double(nugget), nugget.flag=nugget.flag, nugget.estim=nugget.estim, known.covparam=known.covparam) 
									
	if (!scaling) {					
		
		covStruct@range.names  = "theta"
		
		if (is.element(covtype, covsetI)) {
			covStruct@paramset.n <- as.integer(1)
			if (iso) {
				covStruct@param.n <- as.integer(1)
			} else {
				covStruct@param.n <- as.integer(d)
				covStruct@range.n <- as.integer(d)
			}	
		} else {	
			covStruct@paramset.n <- as.integer(2)
			covStruct@param.n <- as.integer(2*d)
			covStruct@range.n <- as.integer(d)
			covStruct@shape.n <- as.integer(d)
			covStruct@shape.names <- "p"
		}

		if (length(coef.cov)>0) covStruct <- vect2covparam(coef.cov, covStruct)

	} else {
   	   	
		covStruct@paramset.n <- as.integer(1)
		covStruct@param.n <- as.integer(2*d)
		covStruct@knots <- c(0,1)
		if (!is.null(coef.cov)) covStruct@eta <- coef.cov
    
	}
			 	
	return(covStruct)	
}

