`covStruct.create` <- 
function(covtype, d, coef.cov=NULL, coef.var=NULL, nugget=NULL) {
	
	covsetI <- c("gauss", "exp", "matern3_2", "matern5_2")
	covsetII <- c("powexp")
	if (identical(covtype, "powexp")) {
		shape.names <- "p"
	} else shape.names <- character(0)
	
	if (is.element(covtype, covsetI)) {

		covStruct <- new("covTensorProduct", d = as.integer(d), name=as.character(covtype),
									  paramset.n   = as.integer(1),
		                        param.n      = as.integer(d),
		                        sd2          = as.double(coef.var),
		                        nugget.flag  = !is.null(nugget),
		                        nugget.estim = FALSE,
		                        nugget	      = as.double(nugget), 
		                        range.n      = as.integer(d),
		                        range.names  = "theta") 

	} else {	

		covStruct <- new("covTensorProduct", d = as.integer(d), name=as.character(covtype),
										paramset.n   = as.integer(2),
										param.n      = as.integer(2*d),
										sd2          = as.double(coef.var),
  		                         nugget.flag  = !is.null(nugget),
  		                         nugget.estim = FALSE,
		                         nugget       = as.double(nugget), 
										range.n      = as.integer(d),
										range.names  = "theta",
		                        	shape.n      = as.integer(d),
		                        	shape.names  = shape.names)
   	}
	
	covStruct <- vect2covparam(coef.cov, covStruct)
	
	return(covStruct)
	
}

