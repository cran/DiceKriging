`km1Nugget` <-
function(model, envir) {
	
	X <- model@X
	y <- model@y
	F <- model@F
	
	fn <- logLikFun
	if (model@gr==TRUE) {
		gr <- logLikGrad}
	else gr <- NULL
	
	
		# getting starting points and boundaries
	model <- km1Nugget.init(model)

	lower <- model@lower
	upper <- model@upper
	parinit <- model@parinit
	lp <- length(parinit)
	
		# printing
	if (control$trace!=0) {
		cat("\n")
		cat("optimisation start\n")
		cat("------------------\n")
		cat("* optimisation method :", model@optim.method, "\n")
		cat("* analytical gradient :")
		if (is.null(gr)) {
			cat(" not used\n")}
		else cat(" used\n")
		cat("* trend model : "); print(model@trend.formula)
		cat("* covariance model : \n")
		cat("  - type : ",  model@covariance@name, "\n")
	
		cat("  - nugget : unknown homogenous nugget effect ")
		#if (!is.null(nugget)) cat("with initial value : ", nugget)
		cat("\n")
		cat("  - parameters lower bounds : ", lower[1:(lp-1)], "\n")
 		cat("  - parameters upper bounds : ", upper[1:(lp-1)], "\n")
		cat("  - best initial point among ", model@control$pop.size, " :\n")
 		cat("         coef.     : ", covparam2vect(model@covariance), "\n")
 		cat("         alpha.    : ", parinit[lp], "\n")
 		cat("\n")
 	}

	
		# optimization
			
	if (model@optim.method=="BFGS") {	
		o <- optim(parinit, fn=fn, gr=gr,
    		  method = "L-BFGS-B", lower = lower, upper = upper,
    	  	  control = list(fnscale=-1, trace=control$trace, REPORT=1, maxit=200),  
    	  	  hessian = FALSE, model, envir=envir)
    }
   	
   	if (model@optim.method=="gen") {
   		domaine <- cbind(lower, upper)
   		o <- genoud(fn, nvars=lp, max=TRUE, 
   	 		  pop.size=control$pop.size, 
   	 		  max.generations=control$max.generations, 
   	 		  wait.generations=control$wait.generations,
              hard.generation.limit=TRUE, starting.values=parinit, MemoryMatrix=TRUE, 
              Domains=domaine, default.domains=10, solution.tolerance=0.001,
              gr=gr, boundary.enforcement=2, lexical=FALSE, gradient.check=FALSE, BFGS=TRUE,
              data.type.int=FALSE, hessian=TRUE, unif.seed=812821, int.seed=53058,
              print.level=control$trace, share.type=0, instance.number=0,
              output.path="stdout", output.append=FALSE, project.path=NULL,
              P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
              P9mix=NULL, 
              BFGSburnin=control$BFGSburnin, 
              BFGSfn=NULL, BFGShelp=NULL,
              cluster=FALSE, balance=FALSE, debug=FALSE, model=model, envir=envir)	 }
                 		      		    
	param <- o$par
	lp <- length(param)
#	model@covariance$coef <- param[1:(lp-2)]
	alpha <- param[lp]
	
	#model@covariance@sd2 <- param[lp-1]
	#model@covariance@nugget <- param[lp]
	
	# beta estimator
	get("T", envir=envir)
	get("z", envir=envir)
	
	x <- backsolve(t(T),y, upper.tri=FALSE)		# x:=(T')^(-1)*y
	M <- backsolve(t(T),F, upper.tri=FALSE)		# M:=(T')^(-1)*F
	
	if (!identical(model@known.param, "Trend")) {
		l <- lm(x ~ M-1)
		beta.hat <- as.numeric(l$coef)
		model@trend.coef <-beta.hat
	}
	
	get("v", envir=envir)
	v <- as.numeric(v)
	s <- sqrt(v)
	model@T <- T * s
	model@z <- as.numeric(z / s)
	model@M <- M / s
		
	model@covariance@sd2 <- alpha*v
	model@covariance@nugget <- (1-alpha)*v
	model@covariance <- vect2covparam(param[1:(lp-1)], model@covariance)
	
	model@logLik <- o$value
	
	model@lower <- model@lower[1:(lp-1)]    
	model@upper <- model@upper[1:(lp-1)]
	model@parinit <- model@parinit[1:(lp-1)]
	
	return(model)
}

