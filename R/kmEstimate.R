`kmEstimate` <-
function(model, envir) {
	
	X <- model@X
	y <- model@y
	F <- model@F

  if (model@case=="NoNugget") case <- "Default"
  if (model@case=="1Nugget") case <- "Nugget"
  if (model@case=="Nuggets") case <- "Noisy"
  
  if (model@case=="Nuggets") {
    if (model@covariance@nugget.flag & !model@covariance@nugget.estim) nugget <- model@covariance@nugget
    if (model@noise.flag) nugget <- model@noise.var
  }
  
  if (is.element(model@method, c("MLE", "PMLE"))){
    fn <- logLikFun
  } else if (model@method=="LOO") {
    fn <- leaveOneOutFun
  } 
	if (model@gr==TRUE) {
	  if (is.element(model@method, c("MLE", "PMLE"))){
	    gr <- logLikGrad
	  } else if (model@method=="LOO") {
	    gr <- leaveOneOutGrad
    } 
	} else {
		gr <- NULL
	}
	
		# initialization: starting points and boundaries
  
  model <- switch(case,
    Default = kmNoNugget.init(model),
    Nugget = km1Nugget.init(model),
    Noisy = kmNuggets.init(model))
  
	parinit <- model@parinit
	lp <- length(parinit)
  
  lower <- model@lower
  upper <- model@upper
	
		# printing
	control <- model@control
	if (control$trace!=0) {
		cat("\n")
		cat("optimisation start\n")
		cat("------------------\n")
    cat("* estimation method   :", model@method, "\n")
		cat("* optimisation method :", model@optim.method, "\n")
		cat("* analytical gradient :")
		if (is.null(gr)) {
			cat(" not used\n")}
		else cat(" used\n")
		cat("* trend model : "); print(model@trend.formula, showEnv=FALSE)
		cat("* covariance model : \n")
		cat("  - type : ",  model@covariance@name, "\n") 
		if (case=="Default") { 
      cat("  - nugget : NO\n")
      cat("  - parameters lower bounds : ", lower, "\n")
   	 cat("  - parameters upper bounds : ", upper, "\n")
 		 cat("  - best initial point among ", model@control$pop.size, " : ", parinit, "\n")
		} else if (case=="Nugget") {
      cat("  - nugget : unknown homogenous nugget effect \n")
      #if (!is.null(nugget)) cat("with initial value : ", nugget)
		  cat("  - parameters lower bounds : ", lower[1:(lp-1)], "\n")
 		  cat("  - parameters upper bounds : ", upper[1:(lp-1)], "\n")
		  cat("  - upper bound for alpha   : ", model@upper[lp],"\n")
		  cat("  - best initial point among ", model@control$pop.size, " :\n")
 		  cat("         coef.     : ", covparam2vect(model@covariance), "\n")
 		  cat("         alpha.    : ", parinit[lp], "\n")
		} else if (case=="Noisy") {  # technically: includes the "known nugget" case
      if (model@covariance@nugget.flag) {
        cat("  - nugget :", model@covariance@nugget, "\n")
      } else {
        cat("  - noise variances :\n")
        print(model@noise.var)
      }
		  cat("  - parameters lower bounds : ", lower[1:(lp-1)], "\n")
 		  cat("  - parameters upper bounds : ", upper[1:(lp-1)], "\n")
 		  cat("  - variance bounds : ", c(lower[lp], upper[lp]), "\n")
 		  cat("  - best initial point among ", model@control$pop.size, " :\n")
 		  cat("          coef.     : ", covparam2vect(model@covariance), "\n")
 		  cat("          variance  : ", model@covariance@sd2, "\n")
		}
    if (model@optim.method=="BFGS") cat("\n")     
	} # end printing
	
		# optimization
	
	if (model@optim.method=="BFGS") {
    BFGSargs <- c("trace", "parscale", "ndeps", "maxit", "abstol", "reltol", "REPORT", "lnm", "factr", "pgtol")
    commonNames <- intersect(BFGSargs, names(control))
    controlChecked <- control[commonNames]
    if (length(control$REPORT)==0) {
      controlChecked$REPORT <- 1
    }
    if (is.element(model@method, c("MLE", "PMLE"))){
      fnscale <- -1
    } else if (model@method=="LOO"){
      fnscale <- 1
    }
    forced <- list(fnscale=fnscale)
    
    controlChecked[names(forced)] <- forced
    
    o <- optim(par=parinit, fn=fn, gr=gr,
    		  	method = "L-BFGS-B", lower = lower, upper = upper,
    	  		control = controlChecked, hessian = FALSE, model, envir=envir)
  }
   
  if (model@optim.method=="gen") {       
    genoudArgs <- formals(genoud)
    commonNames <- intersect(names(genoudArgs), names(control))
    genoudArgs[commonNames] <- control[commonNames]
    if (length(control$print.level)==0) {
      genoudArgs$print.level <- control$trace
    } else {
      genoudArgs$print.level <- control$print.level
    }
    
    if (is.element(model@method, c("MLE", "PMLE"))){
      max.goal <- TRUE
    } else if (model@method=="LOO") {
      max.goal <- FALSE
    }
    forced <- list(fn=fn, nvars=length(parinit), max=max.goal, starting.values=parinit, 
            Domains=cbind(lower, upper), gr=gr, gradient.check=FALSE, boundary.enforcement=2, 
            hessian=TRUE, optim.method="L-BFGS-B", model=model, envir=envir)
       
    genoudArgs[names(forced)] <- forced
    genoudArgs$... <- NULL
       
    o <- do.call(genoud, genoudArgs)   
    	
	}
   
  model@logLik <- o$value
  
  
	if (model@method=="LOO"){
	  
	  model@covariance <- vect2covparam(model@covariance, o$par)
	  
	  errorsLOO <- envir$errorsLOO
	  sigma2LOO <- envir$sigma2LOO
    sigma2.hat <- mean(errorsLOO^2/sigma2LOO)
    sigma.hat <- sqrt(sigma2.hat)
	  model@covariance@sd2 <- sigma2.hat
	  
    R <- envir$R
	  T <- chol(R)
	  x <- backsolve(t(T), y, upper.tri=FALSE)		# x:=(T')^(-1)*y
	  M <- backsolve(t(T), F, upper.tri=FALSE)		# M:=(T')^(-1)*F
	  
	  if (identical(model@known.param, "Trend")) {
	    z <- x - M %*% model@trend.coef
	  } else {
	    l <- lm(x ~ M-1)
	    beta.hat <- as.numeric(l$coef)
#       z <- backsolve(t(T), y - F%*%beta.hat, upper.tri=FALSE)
	    model@trend.coef <-beta.hat
	    Q <- qr.Q(qr(M))
	    H <- Q %*% t(Q)
	    z <- x - H %*% x
	  }
	  
	  model@T <- T*sigma.hat
	  model@z <- as.numeric(z/sigma.hat)
	  model@M <- M/sigma.hat
	  
    return(model)   # Fin cas LOO
	} else {
	
  # beta.hat ?
  T <- envir$T
  z <- envir$z
  
	x <- backsolve(t(T), y, upper.tri=FALSE)		# x:=(T')^(-1)*y
	M <- backsolve(t(T), F, upper.tri=FALSE)		# M:=(T')^(-1)*F
	
	if (!identical(model@known.param, "Trend")) {
		l <- lm(x ~ M-1)
		beta.hat <- as.numeric(l$coef)
		model@trend.coef <-beta.hat
	}
	  
  if (case=="Default") {
    sigma2.hat <- envir$sigma2.hat
    sigma2.hat <- as.numeric(sigma2.hat)
	  sigma.hat <- sqrt(sigma2.hat)
	  model@T <- T*sigma.hat
	  model@z <- as.numeric(z/sigma.hat)
	  model@M <- M/sigma.hat
	  model@covariance@sd2 <- as.numeric(sigma2.hat)		
	  model@covariance <- vect2covparam(model@covariance, o$par)  
  } else if (case=="Nugget") {
    v <- envir$v 
    v <- as.numeric(v)
	  s <- sqrt(v)
	  model@T <- T * s
	  model@z <- as.numeric(z / s)
	  model@M <- M / s
		
    param <- o$par
    lp <- length(param)
    alpha <- param[lp]
    model@covariance@sd2 <- alpha*v
	  model@covariance@nugget <- (1-alpha)*v
	  model@covariance <- vect2covparam(model@covariance, param[1:(lp-1)])
    
    model@lower <- model@lower[1:(lp-1)]    
    model@upper <- model@upper[1:(lp-1)]
	  model@parinit <- model@parinit[1:(lp-1)]    
  }
  if (case=="Noisy") {
    model@T <- T
    model@z <- as.numeric(z)
    model@M <- M
    
    param <- o$par
    lp <- length(param)
	  model@covariance@sd2 <- param[lp]
	  model@covariance <- vect2covparam(model@covariance, param[1:(lp-1)])
    
    model@lower <- model@lower[1:(lp-1)]    
	  model@upper <- model@upper[1:(lp-1)]
	  model@parinit <- model@parinit[1:(lp-1)]
  }
	
	return(model)   # fin cas MLE, PMLE
	}
}