`kmNuggets.init` <-
function(model) {

	n <- nrow(model@X)
	parinit <- model@parinit

	if (model@covariance@nugget.flag & !model@covariance@nugget.estim) nugget.aux <- rep(model@covariance@nugget, n)
	if (model@noise.flag) nugget.aux <- model@noise.var
	
		# variance standard estimate (biased negatively)
	trend.estimate <- lm(model@y~model@F-1)
	random.part.estimate <- trend.estimate$residuals
	varinit.total <- var(random.part.estimate)
	varinit.standard <- varinit.total - mean(nugget.aux)
	if (varinit.standard <=1e-20) varinit.standard <- 1/2*varinit.total
	
		# variance estimate using the variogram (biased negatively, but less)
	x.dist <- dist(model@X)
	y.dist <- dist(random.part.estimate)
	I <- (x.dist > quantile(x.dist, 0.5))
	matrix.nugget.aux <- matrix(nugget.aux, n, n)
	matrix.sym.nugget.aux <- (matrix.nugget.aux + t(matrix.nugget.aux))/2
	nugget.aux.I <- matrix.sym.nugget.aux[I]
	varinit.vario.total <- 1/2*mean(y.dist[I]^2)
	varinit.vario <- varinit.vario.total - mean(nugget.aux.I)
	if (varinit.vario<=1e-20) varinit.vario <- 1/2*varinit.vario.total
	
		# final choice
	varinit <- (varinit.standard + varinit.vario)/2
	
	
		# boundaries    
	varinit.upper <- varinit.total - min(nugget.aux) 
	varinit.lower <- varinit.total - max(nugget.aux)
	if (varinit.upper<=1e-20) varinit.upper <- varinit.total
	if (varinit.lower<=1e-20) varinit.lower <- 1e-20
			
	varinit.vario.upper <- varinit.vario.total - min(nugget.aux.I) 
	varinit.vario.lower <- varinit.vario.total - max(nugget.aux.I)
	if (varinit.vario.upper<=1e-20) varinit.vario.upper <- varinit.vario.total
	if (varinit.vario.lower<=1e-20) varinit.vario.lower <- 1e-20
	
		# final choice
	varinit.lower <- 1/10*min(varinit.lower, varinit.vario.lower)
	varinit.upper <- 10*max(varinit.upper, varinit.vario.upper)
		
		
	lower <- model@lower
	upper <- model@upper
	ninit <- model@control$pop.size
	param.n <- model@covariance@param.n
	
	if (identical(parinit, numeric(0))) {
			# sample ninit design points, generated from uniform [lower, upper]
		matrixinit <- matrix(runif(ninit*param.n), param.n, ninit)
    if ((!is(model@covariance, "covAffineScaling")) & (!is(model@covariance, "covScaling"))) {
  	  matrixinit <- lower + matrixinit*(upper - lower)
    } else {
      matrixinit <- 1/upper + matrixinit*(1/lower - 1/upper)
      matrixinit <- 1/matrixinit
    }    
	} else matrixinit <- matrix(parinit, param.n, ninit) 
			
	varinit.sim <- runif(n=ninit, min=1/2*varinit, max=3/2*varinit)       
				
		# take the best point				 
	matrixinit <- rbind(matrixinit, varinit.sim)     
	logLikinit <- apply(matrixinit, 2, logLikFun, model)
	parinit <- matrixinit[, which.max(logLikinit)]     
		
		# result
	model@parinit <- as.numeric(parinit)
	
	lp <- length(parinit)
	var.init <- parinit[lp]
	parinit <- parinit[1:(lp-1)]
   	   	
   	model@covariance <- vect2covparam(parinit, model@covariance)
	model@covariance@sd2 <- var.init
	
	model@lower <- c(lower, varinit.lower)
	model@upper <- c(upper, varinit.upper)
  	
	return(model)
}

