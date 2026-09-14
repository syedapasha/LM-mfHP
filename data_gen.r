#----------------------------------------
# simulate a marked finite Hawkes process
#----------------------------------------
 
sim_mfHP = function(T, pars) {

	cc <- pars$c
	a  <- pars$a
	B  <- pars$B
	
	bet    <- pars$bet
	lam_e  <- pars$lam

	p <- length(a)
	m <- dim(B)[2]

	
	pp <- NULL								# event times 
	v  <- NULL								# event covariates 

	#*** initialization ***
	lam <- cc								# intensity at time 0
	lam_s <- lam

	#*** first event ***
	s <- -log(runif(1)) / lam_s
	w <- -log(runif(1)) / lam_e
	
	pp <- c(pp, s)
	v <-  c(v, w)
	

	bProceed <- 1
	while (s<T & bProceed==1) {
	
		bProceed <- 0	
		
		#*** general routine ***
		#*** update maximum intensity ***		

		lam <- si(pp[pp<s], v[pp<s], pars, s)
		v_t <- tail(v,1)
		lam_s <- lam + (a[1] + sum(B[1,]*sapply(1:m, function(el) v_t^(2*el))))*bet
		
		while (s<T & bProceed==0) {

			#*** new event ***
			s <- s - log(runif(1)) / lam_s
			w <- -log(runif(1)) / lam_e

			#*** rejection test ***

			lam <- si(pp, v, pars, s)

			D <- runif(1)

			if (s<T & D <= lam/lam_s) {
			
				pp <- c(pp, s)
				v <-  c(v, w)
		
				bProceed <- 1
			}

			else lam_s = lam
		}
	}

	return(rbind(pp, v))
}



#-------------------------------
# calculate stochastic intensity
#-------------------------------

si = function(pp, v, pars, s) {

	cc <- pars$c
	a  <- pars$a
	B  <- pars$B

	bet  <- pars$bet

	p <- length(a)
	m <- dim(B)[2]

	dT <- s-pp
	xi <- sapply(1:p, function(k) bet*sum(exp(-bet*dT)*((bet*dT)^(k-1)))/gamma(k))

	X <- sapply(1:m, function(el) sapply(1:p, function(k) bet*sum((v^(2*el))*exp(-bet*dT)*((bet*dT)^(k-1)))/gamma(k)))
	
	return(cc + sum(a*xi) + sum(B*X))		

}



#---------------------------
# simulate a Poisson process
#---------------------------

sim_pois = function(T, lam) {
 
	pp <- 0
	len <- 1

	while (pp[len] < T) {
		U <- runif(1)
		pp[len+1] <- pp[len] - log(U)/lam
		len <- len + 1
	}
	
	pp <- pp[2:(length(pp)-1)]
	
	return(pp)	
}
