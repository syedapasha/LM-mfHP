#----------------------------------------------
# EM algorithm for marked finite Hawkes Process 
#----------------------------------------------

em_mfHP = function(p, pp, bet, T, eps) {

	n <- length(pp)
	
	L_arr <- NULL

	
	# compute xi(t) at event times
	xi <- cbind(rep(0,p), sapply(2:n, function(r) sapply(1:p, function(k) bet*sum(exp(-bet*(pp[r]-pp[1:(r-1)]))*((bet*(pp[r]-pp[1:(r-1)]))^(k-1)))/gamma(k))))

	# compute integral on [0,T] of xi(t)
	dT <- T - pp
	I_xi <- colSums(sapply(1:p, function(k) (1 - exp(-bet*dT)*.rowSums(sapply(0:(k-1), function(j) ((bet*dT)^j)/factorial(j)), n, k))))


	#*** EM algorithm ***#

	# starting values of parameter vector (c, a)
	cc <- runif(1, 0, 0.1)
	a <- rep(1/(p+1), p)

	
	lam <- cc + a %*% xi									# intensity function under null hypothesis
	
	L <- sum(log(lam)) - (cc*T + sum(a * I_xi))				# log-likelihood (ratio)
	
	res <- NULL
	err <- 1												# relative error

	while (err > eps) {
		
		#*** multiplicative updates ***#

		# c-update 
		cc <- cc * sum(1 / lam) / T
		
		# a-update 		
		a <- a * sapply(1:p, function(k) sum(xi[k,]/lam)) / I_xi

		lam <- cc + a %*% xi								# update intensity function 

		
		# compute log-likelihood (ratio) iterate
		L_old <- L
		L <- sum(log(lam)) - (cc*T + sum(a * I_xi))
	
		err <- abs((L - L_old) / L)
		
		res$c <- c(res$c, cc)
		res$a <- cbind(res$a, a)
		res$L <- c(res$L, L)								# log-likelihood estimates

	}
	
	return(res)	
}

