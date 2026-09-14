#-------------------------------------------------------
# time transformation for marked Hawkes-Laguerre process
#-------------------------------------------------------

residual_mfHP <- function(mpp, pars, bet) {

	cc <- pars$c
	a  <- pars$a
	B  <- pars$B

	p <- length(a)
	m <- dim(B)[2]

	n <- dim(mpp)[2]

	pp <- mpp[1,]
	v  <- mpp[2,]
	
	
	res1 <- sapply(2:n, function(i) sapply(1:p, function(k) sum(1 - exp(-bet*(pp[i]-pp[1:(i-1)]))*.rowSums(sapply(0:(k-1), function(j) ((bet*(pp[i]-pp[1:(i-1)]))^j)/factorial(j)), i-1, k))))

	res2 <- sapply(2:n, function(i) sapply(1:m, function(el) sapply(1:p, function(k) sum((v[1:(i-1)]^(2*el))*(1 - exp(-bet*(pp[i]-pp[1:(i-1)]))*.rowSums(sapply(0:(k-1), function(j) ((bet*(pp[i]-pp[1:(i-1)]))^j)/factorial(j)), i-1, k))))), simplify="array")

	
	tau <- c(cc*pp[1], cc*pp[2:n] + a %*% res1 + sapply(1:(n-1), function(i) sum(B * res2[,,i]))) 


	return(tau)
}
