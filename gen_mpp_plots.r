#---------------
# Generate plots 
#---------------
library(latex2exp)


#*** incremental counts
nbins = 100						# number of time bins
sz <- ceiling(T / nbins)		# bin size
pp <- mpp[1,]
len <- length(pp)
N <- ceiling(T / sz)			# number of partitions

dN <- rep(0, N)					# incremental counts

for (r in 1:len) {

	ix <- round(as.numeric(pp[r]) / sz)
	dN[ix] = dN[ix] + 1
}	

plot_inc(dN, sz)


#*** plot event covariate data
plot_EC(mpp[2,], 25 ,"Event Covariate")



#*** stochastic intensity  
# lam <- c(pars$c, sapply(2:len, function(r) si(pp[1:(r-1)], mpp[2,1:(r-1)], pars, pp[r])))


# # plot intensity function
# plot_si(mpp[1,], lam, "si_sim")

