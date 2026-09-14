# plots
plot_eps = function(fname) {
	postscript(file=paste(fname,".eps",sep=""), onefile=FALSE, horizontal=FALSE, width=4, height=4, paper="special", family="Times")
	# cairo_ps(file=paste(fname,".eps",sep=""), onefile=FALSE, width=4, height=4, family="Times")
	# Trim off excess margin space (bottom, left, top, right)
	par(mar=c(3, 2.9, 0.2, 0.7))
	# Trim off excess outer margin space (bottom, left, top, right)
	par(oma=c(0,0,0,0))
	# Trim off excess space for label and ticks (label, ticks, line)
	par(mgp=c(1.8,0.6,0))
	# lty: line styles (1=solid, 2=dash, 3=dot, 4=dash-dot)
	# lab: (# of x-ticks, # of y-ticks, len of ticks), approximately
	# lwd: line-width
	# cex.lab: fontsize scaling-factor for labels
}
# plot(xxyyzz, xlab="x-label", ylab="y-label", xlim=c(0, 120), ylim=c(0, 50), lty=1:4, lab=c(10, 7, 5), lwd=2, cex.lab=1.3)
# # cex: fontsize scaling-factor for legends
# legend("topright", c("Legend 1", "Legend 2", "Legend 3", "Legend 4"), lty=1:4, lwd=2, cex=1.05)
# dev.off()




# plot incremental counts
plot_inc = function(dN, del) {

    plot_eps("inc")
    par(mar=c(2.9, 1.8, 1.5, 0.6))

	N <- length(dN)
    plot(del*(1:N), dN, type="l", cex.axis=1.5, cex.lab=1.6, cex.main=1.8, xlab="Time (s)", ylab="", main="Incremental Counts")
	
    dev.off()
}


# plot covariate data
plot_EC = function(w, br='FD', mt) {

    plot_eps("EC")
    par(mar=c(1.8, 1.8, 1.5, 0.9))

    hist(w, breaks=br, cex.axis=1.5, cex.lab=1.6, cex.main=1.8, xlab="", ylab="", main=mt)
	
	box()
    dev.off()
}



# plot stochastic intensity
plot_si = function(Tr, lam, fname="si_sim") {

    plot_eps(fname)
    par(mar=c(2.9, 1.8, 1.7, 0.7))

    plot(Tr, lam, type="l", cex.axis=1.5, cex.lab=1.6, cex.main=1.8, xlab="Time (s)", ylab="", main="Stochastic Intensity")
	
    dev.off()
}



plot_L = function(L) {

	plot_eps("LR")
	par(mar=c(2.8, 1.8, 1.5, .3)) 

	plot(L, type="l", cex=0.6, pch=15, lty=1, lwd=2, xlab=("Iteration #"), ylab="", , main="log-Likelihood", cex.axis=1.5, cex.lab=1.6, cex.main=1.8)

	dev.off()
}



# scalar QQ-plot
plot_QQ = function(pp, fname, title="Q-Q Plot") {

	plot_eps(fname)
	par(mar = c(3, 3, 1.5, 1.2))

	s <- seq(0.01, 0.99, 0.01)
	len <- length(s)

	# simulate unit rate Poisson process
	poiss <- sim_pois(pp[length(pp)], 1)		
	qp <- quantile(poiss, probs = s)
	qpp <- quantile(pp, probs = s) 
	
	# browser()
	
	plot(qp, qpp, type="l", lwd=1, cex.axis=1.5, cex.lab=1.6, cex.main=1.8, xlab="Theoretical Quantiles", ylab="Sample Quantiles", main=title)
	lines(0:tail(poiss,1), 0:tail(poiss,1), type="l", lty=2, lwd=1)  
	
	dev.off()
}


# null distribution QQ-plot
plot_nullQQ = function(S, df, fname, title="Q-Q Plot") {

	plot_eps(fname)
	par(mar = c(3, 3, 1.5, 1.2))

	s <- seq(0.01, 0.99, 0.01)
	len <- length(s)

	# chi-squared (degree of freedom df) distribution 
	rv <- rchisq(length(S), df)		
	qrv <- quantile(rv, probs = s)
	qLM <- quantile(S, probs = s) 
	
	plot(qrv, qLM, type="l", lwd=1, cex.axis=1.5, cex.lab=1.6, cex.main=1.8, xlab="Theoretical Quantiles", ylab="Sample Quantiles", main=title)
	lines(0:tail(qLM,1), 0:tail(qLM,1), type="l", lty=2, lwd=1)  
	
	dev.off()
}



# plot LM histogram
plot_LM_hist = function(LMo,LMa,mt) {
  
	plot_eps("LM_hist")
	par(mar=c(1.8, 1.8, 1.5, 0.3))

	hist(LMo, breaks=12, xlim=c(min(min(LMo),min(LMa)),max(max(LMo),max(LMa))), col=rgb(1,0,0,.5), cex.axis=1.5, cex.lab=1.6, cex.main=1.8, xlab="", ylab="", main=mt)
	hist(LMa, breaks=40, col=rgb(0,0,1,.5), add=T)
	legend("topright", c("null", "alternative"), col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), cex=1.6, , pt.cex=2, pch=15)
	
	box()	
	dev.off()
}



# HIR plot
plot_hir = function(u, h) {
  
	plot_eps("hir")
	par(mar=c(2.9, 1.8, 1.5, 0.2))
  
	matplot(u, h, type="l", lty=1, lwd=c(1,2), col="black", ylim=c(0,max(h)), cex.axis=1.5, cex.lab=1.6, cex.main=1.8, xlab="Time (s)", ylab="", main="mHIR")
	legend("topright", c(TeX(r'(true)'), TeX(r'(estimate)')), cex=1.6, lty=1, lwd=c(1,2))

	dev.off()	
}



plot_roc = function(alf, S1, S2) {
	plot_eps("roc")
	par(mar=c(2.8, 3, 1.5, 0.2))
	
	plot(c(0,1), c(0,1), type="l", lty=3, lwd=1, cex.axis=1.5, cex.lab=1.6, cex.main=1.8, xlab="False Positive Rate", ylab="True Positive Rate", main="ROC")
	lines(alf, S1, type="l", lty=1, lwd=2)
	lines(alf, S2, type="l", lty=1, lwd=1)
	legend("bottomright", c("mHIR","HLE-MM"), cex=1.6, lty=1, lwd=c(2,1))
  
	dev.off()
}



