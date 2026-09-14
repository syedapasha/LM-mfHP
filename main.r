#-----------------------------------------------------------------------
# Lagrange multiplier test for event covariates in marked Hawkes process 
#-----------------------------------------------------------------------

rm(list = ls(all = TRUE))
# include packages and source files

# library(doParallel)
# library(foreach)

source("data_gen.r")
source("em.r")
source("residual.r")
source("plots.r")


# simulation setting 
p <- c(4)									# number of Laguerre/Erlang basis terms
m <- c(3)									# number of polynomial basis terms


cc <- c(0.1)								# background rate
rho <- c(.9)								# memory contribution

T <- c(5000)								# simulation time

eps <- c(1e-9)								# EM stopping threshold 


# mHIR specification
bet <- c(1)									# reciprocal of time constant
a <- sapply(1:p, function(k) 1.2^k)
B <- matrix(0,p,m)


# *** uncomment for simulation under alternative hypothesis ***
# ---------------------------
# B[,1] <- (1/1.3)^(1:p)
# B[,2] <- (1/1.2)^(1:p)
# B[,3] <- (1/1.1)^(1:p)
# B <- .4*B
# ---------------------------


# moments of event covariates 
# ---------------------------
# *** for Exp(4) distribution ***
lam_exp <- 4
mu <- sapply(1:m, function(el) factorial(2*el)/(lam_exp^(2*el)))

# *** for standard normal distribution ***
# m_even <- c(1,3,15,105,945,10395)
# mu <- m_even[1:m]


# ensure Hawkes stability condition
Hv <- sum(a) + sum(B %*% mu)
a  <- rho * a / Hv
B  <- rho * B / Hv



pars <- NULL
pars$c <- cc
pars$a <- a
pars$B <- B
pars$bet <- bet
pars$lam <- lam_exp


#-----------------------------------------
#*** simulation and residual analysis ***#
#-----------------------------------------
# cat("\n Starting simulation...")
mpp <- sim_mfHP(T, pars)

# source("gen_mpp_plots.r")



#------------------------------------------------------------------
#*** model fitting under null hypothesis and residual analysis ***#
#------------------------------------------------------------------

# p <- c(4)									# number of Laguerre/Erlang basis for estimation
# m <- c(3)									# number of polynomial basis for estimation



cat("\n Starting EM algorithm...")
res <- em_mfHP(p, mpp[1,], bet, T, eps)

nIter <- length(res$c)						# number of iterations

# maximum likelihood estimates
em <- NULL
em$c <- res$c[nIter]
em$a <- res$a[,nIter]
em$B <- matrix(rep(0,p*m),p,m)    # for residual analysis



# residual analysis for goodness-of-fit
#--------------------------------------
tau <- residual_mfHP(mpp, em, bet)


# QQ-plot 
plot_QQ(tau, "qq_fit")



# plot log-likelihood (ratio) iterates
#-------------------------------------
plot_L(res$L[1:200])



#---------------------------------
#*** Lagrange multiplier test ***#
#---------------------------------
pp <- mpp[1,]
v  <- mpp[2,]

n <- length(pp)


# vec(X(t)) 
vec_X <- cbind(rep(0,p*m), sapply(2:n, function(i) sapply(1:m, function(el) sapply(1:p, function(k) bet*sum((v[1:(i-1)]^(2*el))*exp(-bet*(pp[i]-pp[1:(i-1)]))*((bet*(pp[i]-pp[1:(i-1)]))^(k-1)))/gamma(k) ))))


# compute intensity function from model fit under null
# ----------------------------------------------------
xi <- cbind(rep(0,p), sapply(2:n, function(r) sapply(1:p, function(k) bet*sum(exp(-bet*(pp[r]-pp[1:(r-1)]))*((bet*(pp[r]-pp[1:(r-1)]))^(k-1)))/gamma(k))))

lam_o <- em$c + em$a %*% xi



# compute L_pi 
# ------------
# summand in first term 
X <- sapply(1:(p*m), function(r) vec_X[r,]/lam_o)

# second term 
dT <- T - pp
c_T <- c(colSums(sapply(1:m, function(el) sapply(1:p, function(k) (v^(2*el))*(1 - exp(-bet*dT)*.rowSums(sapply(0:(k-1), function(j) ((bet*dT)^j)/factorial(j)), n, k))), simplify="array")))

L_pi <- colSums(X) - c_T


# compute Hessians
# ----------------
eta <- rbind(rep(1,n), xi)
Z <- sapply(1:(p+1), function(r) eta[r,]/lam_o)

H_pi_pi <- t(X) %*% X
H_pi_psi <- t(X) %*% Z
H_psi_psi <- t(Z) %*% Z



# compute QR decompositions of Z and X
# ------------------------------------
qrZ <- qr(Z)
Qz <- qr.Q(qrZ)
Rz <- qr.R(qrZ)

qrX <- qr(X)
Qx <- qr.Q(qrX)
Rx <- qr.R(qrX)



# compute LM
# ----------
w_T <- colSums(Qx) - solve(t(Rx), c_T)
Qxz <- t(Qx) %*% Qz
W <- diag(p*m) - Qxz %*% t(Qxz)

LM <- c(t(w_T) %*% solve(W, w_T))

print(LM)



# h(u,v) for simulated and fitted marked Hawkes impulse response 
#---------------------------------------------------------------
u <- seq(0, 10, by=.1)
phi <- sapply(1:p, function(k) bet*exp(-bet*u)*((bet*u)^(k-1)) / gamma(k))
h <- phi %*% cbind(a, em$a)


# HIR plot
plot_hir(u, h)


# stochastic intensity based on MLE
#----------------------------------
# plot intensity function
plot_si(mpp[1,], lam_o, "si_fit")

