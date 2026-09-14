# construct ROC curve
# -------------------
# FPR values
m <- 50										# number of FPR grid points
alf <- c(1:m) / m							# FPR values


load("res2.RData")
LM_n <- sort(LM_n)
LM_a <- sort(LM_a)
nR <- 1000

r <- nR * (1 - alf) + 1						# index
t <- LM_n[r]								# threshold

S <- 1 - sapply(t, function(t_j) which(LM_a > t_j, arr.ind = TRUE)[1]) / nR


load("resMM2.RData")
LM_MM_n <- sort(LM_MM_n)
LM_MM_a <- sort(LM_MM_a)
nR <- 1000

r <- nR * (1 - alf) + 1						# index
t <- LM_MM_n[r]								# threshold

S_MM <- 1 - sapply(t, function(t_j) which(LM_MM_a > t_j, arr.ind = TRUE)[1]) / nR


plot_roc(alf, S, S_MM)
