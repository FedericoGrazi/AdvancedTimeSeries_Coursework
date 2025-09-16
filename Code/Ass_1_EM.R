library(TSA)
data("Nile")

y <- Nile


check_both_true <- function(result) {
  if (length(result) == 2 && all(result == TRUE)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

KF_LLM <- function(params){
  yt <- params[1]
  at <- params[2]
  Pt <- params[3]
  sigma_eps <- params[4]
  sigma_eta <- params[5]
  
  vt <- yt-at
  Ft <- Pt + sigma_eps
  Kt <- Pt / Ft
  
  at_upd <- at + Kt * vt
  Pt_upd <- Pt * (1 - Kt) + sigma_eta
  
  names(at_upd) <- names(Pt_upd) <- NULL
  upd_params <- c("at" = at_upd, "Pt" = Pt_upd)
  
  return(upd_params)
}

SS_LLM <- function(params){
  
  yt <- params[1]
  at <- params[2]
  Pt <- params[3]
  sigma_eps <- params[4]
  sigma_eta <- params[5]
  rt <- params[6]
  Nt <- params[7]
  
  vt <- yt-at
  Ft <- Pt + sigma_eps
  Kt <- Pt / Ft
  Lt <- 1-Kt

  rt_1 <- vt / Ft + Lt * rt
  Nt_1 <- Ft^-1 + Lt^2 * Nt
  a_hat <- at + Pt * rt_1
  
  names(rt_1) <- names(Nt_1) <- NULL
  upd_params <- c("rt" = rt_1,"Nt" = Nt_1, "at" = a_hat)
  
  
  return(upd_params)
}

DS_LLM <- function(params){
  yt <- params[1]
  at <- params[2]
  Pt <- params[3]
  sigma_eps <- params[4]
  sigma_eta <- params[5]
  rt <- params[6]
  Nt <- params[7]
  
  vt <- yt-at
  Ft <- Pt + sigma_eps
  Kt <- Pt / Ft
  
  
  Dt <- Ft^-1 + Kt^2 * Nt
  ut <- Ft^-1 * vt - Kt * rt  # Obtain value of u_t
  names(Dt) <- names(ut) <- NULL
  upd_params <- c("ut" = ut, "Dt" = Dt)
  
  return(upd_params)
}

############

# theta[1] = sigma epsilon
# theta[2] = sigma eta
EM_LLM <- function(theta){
  
  a1 <- 0
  P1 <- 1e7
  n <- length(y)
  
  sigma_eps <- theta[1]
  sigma_eta <- theta[2]
  
  a_keep <- P_keep <- NULL
  for(i in 1:length(y)){
    
    at <- a1
    Pt <- P1
    yt <- y[i]
    
    par_KF <- c(yt, at, Pt, sigma_eps,sigma_eta)
    upd_KF <- KF_LLM(par_KF)                # apply recursion
    
    a1 <- upd_KF[1]
    P1 <- upd_KF[2]
    
    a_keep <- c(a_keep, a1)
    P_keep <- c(P_keep, P1)
  }
  
  names(a_keep) <- names(P_keep) <- NULL
  
  rt <- Nt <- double(length(y))
  for(i in length(y):1){
    
    a1 <- a_keep[i]
    P1 <- P_keep[i]
    
    at <- a1
    Pt <- P1
    yt <- y[i]
    cur_rt <- rt[i]
    cur_Nt <- Nt[i]
    
    params <- c(yt, at, Pt,sigma_eps, sigma_eta, cur_rt, cur_Nt)
    upd_recursion <- SS_LLM(params)
    
    rt[i-1] <- upd_recursion[1]
    Nt[i-1] <- upd_recursion[2]
  }  
  
  
  
  ut <- Dt <- double(length(y))
  for(i in 1:length(y)){
    a1 <- a_keep[i]
    P1 <- P_keep[i]
    
    at <- a1
    Pt <- P1
    yt <- y[i]
    cur_rt <- rt[i]
    cur_Nt <- Nt[i]
    
    params <- c(yt, at, Pt,sigma_eps, sigma_eta, cur_rt, cur_Nt)
    upd_smootherr <- DS_LLM(params)
    
    ut[i] <- upd_smootherr[1]
    Dt[i] <- upd_smootherr[2]
  }
  
  em_vect <- c(
    "upt eps" = sigma_eps + 1/n * sigma_eps^2 * sum(ut^2 - Dt),
    "upd eta" = sigma_eta + 1/(n-1) * sigma_eta^2 * sum(rt^2 - Nt)
  )
  
  return(em_vect)
}

extract_loglik <- function(theta){
  n <- length(y)
  a1 <- 0
  P1 <- 1e7
  sigma_eps <- theta[1]
  sigma_eta <- theta[2]
  vt_keep <- Ft_keep <- NULL
  for(i in 1:length(y)){
    
    at <- a1
    Pt <- P1
    yt <- y[i]
    
    vt <- yt-at
    Ft <- Pt + sigma_eps
    Kt <- Pt / Ft
    
    a1 <- at + Kt * vt
    P1 <- Pt * (1 - Kt) + sigma_eta
    
    vt_keep <- c(vt_keep, vt)
    Ft_keep <- c(Ft_keep, Ft)
  }
  
  constant <- -n/2 * log(2*pi)
  logsum <- -.5 * sum( log(Ft_keep[2:n]) + vt_keep[2:n]^2/Ft_keep[2:n])
  
  loglik <- -( constant  + logsum)
  return(loglik)
}

theta_init <- c(1000,1000)
matrixTheta <- matrix(c(theta_init), ncol = 2, byrow = T)
crit <- 1000
loglik <- extract_loglik(theta_init)

i = 0
while(crit>0.001){
  i = i+1
  theta_old <- as.numeric(matrixTheta[i,])
  
  theta_new <- EM_LLM(theta_old)
  loglik <- c(loglik,extract_loglik(theta_new))
  
  matrixTheta <- rbind(matrixTheta, theta_new)
  crit <- abs(loglik[i+1]-loglik[i])
  print(paste0("it: ",i))
}
matrixTheta[nrow(matrixTheta),]; matrixTheta[nrow(matrixTheta),2]/matrixTheta[nrow(matrixTheta),1];print(paste0("iteration: ",i))

# RESULTS 
sel <- round(seq(1,i,length.out = 6))

df <- data.frame(
  iteration = sel,
  q = matrixTheta[sel,2]/matrixTheta[sel,1],
  psi = log(matrixTheta[sel,2]/matrixTheta[sel,1]),
  loglik = -loglik[sel]
)

par(mfrow = c(1,2))
plot(ts(loglik,start = 1871), type = "l", ylab = "diff log Ld" )
plot(diff(loglik), type = "l", ylab = "diff log Ld" ,xlab = "Iteration")
abline(h = 0.001, col = "red", lty = 2, lwd = 2)
legend(x = 120,y = -15, legend = "ϵ", col = "red", lty = 2, lwd = 2, box.lty  =  0, cex = 1.3)


# KF

a1 <- 0
P1 <- 1e7
n <- length(y)

sigma_eps <- matrixTheta[nrow(matrixTheta),1]
sigma_eta <- matrixTheta[nrow(matrixTheta),2]
sigma_eps <- 15099
sigma_eta <- 1469


a_keep <- P_keep <- v_keep  <- NULL
for(i in 1:length(y)){
  
  at <- a1
  Pt <- P1
  yt <- y[i]
  
  par_KF <- c(yt, at, Pt, sigma_eps,sigma_eta)
  upd_KF <- KF_LLM(par_KF)                # apply recursion
  
  a1 <- upd_KF[1]
  P1 <- upd_KF[2]
  
  a_keep <- c(a_keep, a1)
  P_keep <- c(P_keep, P1)
  v_keep<- c(v_keep,yt-at)
}

par(mfrow = c(1,1))
plot(y, bty = "l", type = "p", pch = 19, col = "grey")
lines(ts(a_keep, start = 1871), col = "red", lwd = 2)
lines(ts(a_keep + 1.96 * sqrt(P_keep), start = 1871), col = "red", lty = 2)
lines(ts(a_keep - 1.96 * sqrt(P_keep), start = 1871), col = "red", lty = 2)

names(v_keep) <- 1871:1970
car::qqPlot(v_keep, ylab = "EM estimates for ε")

