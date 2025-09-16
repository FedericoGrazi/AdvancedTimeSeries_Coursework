library(TSA)
data("Nile")

y <- Nile

StructTS(y, "level")


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
  
  names(rt_1) <- names(Nt_1) <- NULL
  upd_params <- c("rt" = rt_1,"Nt" = Nt_1)
  
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
SCORE_LLM <- function(theta){

  a1 <- 0
  P1 <- 1e7
  
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
  
  score_vect <- c(
    "score eps" = sigma_eps * sum( ut^2 - Dt),
    "score eta" = sigma_eta * sum( rt^2 - Nt)
  )
  
  return(score_vect)
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


theta_init <- c(700,700)
theta_1st <- theta_init - SCORE_LLM(theta_init)
theta_old <- theta_1st
theta_old_1 <- theta_init

i = 0
matrixTheta <- matrix(c(theta_init, theta_1st), ncol = 2, byrow = T)
crit <- c(1000,1000)
logik <- NULL
while(!all((abs(diag(crit))<c(0.001,0.0001))==T)){
  i = i+1
  theta_old_1 <- as.numeric(matrixTheta[i,])
  theta_old <- as.numeric(matrixTheta[i+1,])
  crit <- diag(SCORE_LLM(theta_old)-SCORE_LLM(theta_old_1))
  
  
  diff_val <- diag(theta_old-theta_old_1)
  diff_score <- solve(diag(SCORE_LLM(theta_old)-SCORE_LLM(theta_old_1)))
  cur_score <- as.matrix(SCORE_LLM(theta_old))
  
  theta_new <- theta_old - t(cur_score) %*% (diff_val %*% diff_score)
  
  loglik <- c(loglik,extract_loglik(theta_new))
  
  matrixTheta <- rbind(matrixTheta, theta_new)
  if(i %% 50 == 0) print(paste0("Iteration: ", i))  # print or paste the iteration number
  
  
}
matrixTheta[nrow(matrixTheta),]; matrixTheta[nrow(matrixTheta),2]/matrixTheta[nrow(matrixTheta),1];print(paste0("iteration: ",i))



# RESULTS 
sel <- round(seq(1,i,length.out = 6))
scores <- matrix(0,ncol = 2, nrow = 6)
for(i in 1:6) scores[i,] <- SCORE_LLM(matrixTheta[i,])

data.frame(
  iteration = sel,
  q = matrixTheta[sel,2]/matrixTheta[sel,1],
  psi = log(matrixTheta[sel,2]/matrixTheta[sel,1]),
  loglik = -loglik[sel]
)

library(plotly)
N <- 150
grid.eps  <- seq(800,20000,length = N)
grid.eta <- seq(950,1800, length = N)
for(i in 1:N){
  for(j in 1:N){
    th <- c(grid.eps[i], grid.eta[j])
    z[i,j] <- extract_loglik(th)
  }
}
plot_ly(x = grid.eps, y = grid.eta, z = z, type = 'surface')



llik_opt <- optim(c(1000,1000), extract_loglik, 
      method = "L-BFGS-B", 
      lower = c(0,0), 
      control = list(trace = T))
