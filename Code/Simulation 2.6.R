library(TSA)
data("Nile")

y <- Nile

# Obtain the MLE for the variances of the LLM Model
LLM <- StructTS(y, "level")
sigma_eps <- LLM$coef["epsilon"]
sigma_eta <- LLM$coef["level"]


## Kalman Filter
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


# Initialisation of parameters
a1 <- y[1]
P1 <- 1e7
a_keep <- P_keep <- NULL
for(i in 1:length(y)){
  
  at <- a1
  Pt <- P1
  yt <- y[i]
  
  params <- c(yt, at, Pt, sigma_eps,sigma_eta)
  upd_params <- KF_LLM(params)                # apply recursion
  
  a1 <- upd_params[1]
  P1 <- upd_params[2]
  a_keep <- c(a_keep, a1)
  P_keep <- c(P_keep, P1)

}

## State Smoothing

# Backward recursion
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


rt <- Nt <- double(length(y))
a_hat_keep <- NULL
for(i in length(y):1){ # backward recursion from n to 1 
  
  a1 <- a_keep[i]
  P1 <- P_keep[i]
  
  at <- a1
  Pt <- P1
  yt <- y[i]
  cur_rt <- rt[i]
  cur_Nt <- Nt[i]
  
  params <- c(yt, at, Pt,sigma_eps, sigma_eta, cur_rt,cur_Nt)
  upd_recursion <- SS_LLM(params) # apply State Smoothing

  rt[i-1] <- upd_recursion[1]
  Nt[i-1] <- upd_recursion[2]
  a_hat <- upd_recursion[3]
  
  a_hat_keep <- c(a_hat_keep, a_hat)

}

ahat_keep <- ts(a_hat_keep[length(y):1], start = 1871)
plot(y, bty = "l", ylab = "Nile")
lines(ahat_keep, type = "l", lwd = 2, col = "blue")
lines(ts(a_keep,start = 1871), type = "l", lwd = 2, col = "red")
legend("topright", lty = 1, col = c("red","blue"),lwd = 2, legend = c("Filtered State from Kalman Filter", "Smoothed State"), box.lty   = 0) 
## Disturbance Smoothing

# Smoothing recursion
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
  
eps_hat <- eta_hat <- double(length(y))
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
  
  # Apply smoothing error
  eps_hat[i] <- sigma_eps * upd_smootherr[1] 
  eta_hat[i] <- sigma_eta * rt[i] 
}


## Simulation 2.6
a1_plus <- y[1]
eps_plus <- eta_plus <- y_plus <- a_plus <- double(length(y))

# Generate the data
for(i in 1:length(y)){
  eps_plus[i] <- rnorm(1,0,sqrt(sigma_eps))
  eta_plus[i] <- rnorm(1,0,sqrt(sigma_eta))
  
  a_plus[i] <- a1_plus
  y_plus[i] <- a_plus[i] + eps_plus[i]
  a1_plus <- a_plus[i] + eta_plus[i]
  
}

### KF for generated Data
a1 <- a_plus[1]
P1 <- 1e7
a_keep_plus <- P_keep_plus <- NULL
for(i in 1:length(y)){
  
  at <- a1
  Pt <- P1
  yt <- y_plus[i]
  
  params <- c(yt, at, Pt, sigma_eps,sigma_eta)
  upd_params <- KF_LLM(params)
  
  a1 <- upd_params[1]
  P1 <- upd_params[2]
  a_keep_plus <- c(a_keep_plus, a1)
  P_keep_plus <- c(P_keep_plus, P1)
  
}

# Use KF to obtain the values of at and Pt in the generated data

### SS for generated Data

rt_plus <- Nt_plus <- double(length(y))
for(i in length(y):1){
  
  a1 <- a_keep_plus[i]
  P1 <- P_keep_plus[i]
  
  at <- a1
  Pt <- P1
  yt <- y_plus[i]
  cur_rt <- rt_plus[i]
  cur_Nt <- Nt_plus[i]
  
  params <- c(yt, at, Pt,sigma_eps, sigma_eta, cur_rt, cur_Nt)
  upd_recursion <- SS_LLM(params)
  
  rt_plus[i-1] <- upd_recursion[1]
  Nt_plus[i-1] <- upd_recursion[2]

}

# Apply SS to obtain the value of rt for the generated data

#### DS for generated Data

eps_plus_hat <- double(length(y))
for(i in 1:length(y)){
  a1 <- a_keep_plus[i]
  P1 <- P_keep_plus[i]
  
  at <- a1
  Pt <- P1
  yt <- y_plus[i]
  cur_rt <- rt_plus[i]
  cur_Nt <- Nt_plus[i]
  
  params <- c(yt, at, Pt,sigma_eps, sigma_eta, cur_rt, cur_Nt)
  upd_smootherr <- DS_LLM(params)
  
  eps_plus_hat[i] <- sigma_eps * upd_smootherr[1] 
}
# apply the actual smoother to the generated data

eps_tilde <- eps_plus - eps_plus_hat + eps_hat
a_tilde <- y - eps_tilde
eta_tilde <- diff(a_tilde)

## Plots
par(mfrow = c(2,2),
    mar = c(3,3,3,3))

# (i)
plot(ts(a_plus,start = 1871), type = "p",pch = 19, col = "grey", ylim = range(c(a_plus, a_hat_keep)), bty = "l", ylab = NULL, xlab = NULL)
lines(ts(a_hat_keep[length(y):1],start = 1871), lwd = 2)
title(main = "(i)", adj = 0)


# (ii)
plot(ts(a_tilde,start = 1871), type = "p",pch = 19, col = "grey", ylim = range(c(a_tilde, a_hat_keep)), bty = "l", ylab = NULL, xlab = NULL)
lines(ts(a_hat_keep[length(y):1],start = 1871), lwd = 2)
title(main = "(ii)", adj = 0)

# (iii)
plot(ts(eps_hat,start = 1871), type = "p",pch = 19, col = "grey", ylim = range(c(eps_hat, eps_tilde)), bty = "l", ylab = NULL, xlab = NULL)
lines(ts(eps_tilde,start = 1871), lwd = 2)
abline(h = 0)
title(main = "(iii)", adj = 0)

# (iv)
plot(ts(eta_tilde,start = 1871), type = "p",pch = 19, col = "grey", ylim = range(c(eta_hat, eta_tilde)), bty = "l", ylab = NULL, xlab = NULL)
lines(ts(eta_hat,start = 1871), lwd = 2)
abline(h = 0)
title(main = "(iv)", adj = 0)

par(mfrow = c(2,2),
    mar = c(2,2,2,2))

plot(ts(Nt, start = 1871))
title(main = expression("(i) - N"[t]), adj = 0)
plot(ts(Nt_plus, start = 1871))
title(main = expression("(ii) - N"[t]^`+`), adj = 0)
plot(ts(P_keep, start = 1871))
title(main = expression("(iii) - P"[t]), adj = 0)
plot(ts(P_keep_plus, start = 1871))
title(main = expression("(iv) - P"[t]^`+`), adj = 0)
