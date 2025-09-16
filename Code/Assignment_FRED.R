data <- read.csv(file.choose())

ts <- ts(data[1:292,2], start = c(1947,2), frequency = 4)
ts


# TS
par(mar = c(4.5,5,2,4))
plot(ts,bty = "l", ylab = " Real Gross Domestic Product \n (bilions of $)")
dev.off()
acf(ts,lag.max = 292,main = "Time Series")

forecast::auto.arima(ts, max.D = 0, max.P = 0, max.Q = 0)


# What do Local Linear Trend Model look like with these parameters?
fit <- StructTS(ts, type = "trend")
plot(ts)
lines(fit$fitted[,1], col = "red", lty = 2)
plot(fit$fitted[,2], col = "blue")
params <- fit$coef

sim_h2 <- function(n,sigma_eps,sigma_eta,sigma_csi,b = 4, plot = c(1,2)){
  n = n
  
  # epsilon
  sigma_eps <- sigma_eps
  eps <- sigma_eps * rnorm(n)
  
  # eta
  sigma_eta <- sigma_eta
  eta <- sigma_eta * rnorm(n)
  
  # csi
  sigma_csi <- sigma_csi
  csi <- sigma_csi * rnorm(n)
  
  
  m = 3
  b = b
  
  mu <- beta <- double(n)
  mu[1] <- rnorm(1,0,m)
  beta[1] <- rnorm(1,0,b)
  y = 0
  
  for(t in 1:(n-2)){
    beta[t+1] = beta[t] + csi[t]   
    mu[t+2] = beta[t+1] + mu[t+1] + eta[t+1]     
    y[t+2] =  mu[t+2] + eps[t+2]  
  }
  
  if(plot == 1){
    ts.plot(y)
    lines(mu,col="red")
    # legend("topright", col = c("black", "red"), lty = c(1,1), legend = c("Observed Component (y)", "Unobserved Component (μ)"))
  }else{
    par(mfrow= c(3,1), mar = c(1,1,1,1))
    ts.plot(y,col="black", xlab = NULL)
    legend("topleft", lty = 1,col = "black" , legend = "y", box.lwd = 0, box.lty = 1)
    ts.plot(mu,col="red", xlab = NULL)
    legend("topleft", lty = 1,col = "red" , legend = "μ", box.lwd = 0, box.lty = 1)
    ts.plot(beta,col="blue", xlab = NULL)
    legend("topleft", lty = 1,col = "blue" , legend = "β", box.lwd = 0, box.lty = 1)
    par(mfrow= c(1,1))
  }
}
sim_h2(length(ts), sigma_eps = params[3], sigma_eta = params[1], sigma_csi = params[2], plot = 1 )
sim_h2(length(ts), sigma_eps = params[3], sigma_eta = params[1], sigma_csi = params[2], plot = 2 )
dev.off()

# DIFF 1
diff.ts <- diff(ts)
par(mfrow= c(3,1), mar = c(3,4,1,3))
plot(diff.ts)
acf(diff.ts, lag.max = 400,drop.lag.0 = T, main = "")
pacf(diff.ts, lag.max = 400,drop.lag.0 = T,main = "")

# DIFF 2
diff2.ts <- diff(ts,differences = 2)
plot(diff2.ts)
acf(diff2.ts, lag.max = 400,drop.lag.0 = T,main = "")
pacf(diff2.ts, lag.max = 400,drop.lag.0 = T,main = "")
par(mfrow = c(1,1))
# MA(1) behaviour 
# yt is a LLM


# ARIMA(2,1,1) on Δy
arima211 <- arima(diff.ts,order = c(2,1,1))
plot(diff.ts, type = "p",bty = "l", ylab = expression("Δy"[t]))
lines(diff.ts-residuals(arima211), col = "red", lwd = 2)
 
# ARIMA(0,1,1) on Δy
arima011 <- arima(diff.ts,order = c(0,1,1))
lines(diff.ts-residuals(arima011), col = "green", lwd = 2)

# Local Level Model on Δy
llm <- StructTS(diff.ts, type = "level")
lines(llm$fitted, col = "blue",lwd = 2)
legend("bottomright", legend = c("ARIMA(2,1,1)", "ARIMA(0,1,1)","LLM"), col = c("red","green","blue"), lwd = 2, box.lty  =  0, cex = 1)


# Local Linear Trend Model on y
lltm <- StructTS(ts, type = "trend")
lltm

# sigma_eps = 0 !!!!



## Forecasting 


a1 <- matrix(0,2,1)
P1 <- matrix(c(1e2,0,0,1e2), ncol = 2, nrow = 2, byrow = T)
Zt <- matrix(c(1,0), ncol = 2, nrow = 1)
Tt <- matrix(c(1,1,0,1), ncol = 2, nrow = 2, byrow = T)
sigma_eps <- lltm$coef[3]
sigma_eta <- lltm$coef[1]
sigma_zeta <- lltm$coef[2]
Qt <- matrix(c(sigma_eta,0,0,sigma_zeta), ncol = 2, byrow = T)

# KALMAN FILTER
mu_keep <- beta_keep <- NULL
for( i in 1:length(ts)){
  yt <- ts[i]
  
  Pt <- P1
  at <- a1
  
  vt <- yt-Zt %*% at
  
  Ft <- Zt %*% Pt %*% t(Zt) + sigma_eps
  Kt <- Tt %*% Pt %*% t(Zt) %*% solve(Ft)
  
  P1 <- Tt %*% Pt %*% t(Tt) + Qt - Kt %*% Ft %*% t(Kt)
  a1 <- Tt %*% at + Kt %*% vt
  
  mu_keep <- c(mu_keep, a1[1,1])
  beta_keep <- c(beta_keep, a1[2,1])
  
}

mu_keep <- ts(mu_keep,start = c(1947,2), frequency = 4)
beta_keep <- ts(beta_keep,start = c(1947,2), frequency = 4)

# ####
# par(mfrow = c(2,2))
# plot(ts, type = "l",bty = "l", ylab = expression("y"[t]))
# title(main = "StructTS results", adj = 0)
# lines(lltm$fitted[,1], col = "red")
# plot(lltm$fitted[,2], col = "blue", bty = "l", ylab = expression("β"[t]))
# 
# plot(ts, type = "l",bty = "l", ylab = expression("y"[t]))
# title(main = "Kalman Filter results", adj = 0)
# lines(mu_keep, col = "red")
# plot(beta_keep, col = "blue", bty = "l",ylab = expression("β"[t]))
# par(mfrow = c(1,1))
# ####

# Kalman Filter for forecasting
y_pred <- P_pred <-  NULL
ynew <- mu_keep[length(mu_keep)-1]
j = 0
Kt <- matrix(c(0,0),nrow = 2)
while(j < 20){
  j = j+1
  yt1ahead <- ynew
  
  Pt <- P1
  at <- a1
  
  vt <- yt1ahead - Zt %*% at
  
  Ft <- Zt %*% Pt %*% t(Zt) + 200

  P1 <- Tt %*% Pt %*% t(Tt) + Qt - Kt %*% Ft %*% t(Kt)
  a1 <- Tt %*% at + Kt %*% vt
  
  ynew <- Zt %*% a1
  y_pred <- c(y_pred, ynew)
  P_pred <- c(P_pred, P1[1,1])
  
}
par(mfrow = c(1,1),mar = c(4,4,4,4))
plot(ts, type = "l",bty = "l", ylab = expression("y"[t]), xlim = c(2000, 2025) ,ylim = range(c(ts[212],y_pred)))
lines(ts(y_pred, start= c(2020,2), frequency = 4 ), col ="red", lwd = 2)
lines(ts(y_pred + 1.96 * sqrt(P_pred), start= c(2020,2), frequency = 4 ), col ="purple", lty = 2, lwd = 2)
lines(ts(y_pred - 1.96 * sqrt(P_pred), start= c(2020,2), frequency = 4 ), col ="purple", lty = 2, lwd = 2)
lines(ts(data[293:309, 2],start =c(2020,2),frequency = 4), col = "blue", lwd = 2)
legend("bottomright", col = c("red", "purple","blue"), box.lty = 0, lty = c(1,2,1), legend = c("Forecast", "95% CI", "Actual Data"))
title(main = "Forecast with Kalman filter on LLTM", asdj = 0)
###

arima221 <- arima(ts, order = c(2,2,1))
arima.preds <- predict(arima221, 20)

plot(ts, type = "l",bty = "l", ylab = expression("y"[t]), xlim = c(2000, 2025) ,ylim = range(c(ts[212],arima.preds$pred)))
lines(ts(arima.preds$pred, start= c(2020,2), frequency = 4 ), col ="red", lwd = 2)
lines(ts(arima.preds$pred + 1.96 * arima.preds$se, start= c(2020,2), frequency = 4 ), col ="purple", lty = 2, lwd = 2)
lines(ts(arima.preds$pred - 1.96 * arima.preds$se, start= c(2020,2), frequency = 4 ), col ="purple", lty = 2, lwd = 2)
lines(ts(data[293:309, 2],start =c(2020,2),frequency = 4), col = "blue", lwd = 2)
legend("bottomright", col = c("red", "purple","blue"), box.lty = 0, lty = c(1,2,1), legend = c("Forecast", "95% CI", "Actual Data"))
title(main = "Forecast with ARIMA(2,2,1)", adj = 0)

####
llm.ts <- StructTS(ts,"level")
llm.preds <- predict(llm.ts,n.ahead = 20)

plot(ts, type = "l",bty = "l", ylab = expression("y"[t]), xlim = c(2000, 2025) ,ylim = range(c(ts[212],llm.preds$pred + 1.96 * llm.preds$se)))
lines(ts(llm.preds$pred, start= c(2020,2), frequency = 4 ), col ="red", lwd = 2)
lines(ts(llm.preds$pred + 1.96 * llm.preds$se, start= c(2020,2), frequency = 4 ), col ="purple", lty = 2, lwd = 2)
lines(ts(llm.preds$pred - 1.96 * llm.preds$se, start= c(2020,2), frequency = 4 ), col ="purple", lty = 2, lwd = 2)
lines(ts(data[293:309, 2],start =c(2020,2),frequency = 4), col = "blue", lwd = 2)
legend("bottomright", col = c("red", "purple","blue"), box.lty = 0, lty = c(1,2,1), legend = c("Forecast", "95% CI", "Actual Data"))
title(main = "Forecast with LLM", adj = 0)
