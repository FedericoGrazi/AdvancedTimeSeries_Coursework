library(TSA)
data("Nile")

y <- Nile

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

llik_optg <- optimg::optimg(c(1000,1000), extract_loglik, 
                            method = "STGD",
                            full = T)

llik_iter <- rbind(c(1000,1000),llik_optg$Estimates)

sel_optg <- c(1,5,10,15,20,26)
q <- llik_iter[,2]/llik_iter[,1]
psi <- log(q)
loglik <- llik_optg$Cost

opt_df <- data.frame(
  sel_optg,
  q = q[sel_optg],
  psi = psi[sel_optg],
  loglik = -llik_optg$Cost[sel_optg]
)

library(plotly)
library(dplyr)
N <- 150
grid.eps  <- seq(1000,20000,length = N)
grid.eta <- seq(1000,20000, length = N)
z <- matrix(0,150,150)
for(i in 1:N){
  for(j in 1:N){
    th <- c(grid.eps[i], grid.eta[j])
    z[i,j] <- -extract_loglik(th)
  }
}
p <- plot_ly(x = grid.eps, y = grid.eta, z = z, type = 'surface')
max_z <- max(z)
max_z_index <- which(z == max_z, arr.ind = TRUE)
max_x <- grid.eps[max_z_index[1]]
max_y <- grid.eta[max_z_index[2]]
line_df <- data.frame(x = c(max_y), y = c(max_x), z = c(-1000,max_z+10))
p <- add_paths(p, data = line_df, x = ~x, y = ~y, z = ~z, color = I("black"))
p <- add_text(p,text = ~paste("sigma_eta:", round(x)[2], "<br>sigma_eps:", round(y)[2], "<br>loglik:", round(z))[2], 
         x = ~x[2], y = ~y[2], z = ~z[2], 
         hoverinfo = "text+x+y+z")
axx <- list(title = "Sigma Eta")
axy <- list(title = "Sigma Epsilon")
axz <- list(title = "Loglikelihood")
p <- layout(p, scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
print(p)