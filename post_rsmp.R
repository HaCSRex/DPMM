library(magrittr)
library(Rcpp)

sourceCpp("post_rsmp.cpp")
set.seed(158501305)
y_mu <- 2
y_sigma <- 1
mu_prior <- 0
sigma_prior <- 1
yn <- rnorm(10, y_mu, y_sigma)

B <- 1000L
N <- 1000L

set.seed(73216911)
mu_postpred <- post_rsmp(yn = yn, y_sigma = y_sigma, mu_prior = mu_prior, sigma_prior = sigma_prior, B = B, N = N)

# posterior function
mu_post <- function(x, x_sigma, mu_prior, sigma_prior){
  n = length(x)
  a = 1/(1/sigma_prior**2 + n/x_sigma**2)
  b = mu_prior/sigma_prior**2 + sum(x)/x_sigma**2
  return(a*b)
}

sigma_post <- function(x, x_sigma, sigma_prior){
  n = length(x)
  a = 1/sigma_prior**2 + n/x_sigma**2
  return(1/a)
}


png("post_rsmp.png", width = 1200, height = 800)

{
  par(mar = c(3, 4, 3, 2), oma = c(0,0,0,0), bg = "#F0F0F0", 
    xpd = FALSE, xaxs = "r", yaxs = "i", mgp = c(2.1,.3,0), 
    las = 1, col.axis = "#434343", col.main = "#343434", 
    cex.lab = 2,
    tck = 0, lend = 1, mfrow = c(1,2))
  
  ydomin <- c(.3, 2.8)#c(min(mu_postpred)*.8, max(mu_postpred)*1.2)
  
  plot(0, 0, 
       type = "n",
       frame.plot = FALSE,
       xlim = c(0,N), ylim = ydomin,
       ylab = expression("Posterior Mean " ~ bar(theta)[n+i]), xlab = "Forward step i")
  
  for(b in 1:B){
    lines(0:N, mu_postpred[b,], col = "#4682B466")
  }
   
  plot(0, 0, 
       type = "n",
       frame.plot = FALSE,
       xlim = c(0,1.5), ylim = ydomin,
       xlab = "Density", ylab = "")
  
  lines(y = seq(ydomin[1], ydomin[2], len = 500), 
        x = dnorm(seq(ydomin[1], ydomin[2], len = 500), 
                  mu_post(yn, y_sigma, mu_prior = mu_prior, sigma_prior = sigma_prior),
                  sqrt(sigma_post(yn, y_sigma, sigma_prior = sigma_prior))), 
        lty = 2, lwd = 2)
  
  density(mu_postpred[, N+1]) %$% lines(y = x, x = y, type = "l", col = "steelblue", lwd = 2)  
  
  legend("top", legend = c(expression("Analytical posterior density " ~ pi(theta ~ "|" ~ y[1:n])), 
                           expression("Kernel density estimate of " ~ bar(theta)[N])), 
         lty = c(2,1), col = c("#000000", "#4682B4"), seg.len = 1, lwd = 2,
         cex = 1.2, text.width = .8)
}

dev.off()

# render a gif
require(animation)
saveGIF({
  for(i in 1:(N/4)){
    par(mar = c(3, 4, 3, 2), oma = c(0,0,0,0), bg = "#F0F0F0", 
        xpd = FALSE, xaxs = "r", yaxs = "i", mgp = c(2.1,.3,0), 
        las = 1, col.axis = "#434343", col.main = "#343434", 
        tck = 0, lend = 1, mfrow = c(1,2))
    
    ydomin <- c(.3, 2.8)
    
    plot(0, 0, 
         type = "n",
         frame.plot = FALSE,
         xlim = c(0,N), ylim = ydomin,
         ylab = expression("Posterior Mean " ~ bar(theta)[n+i]), xlab = "Forward step i")
    
    for(b in 1:B){
      lines(0:(i*4), mu_postpred[b, 1:(i*4+1)], col = "#4682B466")
    }
    
    plot(0, 0, 
         type = "n",
         frame.plot = FALSE,
         xlim = c(0,1.5), ylim = ydomin,
         xlab = "Density", ylab = "")
    
    lines(y = seq(ydomin[1], ydomin[2], len = 500), 
          x = dnorm(seq(ydomin[1], ydomin[2], len = 500), 
                    mu_post(yn, y_sigma, mu_prior = mu_prior, sigma_prior = sigma_prior),
                    sqrt(sigma_post(yn, y_sigma, sigma_prior = sigma_prior))), 
          lty = 2, lwd = 2)
    
    density(mu_postpred[, i*4+1]) %$% lines(y = x, x = y, type = "l", col = "steelblue", lwd = 2)  
  }
}, "post_rsmp.gif", interval = .1, ani.height  = 500, ani.width = 500)



  