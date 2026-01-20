####packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse","data.table","latex2exp", "Cairo", "gamlss", "numDeriv")


#------------------------------------------------------------------------------
# BGEV distribution in GAMLSS with xi=0
#

BGEV0 <- function(mu.link = "identity", sigma.link = "log", nu.link = "log"){
  
  mstats <- checklink("mu.link", "BGEV", substitute(mu.link),
                      c("identity", "log", "inverse", "sqrt", "own"))
  
  dstats <- checklink("sigma.link", "BGEV", substitute(sigma.link),
                      c("log", "inverse", "identity", "sqrt", "own"))
  
  vstats <- checklink("nu.link", "BGEV", substitute(nu.link),
                      c("log", "inverse", "identity", "sqrt", "(0,2]", "logit", "probit", "cloglog",
                      "cauchit", "own"))
  


  structure(
    list(family = c("BGEV0", "Bimodal Generalized Extreme Value distribution (Gumbel)"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE),
         nopar = 3, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         nu.link = as.character(substitute(nu.link)),
         mu.linkfun = mstats$linkfun,
         sigma.linkfun = dstats$linkfun,
         nu.linkfun = vstats$linkfun,
         mu.linkinv = mstats$linkinv,
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         mu.dr = mstats$mu.eta,
         sigma.dr = dstats$mu.eta,
         nu.dr = vstats$mu.eta,
                
         
          ######first order derivates
         dldm =  function(y, mu, sigma, nu){
           dldm = as.vector(attr(gamlss:::numeric.deriv(dBGEV0(y, mu, sigma, nu, log = TRUE), "mu", delta = 1e-04), "gradient"))
           dldm
         },
         
         dldd = function(y, mu, sigma, nu){
           dldd = as.vector(attr(gamlss:::numeric.deriv(dBGEV0(y, mu, sigma, nu, log = TRUE), "sigma", delta = 1e-04), "gradient"))
           dldd
         },
         
         dldv = function(y, mu, sigma, nu){
           dldv = as.vector(attr(gamlss:::numeric.deriv(dBGEV0(y, mu, sigma, nu, log = TRUE), "nu", delta = 1e-04), "gradient"))
           dldv
         },

         #second derivates

          d2ldm2 = function(y, mu, sigma, nu) {
          dldm = as.vector(attr(gamlss:::numeric.deriv(dBGEV0(y, mu, sigma, nu, log = TRUE), "mu", delta = 1e-04), "gradient"))
          d2ldm2 = -dldm^2
          d2ldm2
          },

           d2ldd2 = function(y, mu, sigma, nu) {
           dldd = as.vector(attr(gamlss:::numeric.deriv(dBGEV0(y, mu, sigma, nu, log = TRUE), "sigma", delta = 1e-04), "gradient"))
           d2ldd2 = -dldd^2
           d2ldd2
           },
 
           d2ldv2 = function(y, mu, sigma, nu) {
            dldv = as.vector(attr(gamlss:::numeric.deriv(dBGEV0(y, mu, sigma, nu, log = TRUE), "nu", delta = 1e-04), "gradient"))
            d2ldv2 = -dldv^2
           d2ldv2
           }, 
  
           d2ldmdd = function(y, mu, sigma, nu) {
           dldm = as.vector(attr(gamlss:::numeric.deriv(dBGEV0(y, mu, sigma, nu, log = TRUE), "mu", delta = 1e-04), "gradient"))
           dldd = as.vector(attr(gamlss:::numeric.deriv(dBGEV0(y, mu, sigma, nu, log = TRUE), "sigma", delta = 1e-04), "gradient"))
           
           d2ldmdd = -(dldm*dldd)
           #d2ldmdd = ifelse(d2ldmdd < -1e-15, d2ldmdd,-1e-15)#add
           d2ldmdd
         },
         
         d2ldmdv = function(y, mu, sigma, nu) {
           dldm = as.vector(attr(gamlss:::numeric.deriv(dBGEV0(y, mu, sigma, nu, log = TRUE), "mu", delta = 1e-04), "gradient"))
           dldv = as.vector(attr(gamlss:::numeric.deriv(dBGEV0(y, mu, sigma, nu, log = TRUE), "nu", delta = 1e-04), "gradient"))
           d2ldmdv = -(dldm*dldv)
           #d2ldmdv = ifelse(d2ldmdv < -1e-15, d2ldmdv,-1e-15)#add
           d2ldmdv
         },

         
         d2ldddv = function(y, mu, sigma, nu) {
           dldd = as.vector(attr(gamlss:::numeric.deriv(dBGEV0(y, mu, sigma, nu, log = TRUE), "sigma", delta = 1e-04), "gradient"))
           dldv = as.vector(attr(gamlss:::numeric.deriv(dBGEV0(y, mu, sigma, nu, log = TRUE), "nu", delta = 1e-04), "gradient"))
           d2ldddv = -(dldd*dldv)
           #d2ldddv = ifelse(d2ldddv < -1e-15, d2ldddv,-1e-15)#add
           d2ldddv
         },
         
                 
         #####
         G.dev.incr = function(y, mu, sigma, nu, ...) -2*dBGEV0(y = y, mu = mu, sigma = sigma, nu=nu, log=TRUE),
         rqres = expression(rqres(pfun="pBGEV0", type = "Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
         
         #####Initial values
         
            
         mu.initial = expression(mu <-  y),
         sigma.initial = expression(sigma <- rep(1.0,length(y))),
         nu.initial = expression(nu <- rep(0.5, length(y))),      
                  #####Restrictions
         mu.valid = function(mu) TRUE,
         sigma.valid = function(sigma) all(is.finite(sigma) & sigma > 0),
         nu.valid = function(nu) all(is.finite(nu) & nu > 0),
         y.valid = function(y) TRUE
    ),
    class = c("gamlss.family","family"))
}
# Density function

dBGEV0 <- function(y, mu = 1, sigma = 1, nu = 0.5, log = FALSE) {
  if (any(sigma <= 0)) warning("sigma must be positive")
  if (any(nu <= 0)) warning("nu must be > 0")
   
  eps <- .Machine$double.eps  # ~2.2e-16
  
 # Compute auxiliary variables:

  mub <- mu - (-sigma*log(-log(0.5)))^(1/(nu+1))
  y_mub <- y - mub

  abs_y_mub <- pmax(abs(y_mub), eps)  
  T <- y_mub * abs_y_mub^nu
  T_div_sigma <- pmin(pmax(T / sigma, -700), 700) 
  

  # Derivada de T
  derivate_T <- (nu + 1) * abs_y_mub^nu
  derivate_T_safe <- pmin(pmax(derivate_T, 1e-10), 1e10)

  fy1 <- (1/sigma)*exp(-T_div_sigma)*exp(-exp(-T_div_sigma))*derivate_T_safe 

    if(log == TRUE) fy <- log(fy1) else fy <- fy1
    return(fy)

}


#-------------------------------------------------------------------------------
# Cumulative function
pBGEV0 <- function(q, mu=1, sigma=1, nu=0.5, lower.tail = TRUE, log = FALSE){
  if (any(sigma <= 0)) warning("sigma must be positive")
  if (any(nu <= 0)) warning("nu must be greater or equal to 0")
    
  # Compute auxiliary variables:
  mub <- mu - (-sigma*log(-log(0.5)))^(1/(nu+1))
  T      <- (q-mub)*(abs(q-mub)^nu)
       
  cdf <- exp(-exp(-T/sigma))

  if(lower.tail==TRUE) cdf <- cdf else cdf <- 1 - cdf
  if(log==TRUE) cdf <- log(cdf)
  
  return(cdf)
}


#-------------------------------------------------------------------------------
# Quantile function
qBGEV0 <- function(p, mu=1, sigma=1, nu=0.5, lower.tail = TRUE, log = FALSE){
  if (any(sigma <= 0)) warning("sigma must be positive")
  if (any(nu <= 0)) warning("nu must be greater than 0")

  if (log==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  
# Compute:
   mub <- mu - (-sigma*log(-log(0.5)))^(1/(nu+1))

   quantile <- mub + sign(-sigma*log(-log(p)))*((abs(-sigma*log(-log(p))))^(1/(nu+1))) 
  
  return(quantile)
}

#-------------------------------------------------------------------------------
# Pseudo-random function
rBGEV0 <- function(n, mu=1, sigma=1, nu=0.5){
  if (any(sigma <= 0)) warning("sigma must be positive")
  if (any(nu <= 0)) warning("nu must be greater than 0")
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  
  U <- runif(n)
  rnumber <- qBGEV0(U, mu, sigma, nu)
   
  return(rnumber)
}
