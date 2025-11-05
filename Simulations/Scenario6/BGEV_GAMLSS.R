####packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse","data.table","latex2exp", "Cairo", "gamlss", "numDeriv")
  
  own.linkfun <- function(mu) {
  if (any(mu < -0.5 | mu > 1, na.rm=T)) {
    warning("nu must be in [-0.5, 1]")
    mu <- pmin(pmax(mu, -0.499), 0.999)  
  }
  log((mu + 0.5) / (1 - mu))
  }
    
  own.linkinv <- function(eta) {
     
    (exp(eta)-0.5)/(1 + exp(eta))
  }
  
  own.mu.eta <- function(eta) {
    # Derivate of linkinv
    (1.5 * exp(eta)) / ((1 + exp(eta))^2)
  }
  
  own.valideta <- function(eta) TRUE  # eta in R
  

#------------------------------------------------------------------------------
# BGEV distribution in GAMLSS
#

BGEV <- function(mu.link = "identity", sigma.link = "log", nu.link = "identity", tau.link = "logit"){
  
  mstats <- checklink("mu.link", "BGEV", substitute(mu.link),
                      c("identity", "log", "inverse", "sqrt", "own"))
  
  dstats <- checklink("sigma.link", "BGEV", substitute(sigma.link),
                      c("log", "inverse", "identity", "sqrt", "own"))
  
 # vstats <- checklink("nu.link", "BGEV", substitute(nu.link),
  #                    c("identity", "log", "inverse", "sqrt", "[-1,1]", "own"))
  
  vstats <- if (inherits(nu.link, "link.gamlss")) {
  nu.link
    } else {
       checklink("nu.link", "BGEV", substitute(nu.link),
            c("identity", "log", "inverse", "sqrt", "[-1,1]", "own"))
   }

  tstats <- checklink("tau.link", "BGEV", substitute(tau.link),
                      c("log", "inverse", "identity", "sqrt", "(0,2]", "logit", "probit", "cloglog",
                      "cauchit", "own"))


  structure(
    list(family = c("BGEV", "Bimodal Generalized Extreme Value distribution"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE),
         nopar = 4, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         nu.link = as.character(substitute(nu.link)),
         tau.link = as.character(substitute(tau.link)),
         mu.linkfun = mstats$linkfun,
         sigma.linkfun = dstats$linkfun,
         nu.linkfun = vstats$linkfun,
         tau.linkfun = tstats$linkfun,
         mu.linkinv = mstats$linkinv,
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         tau.linkinv = tstats$linkinv,
         mu.dr = mstats$mu.eta,
         sigma.dr = dstats$mu.eta,
         nu.dr = vstats$mu.eta,
         tau.dr = tstats$mu.eta,
        
         
         ######first order derivates
         dldm =  function(y, mu, sigma, nu, tau){
           dldm = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "mu", delta = 1e-04), "gradient"))
           dldm
         },
         
         dldd = function(y, mu, sigma, nu, tau){
           dldd = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "sigma", delta = 1e-04), "gradient"))
           dldd
         },
         
         dldv = function(y, mu, sigma, nu, tau){
           dldv = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "nu", delta = 1e-04), "gradient"))
           dldv
         },

           dldt = function(y, mu, sigma, nu, tau){
           dldt = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "tau", delta = 1e-04), "gradient"))
           dldt
         },

          #second derivates

          d2ldm2 = function(y, mu, sigma, nu, tau) {
          dldm = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "mu", delta = 1e-04), "gradient"))
          d2ldm2 = -dldm^2
          d2ldm2
          },

           d2ldd2 = function(y, mu, sigma, nu, tau) {
           dldd = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "sigma", delta = 1e-04), "gradient"))
           d2ldd2 = -dldd^2
           d2ldd2
           },
 
           d2ldv2 = function(y, mu, sigma, nu, tau) {
            dldv = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "nu", delta = 1e-04), "gradient"))
            d2ldv2 = -dldv^2
           d2ldv2
           }, 

         d2ldt2 = function(y, mu, sigma, nu, tau) {
           dldt = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "tau", delta = 1e-04), "gradient"))
           d2ldt2 = -dldt^2
           d2ldt2  
         }, 
         
         
           d2ldmdd = function(y, mu, sigma, nu, tau) {
           dldm = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "mu", delta = 1e-04), "gradient"))
           dldd = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "sigma", delta = 1e-04), "gradient"))
           
           d2ldmdd = -(dldm*dldd)
           d2ldmdd
         },
         
         d2ldmdv = function(y, mu, sigma, nu, tau) {
           dldm = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "mu", delta = 1e-04), "gradient"))
           dldv = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "nu", delta = 1e-04), "gradient"))
           d2ldmdv = -(dldm*dldv)
           #d2ldmdv = ifelse(d2ldmdv < -1e-15, d2ldmdv,-1e-15)#add
           d2ldmdv
         },

         d2ldmdt = function(y, mu, sigma, nu, tau) {
           dldm = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "mu", delta = 1e-04), "gradient"))
           dldt = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "tau", delta = 1e-04), "gradient"))
           d2ldmdt = -(dldm*dldt)
           #d2ldmdt = ifelse(d2ldmdt < -1e-15, d2ldmdt,-1e-15)#add
           d2ldmdt
         },

         d2ldddv = function(y, mu, sigma, nu, tau) {
           dldd = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "sigma", delta = 1e-04), "gradient"))
           dldv = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "nu", delta = 1e-04), "gradient"))
           d2ldddv = -(dldd*dldv)
           #d2ldddv = ifelse(d2ldddv < -1e-15, d2ldddv,-1e-15)#add
           d2ldddv
         },
         
         d2ldddt = function(y, mu, sigma, nu, tau) {
           dldd = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "sigma", delta = 1e-04), "gradient"))
           dldt = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "tau", delta = 1e-04), "gradient"))
           d2ldddt = -(dldd*dldt)
           #d2ldddt = ifelse(d2ldddt < -1e-15, d2ldddt,-1e-15)#add
           d2ldddt
         }, 
          
         d2ldvdt = function(y, mu, sigma, nu, tau) {
           dldv = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "nu", delta = 1e-04), "gradient"))
           dldt = as.vector(attr(gamlss:::numeric.deriv(dBGEV(y, mu, sigma, nu, tau, log = TRUE), "tau", delta = 1e-04), "gradient"))
           d2ldvdt = -(dldv*dldt)
           #d2ldvdt = ifelse(d2ldvdt < -1e-15, d2ldvdt,-1e-15)#add
           d2ldvdt
         }, 
      

         
         #####
         G.dev.incr = function(y, mu, sigma, nu, tau, ...){
          -2*dBGEV(y = y, mu = mu, sigma = sigma, nu=nu, tau=tau, log=TRUE)
         },
         rqres = expression(rqres(pfun="pBGEV", type = "Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)),
         
         #####Initial values
         
            
         mu.initial = expression(mu <-  y),
         sigma.initial = expression(sigma <- rep(1.0,length(y))),
         nu.initial = expression(nu <- rep(0.5, length(y))),      
         tau.initial = expression(tau <- rep(0.5, length(y))),  

         #####Restrictions
         mu.valid = function(mu) TRUE,
         sigma.valid = function(sigma) all(is.finite(sigma) & sigma > 0),
         nu.valid = function(nu) {
           valid <- all(is.finite(nu) & nu >= -0.5 & nu <=1)
           if (!valid) warning("Algum valor de nu está fora do intervalo [-0.5, 1]")
           TRUE  # Retorna TRUE para não interromper
         },
         tau.valid = function(tau) all(is.finite(tau) & sigma >= 0),
         y.valid = function(y) TRUE
    ),

    class = c("gamlss.family","family"))
}

# Density function

dBGEV <- function(y, mu = 1, sigma = 1, nu = 0.5, tau = 0.5, log = FALSE) {
  if (any(sigma <= 0)) warning("sigma must be positive")
  if (any(tau < 0)) warning("tau must be >= 0")  
  if (any(nu < -0.5 | nu > 1)) warning("If link.nu=own, nu must be >= -0.5 and <= 1")

  eps <- .Machine$double.eps  # ~2.2e-16

  log05 <- log(0.5)  # = -log(2), estável
  temp <- (sigma / nu) * ((-log05)^(-nu) - 1)
  temp_safe <- sign(temp) * (abs(temp) + eps)^(1 / (tau + 1))  
  mub <- mu - temp_safe
  y_mub <- y - mub
  abs_y_mub <- pmax(abs(y_mub), eps) 
  T <- y_mub * abs_y_mub^tau
  derivate_T <- (tau + 1) * abs_y_mub^tau
  derivate_T_safe <- pmin(pmax(derivate_T, 1e-10), 1e10)
  inner <- 1 + nu * (T / sigma)
  inner_safe <- pmax(pmin(inner, 1e10), 1e-10) 
  exponent <- - (1 / nu) * log(inner_safe)
  exponent <- pmin(pmax(exponent, -700), 700)
  inner_pow <- exp(exponent)  

  log_fy <- -log(sigma) +
            (-1 / nu - 1) * log(inner_safe) -
            inner_pow +
            log(derivate_T_safe)
  if (log) return(log_fy) else return(exp(log_fy))
}

#-------------------------------------------------------------------------------
# Cumulative function
pBGEV <- function(q, mu=1, sigma=1, nu=0.5, tau=0.5, lower.tail = TRUE, log = FALSE){
  if (any(sigma <= 0)) warning("sigma must be positive")
  if (any(nu < -0.5 | nu > 1)) warning("If link.nu=own, nu must be >= -0.5 and <= 1")
  if (any(tau < 0)) warning("tau must be greater or equal to 0")
    
  # Compute auxiliary variables:
  mub <-  mu - sign((sigma/nu)*((-log(0.5))^(-nu) -1))*(abs((sigma/nu)*((-log(0.5))^(-nu) -1))^(1/(tau+1)))
  T      <- (q-mub)*(abs(q-mub)^tau)
       
  cdf <- exp(-((1+nu*(T/sigma))^(-1/nu)))

  if(lower.tail==TRUE) cdf <- cdf else cdf <- 1 - cdf
  if(log==TRUE) cdf <- log(cdf)
  
  return(cdf)
}


#-------------------------------------------------------------------------------
# Quantile function
qBGEV <- function(p, mu=1, sigma=1, nu=0.5, tau=0.5, lower.tail = TRUE, log = FALSE){
  if (any(sigma <= 0)) warning("sigma must be positive")
  if (any(nu < -0.5 | nu > 1)) warning("If link.nu=own, nu must be >= -0.5 and <= 1")
  if (any(tau < 0)) warning("tau must be greater than 0")

  if (log==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  
# Compute:
  mub <-  mu - sign((sigma/nu)*((-log(0.5))^(-nu) -1))*(abs((sigma/nu)*((-log(0.5))^(-nu) -1))^(1/(tau+1)))
 
  quantile <- mub + sign((sigma/nu)*(((-log(p))^(-nu))-1))*((abs((sigma/nu)*(((-log(p))^(-nu))-1)))^(1/(tau + 1))) 
  
  return(quantile)
}

#-------------------------------------------------------------------------------
# Pseudo-random function
rBGEV <- function(n, mu=1, sigma=1, nu=0.5, tau=0.5){
  if (any(sigma <= 0)) warning("sigma must be positive")
  if (any(nu < -0.5 | nu > 1)) warning("If link.nu=own, nu must be >= -0.5 and <= 1")
  if (any(tau < 0)) warning("tau must be greater than 0")

  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  
  U <- runif(n)
  mub <-  mu - sign((sigma/nu)*((-log(0.5))^(-nu) -1))*(abs((sigma/nu)*((-log(0.5))^(-nu) -1))^(1/(tau+1)))
  rnumber <- qBGEV(U, mu, sigma, nu, tau)
   
  return(rnumber)
}

