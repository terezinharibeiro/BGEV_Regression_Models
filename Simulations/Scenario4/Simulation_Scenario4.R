source("BGEV_GAMLSS.R") #Function to fit BGEV regression models with xi=0 using gamlss package
if(!suppressWarnings(require(BBmisc))){suppressWarnings(install.packages("BBmisc"));suppressWarnings(library(BBmisc))} #used for is.error function
if(!suppressWarnings(require(VGAM))){suppressWarnings(install.packages("VGAM"));suppressWarnings(library(VGAM))} #used for starting values
if(!suppressWarnings(require(checkmate))){suppressWarnings(install.packages("checkmate"));suppressWarnings(library(checkmate))} 

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

#qBGEV(0.5, mu=200, sigma=1, nu=1.5, tau=1.5)


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


sink("Results_Scenario4.txt")

names_results1 <- c("Estimates_Scenario4_n50", "Estimates_Scenario4_n100", "Estimates_Scenario4_n500")
names_results2 <- c("SE_Estimates_Scenario4_n50", "SE_Estimates_Scenario4_n100", "SE_Estimates_Scenario4_n500")

########################################################################
#######################Global Variables#################################
########################################################################
 
VN <- c(50, 100, 500) #Sample sizes
VBETA <- c(4, 2, 3) #true beta values
VSIGMA <- 1 #true sigma value
VXI <- -0.25 #true xi value
VDELTA <- 0.3 #true delta value

VTHETA <- c(VBETA, VSIGMA, VXI, VDELTA) # true theta values
NREP <- 10000 #number of Monte Carlo replicates
kk1 <- length(VBETA)

#*************generated covariates******************#
se1 = 19 ; se2 = 94 #random seed
suppressWarnings(set.seed(c(se1,se2), kind="Marsaglia-Multicarry")) #To ensure repeatability of the experiment

x1_ini <- runif(VN[1]); 
#x2_ini <- runif(VN[1]);
x2_ini <- rbinom(VN[1], size=1, prob=0.4)

for(l in 1:length(VN))#Loop for the sample size
{
N = VN[l]
if(N == VN[1]){     
		x1 <- x1_ini
		x2 <- x2_ini
		}
		else if (N == 2*VN[1])
		  {
		  x1 <- c(x1_ini, x1_ini)
		  x2 <- c(x2_ini, x2_ini)
		  }
		  else if (N == 10*VN[1])
		  {
		  x1 <- c(x1_ini, x1_ini, x1_ini, x1_ini, x1_ini, x1_ini, x1_ini, x1_ini, x1_ini, x1_ini)
		  x2 <- c(x2_ini, x2_ini, x2_ini, x2_ini, x2_ini, x2_ini, x2_ini, x2_ini, x2_ini, x2_ini)
		  }

X <- model.matrix(~x1+x2)
eta <- as.vector(X%*%VBETA)	
m <- eta
#summary(m)

#***Initializing vectors to store the estimates ***#
thetahat <- matrix(NA, NREP, kk1 + 3)
se_thetahat <- matrix(NA, NREP, kk1 + 3)

cont <- 0 # Counter for Monte Carlo replicates 
f1 <- 0 # Counter for failure in parameter estimation
f2 <- 0 # Counter for failure in hessian estimation
f3 <- 0 # Counter for error failure in parameter estimation

while(cont < NREP)
{
cont <- cont + 1 
#print(cont)
perc <- cont/NREP
#***To print the progress of the simulation***#
if(perc == 0.25 || perc == 0.5 || perc == 0.75 || perc ==1) cat("Perc. Replic. MC =",perc*100,"%","\n")
      
 y <- rBGEV(n=N, mu=m, sigma=VSIGMA, nu=VXI, tau=VDELTA)                                     
 
#---- Starting values via GEV fit (VGAM) ----#

fitgev <-  suppressWarnings(tryCatch(vglm(y ~ X[,2:kk1], gevff(perc = NULL)), error=function(e) {e} )) #starting values via GEV regression
if(is.error(fitgev)){#starts else in case of error in the starting values
muStart <- y
sigmaStart <- rep(1.0,length(y))
nuStart <- rep(0.5, length(y))
tauStart <- rep(0.5, length(y))
}
else
{
estimates <- as.numeric(coef(fitgev)); beta1hat <- estimates[1]; sigmahat <- exp(estimates[2]); xihat <- exp(estimates[3])-0.5;
betascovhat <- estimates[4:(kk1+2)]
betahat <- c(beta1hat, betascovhat)
muStart <- as.vector(X%*%betahat)
sigmaStart <- rep(sigmahat,length(y))
nuStart <- rep(xihat, length(y))
tauStart <- rep(0.5, length(y))
}

 #---- BGEV Fitted model ----#

 fit   <- tryCatch(
gamlss(formula = y~x1+x2,
family = BGEV(mu.link = "identity", sigma.link = "log", nu.link = "own", tau.link="log"),
method = RS(), control = gamlss.control(n.cyc = 500, trace = FALSE),
mu.start = muStart, sigma.start = sigmaStart, nu.start = nuStart), error=function(e) {e})

      
    if(is.error(fit)){# in case of error in the estimation procedure
        f3 <- f3 + 1    
        cont <- cont - 1    
    } 
    else
        { 

          SE_fit <- tryCatch(as.numeric(sqrt(diag(vcov(fit)))), error=function(e) {e})
          
          if(anyNA(SE_fit)||anyNaN(SE_fit)||is.error(SE_fit)){#to avoid hessian problems
            cont <- cont - 1
            f2 <- f2 + 1
          }
          else
             {
              if(fit$converged=="TRUE"){  #estimates only if there is convergence

              #***************parameter estimates*************************#
                         thetahat[cont,] <- c(as.numeric(coef(fit, what = "mu")),
                                              as.numeric(predict(fit, what= "sigma", type = "response")[1]),
                                              as.numeric(predict(fit, what= "nu", type = "response")[1]),
                                              as.numeric(predict(fit, what= "tau", type = "response")[1]))
 
              #***************asymptotic standard error estimates*************************#
                       
                         # --->  For beta's
                         ep.betahat <- sqrt(diag(vcov(fit)[1:kk1,1:kk1]))

                         # --->  For sigma
                         etahat.sigma <- predict(fit, what= "sigma")
                         var.etahat.sigma <-  vcov(fit)[kk1+1,kk1+1]
                         ep.sigmahat <- sqrt(((exp(etahat.sigma[1]))^2)*var.etahat.sigma[1])

                         # ---> For nu
                         etahat.nu <- predict(fit, what= "nu")
                         var.etahat.nu <-  vcov(fit)[kk1+2,kk1+2]
                         ep.nuhat <- sqrt((((1.5 * exp(etahat.nu[1])) / ((1 + exp(etahat.nu[1]))^2))^2)*var.etahat.nu[1])

                         # --->  For tau
                         etahat.tau <- predict(fit, what= "tau")
                         var.etahat.tau <-  vcov(fit)[kk1+3,kk1+3]
                         ep.tauhat <- sqrt(((exp(etahat.tau[1]))^2)*var.etahat.tau[1])

                         se_thetahat[cont,] <- as.numeric(ep.betahat, ep.sigmahat, ep.nuhat, ep.tauhat)
  

             }
             else
                { 
                  cont <- cont - 1
                  f1 <- f1 + 1 
                }
             }
           }
          }#closes MC replicates

#******************results***************************#

#**************mean of estimates*************#
M_thetahat <- apply(thetahat,2,mean)

#************** bias of estimates*************#
B_thetahat <- (M_thetahat - VTHETA)

#**************mean of se estimates*************#
M_se_thetahat <- apply(se_thetahat,2,mean)

#**************root mean squared error of estimates*************#
RMSE_thetahat <- sqrt(apply(apply( t(thetahat) - VTHETA, 1, function(x) x^2), 2, mean) )

###############################################Outputs################################################

out_estimates_MLE <- cbind(M_thetahat, B_thetahat,RMSE_thetahat, M_se_thetahat)
fail_rate1 <- round((f1/NREP)*100,3)
fail_rate2 <- round((f2/NREP)*100,3)
fail_rate3 <- round((f3/NREP)*100,3)
cat(" ========================================================================================== \n")
cat(" TRUE VALUES ","\n")
cat(" Sample size ",N ,"\n")
out_VTHETA <- rbind(VTHETA);colnames(out_VTHETA)=c("beta1", "beta2", "beta3", "sigma", "xi", "delta");rownames(out_VTHETA)=c("");print(out_VTHETA)
Summ_m <- round(c(min(m), max(m), median(m)),3)
out_Sm <- rbind(Summ_m);colnames(out_Sm)=c("m_min", "m_max", "m_median");rownames(out_Sm)=c("");print(out_Sm)
cat(" ========================================================================================== \n")
cat(" ================================ MLE estimates for BGEV regression ================================ \n")
output_mle = rbind(out_estimates_MLE);colnames(output_mle)=c("Mean estimates","Bias","RMSE", "Mean SE");rownames(output_mle)=c("beta1", "beta2", "beta3", "sigma", "xi", "delta");print(output_mle)
cat(" ========================================================================================== \n")
cat("Percentage of estimation failure=", fail_rate1,"%","\n")
cat("Percentage of hessian failure=", fail_rate2,"%","\n")
cat("Percentage of estimation error failure=", fail_rate3,"%","\n")
cat(" ========================================================================================== \n")

#******Print replicates*****#

write.table(cbind(thetahat), 
file= paste(names_results1[l], ".txt", sep = ""), append = FALSE, sep = " ", dec = ".", row.names = FALSE,
col.names = c("MLE_beta1", "MLE_beta2", "MLE_beta3", "MLE_sigma", "MLE_xi", "MLE_delta"))

write.table(cbind(se_thetahat), 
file= paste(names_results2[l], ".txt", sep = ""), append = FALSE, sep = " ", dec = ".", row.names = FALSE,
col.names = c("se_MLE_beta1", "se_MLE_beta2", "se_MLE_beta3", "se_MLE_sigma", "se_MLE_xi", "se_MLE_delta"))

}

sink()


