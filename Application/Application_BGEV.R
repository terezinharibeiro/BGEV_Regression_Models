
source("BGEV_GAMLSS.R") #Function to fit BGEV regression models using gamlss package, where mu denotes the 0.50 quantile (median) and xi\=0 
source("BGEV0_GAMLSS.R") #Function to fit BGEV regression models using gamlss package, where mu denotes the 0.50 quantile (median) and xi=0
source("BGEV_GAMLSS_q99.R") #Function to fit BGEV regression models using gamlss package, where mu denotes the 0.99 quantile (median) and xi\=0
source("BGEV0_GAMLSS_q99.R") #Function to fit BGEV regression models using gamlss package, where mu denotes the 0.99 quantile (median) and xi=0

library(ggplot2)
#instal.packages("qqboxplot")
library(qqboxplot) 
library(dplyr)

#**** Read data ****#
DTP_data <- read_rds(file = 'DTPdata.rds')
attach(DTP_data)

y <- Temperatura_orvalho # minimum dew point temperature - response variable
HUM <- Umidade_rel - mean(Umidade_rel) #covariate minimum humidity (HUM) centered
WS <- Vento_rajada_max - mean(Vento_rajada_max) #covariate maximum wind speed (WS) centered
S <- estacoes #covariate season (S) -  S = 1 corresponds to the rainy season (October to April),
# and  S = 0 corresponds to the dry season (May to September)
P <- Pressao_media - mean(Pressao_media) #covariate maximum pressure (P) centered


#------- Scatter plots ------- #

# Response versus HUM
dados <- data.frame(Umidade_rel, y)
ggplot(dados, aes(x=Umidade_rel, y=y)) +
  geom_point(size=2, shape=19)+
  xlab("")+ylab("")+
  theme_bw()+
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"))
ggsave("DTPHUM.png", width = 800, height = 500, units = 'px')

# Response versus WS
dados <- data.frame(Vento_rajada_max, y)
ggplot(dados, aes(x=Vento_rajada_max, y=y)) +
  geom_point(size=2, shape=19)+
  xlab("")+ylab("")+
  theme_bw()+
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"))
ggsave("DTPWS.png", width = 800, height = 500, units = 'px')

# Response versus P
dados <- data.frame(Pressao_media, y)
ggplot(dados, aes(x=Pressao_media, y=y)) +
  geom_point(size=2, shape=19)+
  xlab("")+ylab("")+
  theme_bw()+
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"))
ggsave("DTPP.png", width = 800, height = 500, units = 'px')

# Response versus S
dados <- data.frame(S, y)
dados$S <- factor(dados$S, levels = c(0, 1),labels = c("Dry", "Rainy"))
ggplot(dados, aes(x=as.factor(S), y=y)) + 
    geom_boxplot(fill="slateblue", alpha=0.2) +
    xlab("")+ylab("")+
  theme_bw()+
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"))
ggsave("boxplot.png", width = 800, height = 500, units = 'px' )


#*** Descritive to assess potential multicollinearity ***#

DTP_VC <- cbind(y, HUM, WS, P)
cor(DTP_VC) #Pearson correlations

#boxplots of quantitative covariates versus Season
boxplot(HUM~S)
boxplot(WS~S)
boxplot(P~S)

#***** BGEV fitted regression models with xi=0 ****#

#install.packages(VGAM)
require(VGAM)


#**** Specification 1 : ~HUM ****#

#---- Starting values via Gumbel fit (VGAM) ----#

X <- model.matrix(~HUM)
kk1 <- ncol(X);  n <- nrow(X)
fitgumbel <-  suppressWarnings(vglm(y ~ X[,2:kk1], gumbelff(perc = NULL))) #starting values via Gumbel regression
estimates <- as.numeric(coef(fitgumbel)); beta1hat <- estimates[1]; sigmahat <- exp(estimates[2]);
betascovhat <- estimates[3:(kk1+1)]
betahat <- c(beta1hat, betascovhat)
muStart <- as.vector(X%*%betahat)
sigmaStart <- rep(sigmahat,length(y))
nuStart <- rep(0.5, length(y))

AIC_GUM_SP1 <- -2*logLik(fitgumbel) + 2*(kk1+1);  AIC_GUM_SP1


#---- BGEV Fitted model ----#

# *** Median fit - 0.50 quantile fit regression ***#

fit_BGEV0_SP1 <- gamlss(formula = y~HUM,
family = BGEV0(mu.link = "identity", sigma.link = "identity", nu.link = "identity"),
method = CG(), control = gamlss.control(n.cyc = 200, trace = FALSE),
mu.start = muStart, sigma.start = sigmaStart, nu.start = nuStart)
summary(fit_BGEV0_SP1)

#*** Diagnostics ***#
par(mar=c(5.0,5.0,4.0,2.0))
wp(fit_BGEV0_SP1, cex = 1.5, cex.lab = 2.0) #worm plot
plot(fit_BGEV0_SP1) #residuals plots

#--- qqboxplot ---#

rq_BGEV0_SP1 <- resid(fit_BGEV0_SP1) #quantile residuals

data <-  data.frame(y = rq_BGEV0_SP1)
  ggplot(data, aes(y = rq_BGEV0_SP1)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=25),
        axis.title.x = element_text(size=25),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid = element_line(colour = "grey70"))


# *** 0.99 quantile fit regression ***#

muStart <- y
#sigmaStart and nuStart are the same values used previously

fit_BGEV0_q99_SP1 <- gamlss(formula = y~HUM,
family = BGEV0_q99(mu.link = "identity", sigma.link = "identity", nu.link = "identity"),
method = CG(), control = gamlss.control(n.cyc = 200, trace = FALSE),
mu.start = muStart, sigma.start = sigmaStart, nu.start = nuStart)
summary(fit_BGEV0_q99_SP1)


#**** Specification 2: ~HUM+S ****#

#---- Starting values via Gumbel fit (VGAM) ----#

X <- model.matrix(~HUM+S)
kk1 <- ncol(X);  n <- nrow(X)
fitgumbel <-  suppressWarnings(vglm(y ~ X[,2:kk1], gumbelff(perc = NULL))) #starting values via Gumbel regression
estimates <- as.numeric(coef(fitgumbel)); beta1hat <- estimates[1]; sigmahat <- exp(estimates[2]);
betascovhat <- estimates[3:(kk1+1)]
betahat <- c(beta1hat, betascovhat)
muStart <- as.vector(X%*%betahat)
sigmaStart <- rep(sigmahat,length(y))
nuStart <- rep(0.5, length(y))

AIC_GUM_SP2 <- -2*logLik(fitgumbel) + 2*(kk1+1);  AIC_GUM_SP2

#---- BGEV Fitted model ----#

# *** Median fit - 0.50 quantile fit regression ***#

fit_BGEV0_SP2 <- gamlss(formula = y~HUM+S,  
family = BGEV0(mu.link = "identity", sigma.link = "identity", nu.link = "identity"),
method = CG(), control = gamlss.control(n.cyc = 200, trace = FALSE),
mu.start = muStart, sigma.start = sigmaStart, nu.start = nuStart)
summary(fit_BGEV0_SP2)

#*** Diagnostics ***#

par(mar=c(5.0,5.0,4.0,2.0))
wp(fit_BGEV0_SP2, cex = 1.5, cex.lab = 2.0)
plot(fit_BGEV0_SP2)

#--- qqboxplot ---#

rq_BGEV0_SP2 <- resid(fit_BGEV0_SP2) #quantile residuals

data <-  data.frame(y = rq_BGEV0_SP2)
  ggplot(data, aes(y = rq_BGEV0_SP2)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=25),
        axis.title.x = element_text(size=25),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid = element_line(colour = "grey70"))


# *** 0.99 quantile fit regression ***#

muStart <- y
#sigmaStart and nuStart are the same values used previously

fit_BGEV0_q99_SP2 <- gamlss(formula = y~HUM+S,  
family = BGEV0_q99(mu.link = "identity", sigma.link = "identity", nu.link = "identity"),
method = CG(), control = gamlss.control(n.cyc = 200, trace = FALSE),
mu.start = muStart, sigma.start = sigmaStart, nu.start = nuStart)
summary(fit_BGEV0_q99_SP2)


#**** Specification 3: ~HUM+S+P ****#

#---- Starting values via Gumbel fit (VGAM) ----#

X <- model.matrix(~HUM+S+P)
kk1 <- ncol(X);  n <- nrow(X)
fitgumbel <-  suppressWarnings(vglm(y ~ X[,2:kk1], gumbelff(perc = NULL))) #starting values via Gumbel regression
estimates <- as.numeric(coef(fitgumbel)); beta1hat <- estimates[1]; sigmahat <- exp(estimates[2]);
betascovhat <- estimates[3:(kk1+1)]
betahat <- c(beta1hat, betascovhat)
muStart <- as.vector(X%*%betahat)
sigmaStart <- rep(sigmahat,length(y))
nuStart <- rep(0.5, length(y))

AIC_GUM_SP3 <- -2*logLik(fitgumbel) + 2*(kk1+1);  AIC_GUM_SP3

#---- BGEV Fitted model ----#

# *** Median fit - 0.50 quantile fit regression ***#

fit_BGEV0_SP3 <- gamlss(formula = y~HUM+S+P,  
family = BGEV0(mu.link = "identity", sigma.link = "identity", nu.link = "identity"),
method = CG(), control = gamlss.control(n.cyc = 200, trace = FALSE),
mu.start = muStart, sigma.start = sigmaStart, nu.start = nuStart)
summary(fit_BGEV0_SP3)

#*** Diagnostics ***#

par(mar=c(5.0,5.0,4.0,2.0))
wp(fit_BGEV0_SP3, cex = 1.5, cex.lab = 2.0)
plot(fit_BGEV0_SP3)

#--- qqboxplot ---#

rq_BGEV0_SP3 <- resid(fit_BGEV0_SP3) #quantile residuals

data <-  data.frame(y = rq_BGEV0_SP3)
  ggplot(data, aes(y = rq_BGEV0_SP3)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=25),
        axis.title.x = element_text(size=25),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid = element_line(colour = "grey70"))


# *** 0.99 quantile fit regression ***#

muStart <- y
#sigmaStart and nuStart are the same values used previously

fit_BGEV0_q99_SP3 <- gamlss(formula = y~HUM+S+P,  
family = BGEV0_q99(mu.link = "identity", sigma.link = "identity", nu.link = "identity"),
method = CG(), control = gamlss.control(n.cyc = 200, trace = FALSE),
mu.start = muStart, sigma.start = sigmaStart, nu.start = nuStart)
summary(fit_BGEV0_q99_SP3)


#**** Specification 4: ~HUM+S+P+WS ****#

#---- Starting values via Gumbel fit (VGAM) ----#

X <- model.matrix(~HUM+S+P+WS)
kk1 <- ncol(X);  n <- nrow(X)
fitgumbel <-  suppressWarnings(vglm(y ~ X[,2:kk1], gumbelff(perc = NULL))) #starting values via Gumbel regression
estimates <- as.numeric(coef(fitgumbel)); beta1hat <- estimates[1]; sigmahat <- exp(estimates[2]);
betascovhat <- estimates[3:(kk1+1)]
betahat <- c(beta1hat, betascovhat)
muStart <- as.vector(X%*%betahat)
sigmaStart <- rep(sigmahat,length(y))
nuStart <- rep(0.5, length(y))

AIC_GUM_SP4 <- -2*logLik(fitgumbel) + 2*(kk1+1);  AIC_GUM_SP4


#---- BGEV Fitted model ----#

# *** Median fit - 0.50 quantile fit regression ***#

fit_BGEV0_SP4 <- gamlss(formula = y~HUM+S+P+WS,
family = BGEV0(mu.link = "identity", sigma.link = "identity", nu.link = "identity"),
method = CG(), control = gamlss.control(n.cyc = 200, trace = FALSE),
mu.start = muStart, sigma.start = sigmaStart, nu.start = nuStart)
summary(fit_BGEV0_SP4)

#*** Diagnostics ***#

par(mar=c(5.0,5.0,4.0,2.0))
wp(fit_BGEV0_SP4, cex = 1.5, cex.lab = 2.0)
plot(fit_BGEV0_SP4)

#--- qqboxplot ---#

rq_BGEV0_SP4 <- resid(fit_BGEV0_SP4) #quantile residuals

data <-  data.frame(y = rq_BGEV0_SP4)
  ggplot(data, aes(y = rq_BGEV0_SP4)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=25),
        axis.title.x = element_text(size=25),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid = element_line(colour = "grey70"))


# *** 0.99 quantile fit regression ***#

muStart <- y
#sigmaStart and nuStart are the same values used previously

fit_BGEV0_q99_SP4 <- gamlss(formula = y~HUM+S+P+WS,  
family = BGEV0_q99(mu.link = "identity", sigma.link = "identity", nu.link = "identity"),
method = CG(), control = gamlss.control(n.cyc = 200, trace = FALSE),
mu.start = muStart, sigma.start = sigmaStart, nu.start = nuStart)
summary(fit_BGEV0_q99_SP4)



#***** BGEV fitted regression models with xi\=0 ****#


#**** Specification 1 : ~HUM ****#

#---- GEV Fitted model ----#

X <- model.matrix(~HUM)
kk1 <- ncol(X);  n <- nrow(X)
fitgev <-  suppressWarnings(tryCatch(vglm(y ~ X[,2:kk1], gevff(perc = NULL)), error=function(e) {e} )) #starting values via GEV regression
AIC_GEV_SP1 <- -2*logLik(fitgev) + 2*(kk1+2);  AIC_GEV_SP1

#---- BGEV Fitted model ----#

# *** Median fit - 0.50 quantile fit regression ***#

fit_BGEV_SP1 <- gamlss(formula = y~HUM,  
family = BGEV(mu.link = "identity", sigma.link = "log", nu.link = "identity", tau.link = "log"),
method = CG(), control = gamlss.control(n.cyc = 200, trace = FALSE))
summary(fit_BGEV_SP1)

#*** Diagnostics ***#
par(mar=c(5.0,5.0,4.0,2.0))
wp(fit_BGEV_SP1, cex = 1.5, cex.lab = 2.0)
plot(fit_BGEV_SP1)

#--- qqboxplot ---#

rq_BGEV_SP1 <- resid(fit_BGEV_SP1) #quantile residuals

data <-  data.frame(y = rq_BGEV_SP1)
  ggplot(data, aes(y = rq_BGEV_SP1)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=25),
        axis.title.x = element_text(size=25),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid = element_line(colour = "grey70"))


#**** Parameter estimates ****#

sigmahat <- predict(fit_BGEV_SP1, what= "sigma", type = "response"); sigmahat
tauhat <- predict(fit_BGEV_SP1, what= "tau", type = "response"); tauhat

#**** Parameter standard error estimates via Delta method ****#

# --->  For sigma

#Derivate of h(theta) = exp(theta) (inverse of log function) is h'(theta) = exp(theta)

etahat <- predict(fit_BGEV_SP1, what= "sigma")
var.etahat <-  vcov(fit_BGEV_SP1)[kk1+1,kk1+1]
ep.sigmahat <- sqrt(((exp(etahat))^2)*var.etahat); ep.sigmahat

# --->  For tau

etahat <- predict(fit_BGEV_SP1, what= "tau")
var.etahat <-  vcov(fit_BGEV_SP1)[kk1+3,kk1+3]
ep.tauhat <- sqrt(((exp(etahat))^2)*var.etahat); ep.tauhat



# *** 0.99 quantile fit regression ***#

fit_BGEV_q99_SP1 <- gamlss(formula = y~HUM,  
family = BGEV_q99(mu.link = "identity", sigma.link = "log", nu.link = "identity", tau.link = "log"),
method = CG(), sigma.start = 1.2,control = gamlss.control(n.cyc = 200, trace = FALSE))
summary(fit_BGEV_q99_SP1)

#**** Parameter estimates ****#

sigmahat <- predict(fit_BGEV_q99_SP1, what= "sigma", type = "response"); sigmahat
tauhat <- predict(fit_BGEV_q99_SP1, what= "tau", type = "response"); tauhat

#**** Parameter standard error estimates via Delta method ****#

# --->  For sigma

#Derivate of h(theta) = exp(theta) (inverse of log function) is h'(theta) = exp(theta)

obj <- summary(fit_BGEV_q99_SP1, type= "qr")

etahat <- predict(fit_BGEV_q99_SP1, what= "sigma")
ep.sigmahat <- sqrt((exp(etahat))^2)*as.numeric(obj[,2][kk1+1]); ep.sigmahat

# --->  For tau

etahat <- predict(fit_BGEV_q99_SP1, what= "tau")
ep.tauhat <- sqrt((exp(etahat))^2)*as.numeric(obj[,2][kk1+3]); ep.tauhat



#**** Specification 2 : ~HUM+S ****#

#---- GEV Fitted model ----#

X <- model.matrix(~HUM+S)
kk1 <- ncol(X);  n <- nrow(X)
fitgev <-  suppressWarnings(tryCatch(vglm(y ~ X[,2:kk1], gevff(perc = NULL)), error=function(e) {e} )) #starting values via GEV regression
AIC_GEV_SP2 <- -2*logLik(fitgev) + 2*(kk1+2);  AIC_GEV_SP2


#---- BGEV Fitted model ----#

# *** Median fit - 0.50 quantile fit regression ***#

fit_BGEV_SP2 <- gamlss(formula = y~HUM+S, 
family = BGEV(mu.link = "identity", sigma.link = "log", nu.link = "own", tau.link = "log"),
method = CG(), control = gamlss.control(n.cyc = 200, trace = FALSE))
summary(fit_BGEV_SP2)

#*** Diagnostics ***#

par(mar=c(5.0,5.0,4.0,2.0))
wp(fit_BGEV_SP2, cex = 1.5, cex.lab = 2.0)
plot(fit_BGEV_SP2)

#--- qqboxplot ---#

rq_BGEV_SP2 <- resid(fit_BGEV_SP2) #quantile residuals

data <-  data.frame(y = rq_BGEV_SP2)
  ggplot(data, aes(y = rq_BGEV_SP2)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=25),
        axis.title.x = element_text(size=25),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid = element_line(colour = "grey70"))


#**** Parameter estimates ****#

sigmahat <- predict(fit_BGEV_SP2, what= "sigma", type = "response"); sigmahat
nuhat <- predict(fit_BGEV_SP2, what= "nu", type = "response"); nuhat
tauhat <- predict(fit_BGEV_SP2, what= "tau", type = "response"); tauhat

#**** Parameter standard error estimates via Delta method ****#

# --->  For sigma

#Derivate of h(theta) = exp(theta) (inverse of log function) is h'(theta) = exp(theta)

etahat <- predict(fit_BGEV_SP2, what= "sigma")
var.etahat <-  vcov(fit_BGEV_SP2)[kk1+1,kk1+1]
ep.sigmahat <- sqrt(((exp(etahat))^2)*var.etahat); ep.sigmahat


# ---> For nu

#Derivate of h(theta) = (exp(theta)-0.5)/(1 + exp(theta)) (inverse of own function)
# is g'(theta) = (1.5 * exp(theta)) / ((1 + exp(theta))^2)

etahat <- predict(fit_BGEV_SP2, what= "nu")
var.etahat <-  vcov(fit_BGEV_SP2)[kk1+2,kk1+2]
ep.nuhat <- sqrt((((1.5 * exp(etahat)) / ((1 + exp(etahat))^2))^2)*var.etahat); ep.nuhat

z_stat <- as.numeric(nuhat[1]/ep.nuhat[1]); z_stat
pvalue <- 2*pnorm(abs(z_stat), lower.tail=F); pvalue


# --->  For tau

etahat <- predict(fit_BGEV_SP2, what= "tau")
var.etahat <-  vcov(fit_BGEV_SP2)[kk1+3,kk1+3]
ep.tauhat <- sqrt(((exp(etahat))^2)*var.etahat); ep.tauhat


# *** 0.99 quantile fit regression ***#

fit_BGEV_q99_SP2 <- gamlss(formula = y~HUM+S,  
family = BGEV_q99(mu.link = "identity", sigma.link = "log", nu.link = "own", tau.link = "log"),
method = CG(), sigma.start = 1.2, control = gamlss.control(n.cyc = 200, trace = FALSE))
summary(fit_BGEV_q99_SP2)

#**** Parameter estimates ****#

sigmahat <- predict(fit_BGEV_q99_SP2, what= "sigma", type = "response"); sigmahat
nuhat <- predict(fit_BGEV_q99_SP2, what= "nu", type = "response"); nuhat
tauhat <- predict(fit_BGEV_q99_SP2, what= "tau", type = "response"); tauhat

#**** Parameter standard error estimates via Delta method ****#

# --->  For sigma

#Derivate of h(theta) = exp(theta) (inverse of log function) is h'(theta) = exp(theta)

obj <- summary(fit_BGEV_q99_SP2, type= "qr")

etahat <- predict(fit_BGEV_q99_SP2, what= "sigma")
ep.sigmahat <- sqrt((exp(etahat))^2)*as.numeric(obj[,2][kk1+1]); ep.sigmahat

# ---> For nu

#Derivate of h(theta) = (exp(theta)-0.5)/(1 + exp(theta)) (inverse of own function)
# is g'(theta) = (1.5 * exp(theta)) / ((1 + exp(theta))^2)

etahat <- predict(fit_BGEV_q99_SP2, what= "nu")
ep.nuhat <- sqrt(((1.5 * exp(etahat)) / ((1 + exp(etahat))^2))^2)*as.numeric(obj[,2][kk1+2]); ep.nuhat

z_stat <- as.numeric(nuhat[1]/ep.nuhat[1]); z_stat
pvalue <- 2*pnorm(abs(z_stat), lower.tail=F); pvalue


# --->  For tau

etahat <- predict(fit_BGEV_q99_SP2, what= "tau")
ep.tauhat <- sqrt((exp(etahat))^2)*as.numeric(obj[,2][kk1+3]); ep.tauhat



#**** Specification 3 : ~HUM+S+P ****#

#---- GEV Fitted model ----#

X <- model.matrix(~HUM+S+P)
kk1 <- ncol(X);  n <- nrow(X)
fitgev <-  suppressWarnings(tryCatch(vglm(y ~ X[,2:kk1], gevff(perc = NULL)), error=function(e) {e} )) #starting values via GEV regression
AIC_GEV_SP3 <- -2*logLik(fitgev) + 2*(kk1+2);  AIC_GEV_SP3


#---- BGEV Fitted model ----#

# *** Median fit - 0.50 quantile fit regression ***#

fit_BGEV_SP3 <- gamlss(formula = y~HUM+S+P,
family = BGEV(mu.link = "identity", sigma.link = "log", nu.link = "own", tau.link = "logit"),
method = CG(), control = gamlss.control(n.cyc = 200, trace = FALSE))
summary(fit_BGEV_SP3)

#*** Diagnostics ***#

par(mar=c(5.0,5.0,4.0,2.0))
wp(fit_BGEV_SP3, cex = 1.5, cex.lab = 2.0)
plot(fit_BGEV_SP3)

#--- qqboxplot ---#

rq_BGEV_SP3 <- resid(fit_BGEV_SP3) #quantile residuals

data <-  data.frame(y = rq_BGEV_SP3)
  ggplot(data, aes(y = rq_BGEV_SP3)) +                       
  geom_qqboxplot(notch=F, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=25),
        axis.title.x = element_text(size=25),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid = element_line(colour = "grey70"))

#**** Parameter estimates ****#

sigmahat <- predict(fit_BGEV_SP3, what= "sigma", type = "response"); sigmahat
nuhat <- predict(fit_BGEV_SP3, what= "nu", type = "response"); nuhat
tauhat <- predict(fit_BGEV_SP3, what= "tau", type = "response"); tauhat

#**** Parameter standard error estimates ****#

# --->  For sigma

#Derivate of h(theta) = exp(theta) (inverse of log function) is h'(theta) = exp(theta)

etahat <- predict(fit_BGEV_SP3, what= "sigma")
var.etahat <-  vcov(fit_BGEV_SP3)[kk1+1,kk1+1]
ep.sigmahat <- sqrt(((exp(etahat))^2)*var.etahat); ep.sigmahat

# ---> For nu

#Derivate of h(theta) = (exp(theta)-0.5)/(1 + exp(theta)) (inverse of own function)
# is g'(theta) = (1.5 * exp(theta)) / ((1 + exp(theta))^2)

etahat <- predict(fit_BGEV_SP3, what= "nu")
var.etahat <-  vcov(fit_BGEV_SP3)[kk1+2,kk1+2]
ep.nuhat <- sqrt((((1.5 * exp(etahat)) / ((1 + exp(etahat))^2))^2)*var.etahat); ep.nuhat

z_stat <- as.numeric(nuhat[1]/ep.nuhat[1]); z_stat
pvalue <- 2*pnorm(abs(z_stat), lower.tail=F); pvalue

# --->  For tau

#Derivate of h(theta) = exp(theta)/(1+exp(theta)) (inverse of logit function) is 
# h'(theta) = exp(theta)/((1+exp(theta))^2)

etahat <- predict(fit_BGEV_SP3, what= "tau")
var.etahat <-  vcov(fit_BGEV_SP3)[kk1+3,kk1+3]
ep.tauhat <- sqrt(((exp(etahat)/((1+exp(etahat))^2))^2)*var.etahat); ep.tauhat



# *** 0.99 quantile fit regression ***#

fit_BGEV_q99_SP3 <- gamlss(formula = y~HUM+S+P,  
family = BGEV_q99(mu.link = "identity", sigma.link = "log", nu.link = "own", tau.link = "logit"),
method = CG(), sigma.start = 1.3, control = gamlss.control(n.cyc = 200, trace = FALSE))
summary(fit_BGEV_q99_SP3)

#**** Parameter estimates ****#

sigmahat <- predict(fit_BGEV_q99_SP3, what= "sigma", type = "response"); sigmahat
nuhat <- predict(fit_BGEV_q99_SP3, what= "nu", type = "response"); nuhat
tauhat <- predict(fit_BGEV_q99_SP3, what= "tau", type = "response"); tauhat

#**** Parameter standard error estimates ****#

# --->  For sigma

#Derivate of h(theta) = exp(theta) (inverse of log function) is h'(theta) = exp(theta)

obj <- summary(fit_BGEV_q99_SP3, type= "qr")

etahat <- predict(fit_BGEV_q99_SP3, what= "sigma")
ep.sigmahat <- sqrt((exp(etahat))^2)*as.numeric(obj[,2][kk1+1]); ep.sigmahat

# ---> For nu

#Derivate of h(theta) = (exp(theta)-0.5)/(1 + exp(theta)) (inverse of own function)
# is g'(theta) = (1.5 * exp(theta)) / ((1 + exp(theta))^2)

etahat <- predict(fit_BGEV_q99_SP3, what= "nu")
ep.nuhat <- sqrt(((1.5 * exp(etahat)) / ((1 + exp(etahat))^2))^2)*as.numeric(obj[,2][kk1+2]); ep.nuhat

z_stat <- as.numeric(nuhat[1]/ep.nuhat[1]); z_stat
pvalue <- 2*pnorm(abs(z_stat), lower.tail=F); pvalue

# --->  For tau

#Derivate of h(theta) = exp(theta)/(1+exp(theta)) (inverse of logit function) is 
# h'(theta) = exp(theta)/((1+exp(theta))^2)

etahat <- predict(fit_BGEV_q99_SP3, what= "tau")
ep.tauhat <- sqrt((exp(etahat)/((1+exp(etahat))^2))^2)*as.numeric(obj[,2][kk1+3]); ep.tauhat

#--- qqboxplot ---#

rq_BGEV_SP3 <- resid(fit_BGEV_q99_SP3) #quantile residuals

data <-  data.frame(y = rq_BGEV_SP3)
  ggplot(data, aes(y = rq_BGEV_SP3)) +                       
  geom_qqboxplot(notch=F, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=25),
        axis.title.x = element_text(size=25),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid = element_line(colour = "grey70"))



#**** Specification 4 : ~HUM+S+P+WS ****#

#---- GEV Fitted model ----#

X <- model.matrix(~HUM+S+P+WS)
kk1 <- ncol(X);  n <- nrow(X)
fitgev <-  suppressWarnings(tryCatch(vglm(y ~ X[,2:kk1], gevff(perc = NULL)), error=function(e) {e} )) #starting values via GEV regression
AIC_GEV_SP4 <- -2*logLik(fitgev) + 2*(kk1+2);  AIC_GEV_SP4


#---- BGEV Fitted model ----#

# *** Median fit - 0.50 quantile fit regression ***#

fit_BGEV_SP4 <- gamlss(formula = y~HUM+S+P+WS,
family = BGEV(mu.link = "identity", sigma.link = "log", nu.link = "own", tau.link = "logit"),
method = CG(), control = gamlss.control(n.cyc = 200, trace = FALSE))
summary(fit_BGEV_SP4)

#*** Diagnostics ***#

par(mar=c(5.0,5.0,4.0,2.0))
wp(fit_BGEV_SP4, cex = 1.5, cex.lab = 2.0)
plot(fit_BGEV_SP4)

#--- qqboxplot ---#

rq_BGEV_SP4 <- resid(fit_BGEV_SP4) #quantile residuals

data <-  data.frame(y = rq_BGEV_SP4)
  ggplot(data, aes(y = rq_BGEV_SP4)) +                       
  geom_qqboxplot(notch=TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=25),
        axis.title.x = element_text(size=25),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid = element_line(colour = "grey70"))


#**** Parameter estimates ****#

sigmahat <- predict(fit_BGEV_SP4, what= "sigma", type = "response"); sigmahat
nuhat <- predict(fit_BGEV_SP4, what= "nu", type = "response"); nuhat
tauhat <- predict(fit_BGEV_SP4, what= "tau", type = "response"); tauhat

#**** Parameter standard error estimates ****#

# --->  For sigma

#Derivate of h(theta) = exp(theta) (inverse of log function) is h'(theta) = exp(theta)

etahat <- predict(fit_BGEV_SP4, what= "sigma")
var.etahat <-  vcov(fit_BGEV_SP4)[kk1+1,kk1+1]
ep.sigmahat <- sqrt(((exp(etahat))^2)*var.etahat); ep.sigmahat

# ---> For nu

#Derivate of h(theta) = (exp(theta)-0.5)/(1 + exp(theta)) (inverse of own function)
# is g'(theta) = (1.5 * exp(theta)) / ((1 + exp(theta))^2)

etahat <- predict(fit_BGEV_SP4, what= "nu")
var.etahat <-  vcov(fit_BGEV_SP4)[kk1+2,kk1+2]
ep.nuhat <- sqrt((((1.5 * exp(etahat)) / ((1 + exp(etahat))^2))^2)*var.etahat); ep.nuhat

z_stat <- as.numeric(nuhat[1]/ep.nuhat[1]); z_stat
pvalue <- 2*pnorm(abs(z_stat), lower.tail=F); pvalue

# --->  For tau

#Derivate of h(theta) = exp(theta)/(1+exp(theta)) (inverse of logit function) is 
# h'(theta) = exp(theta)/((1+exp(theta))^2)

etahat <- predict(fit_BGEV_SP4, what= "tau")
var.etahat <-  vcov(fit_BGEV_SP4)[kk1+3,kk1+3]
ep.tauhat <- sqrt(((exp(etahat)/((1+exp(etahat))^2))^2)*var.etahat); ep.tauhat



# *** 0.99 quantile fit regression ***#

fit_BGEV_q99_SP4 <- gamlss(formula = y~HUM+S+P+WS,  
family = BGEV_q99(mu.link = "identity", sigma.link = "log", nu.link = "own", tau.link = "logit"),
method = CG(), sigma.start = 1.2, control = gamlss.control(n.cyc = 200, trace = FALSE))
summary(fit_BGEV_q99_SP4)


#**** Parameter estimates ****#

sigmahat <- predict(fit_BGEV_q99_SP4, what= "sigma", type = "response"); sigmahat
nuhat <- predict(fit_BGEV_q99_SP4, what= "nu", type = "response"); nuhat
tauhat <- predict(fit_BGEV_q99_SP4, what= "tau", type = "response"); tauhat

#**** Parameter standard error estimates ****#

# --->  For sigma

#Derivate of h(theta) = exp(theta) (inverse of log function) is h'(theta) = exp(theta)

etahat <- predict(fit_BGEV_q99_SP4, what= "sigma")
var.etahat <-  vcov(fit_BGEV_q99_SP4)[kk1+1,kk1+1]
ep.sigmahat <- sqrt(((exp(etahat))^2)*var.etahat); ep.sigmahat

# ---> For nu

#Derivate of h(theta) = (exp(theta)-0.5)/(1 + exp(theta)) (inverse of own function)
# is g'(theta) = (1.5 * exp(theta)) / ((1 + exp(theta))^2)

etahat <- predict(fit_BGEV_q99_SP4, what= "nu")
var.etahat <-  vcov(fit_BGEV_q99_SP4)[kk1+2,kk1+2]
ep.nuhat <- sqrt((((1.5 * exp(etahat)) / ((1 + exp(etahat))^2))^2)*var.etahat); ep.nuhat

z_stat <- as.numeric(nuhat[1]/ep.nuhat[1]); z_stat
pvalue <- 2*pnorm(abs(z_stat), lower.tail=F); pvalue

# --->  For tau

#Derivate of h(theta) = exp(theta)/(1+exp(theta)) (inverse of logit function) is 
# h'(theta) = exp(theta)/((1+exp(theta))^2)

etahat <- predict(fit_BGEV_q99_SP4, what= "tau")
var.etahat <-  vcov(fit_BGEV_q99_SP4)[kk1+3,kk1+3]
ep.tauhat <- sqrt(((exp(etahat)/((1+exp(etahat))^2))^2)*var.etahat); ep.tauhat


#--- qqboxplot ---#

rq_BGEV_SP4 <- resid(fit_BGEV_q99_SP4) #quantile residuals

data <-  data.frame(y = rq_BGEV_SP4)
  ggplot(data, aes(y = rq_BGEV_SP4)) +                       
  geom_qqboxplot(notch=TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=25),
        axis.title.x = element_text(size=25),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid = element_line(colour = "grey70"))


#*** Histograma of ordinary and quantile residuals under median regression fit ***#

# Specification 1

ordinary_res0 <- y-fitted(fit_BGEV0_SP1)
ordinary_res <- y-fitted(fit_BGEV_SP1)

# BGEV with xi=0
df <- data.frame(x=ordinary_res0)

 ggplot(df,aes(x=x))+
    geom_histogram(aes(y = ..density..),colour='white',fill='#696969', bins=7)+
    geom_density(lwd=2,linetype = 1,colour = 'blue')+
    labs(x='untransformed residual',y='', color='')+
    scale_color_manual(values=colors)+
    theme_bw()+
    theme(axis.title.y=element_text(colour='black',size=25),
          axis.title.x=element_text(colour='black',size=25),
          axis.text=element_text(colour='black',size=25),
          panel.border=element_blank(),
          axis.line=element_line(colour='black'))


# BGEV with xi neq 0

df <- data.frame(x=ordinary_res)

 ggplot(df,aes(x=x))+
    geom_histogram(aes(y = ..density..),colour='white',fill='#696969', bins=7)+
    geom_density(lwd=2,linetype = 1,colour = 'blue')+
    labs(x='untransformed residual',y='', color='')+
    scale_color_manual(values=colors)+
    theme_bw()+
    theme(axis.title.y=element_text(colour='black',size=25),
          axis.title.x=element_text(colour='black',size=25),
          axis.text=element_text(colour='black',size=25),
          panel.border=element_blank(),
          axis.line=element_line(colour='black'))



# Specification 2

ordinary_res0 <- y-fitted(fit_BGEV0_SP2)
ordinary_res <- y-fitted(fit_BGEV_SP2)

# BGEV with xi=0
df <- data.frame(x=ordinary_res0)

 ggplot(df,aes(x=x))+
    geom_histogram(aes(y = ..density..),colour='white',fill='#696969', bins=7)+
    geom_density(lwd=2,linetype = 1,colour = 'blue')+
    labs(x='untransformed residual',y='', color='')+
    scale_color_manual(values=colors)+
    theme_bw()+
    theme(axis.title.y=element_text(colour='black',size=25),
          axis.title.x=element_text(colour='black',size=25),
          axis.text=element_text(colour='black',size=25),
          panel.border=element_blank(),
          axis.line=element_line(colour='black'))


# BGEV with xi neq 0

df <- data.frame(x=ordinary_res)

 ggplot(df,aes(x=x))+
    geom_histogram(aes(y = ..density..),colour='white',fill='#696969', bins=7)+
    geom_density(lwd=2,linetype = 1,colour = 'blue')+
    labs(x='untransformed residual',y='', color='')+
    scale_color_manual(values=colors)+
    theme_bw()+
    theme(axis.title.y=element_text(colour='black',size=25),
          axis.title.x=element_text(colour='black',size=25),
          axis.text=element_text(colour='black',size=25),
          panel.border=element_blank(),
          axis.line=element_line(colour='black'))


# Specification 3

ordinary_res0 <- y-fitted(fit_BGEV0_SP3)
ordinary_res <- y-fitted(fit_BGEV_SP3)


# BGEV with xi=0
df <- data.frame(x=ordinary_res0)

 ggplot(df,aes(x=x))+
    geom_histogram(aes(y = ..density..),colour='white',fill='#696969', bins=7)+
    geom_density(lwd=2,linetype = 1,colour = 'blue')+
    labs(x='untransformed residual',y='', color='')+
    scale_color_manual(values=colors)+
    theme_bw()+
    theme(axis.title.y=element_text(colour='black',size=25),
          axis.title.x=element_text(colour='black',size=25),
          axis.text=element_text(colour='black',size=25),
          panel.border=element_blank(),
          axis.line=element_line(colour='black'))


# BGEV with xi neq 0

df <- data.frame(x=ordinary_res)

 ggplot(df,aes(x=x))+
    geom_histogram(aes(y = ..density..),colour='white',fill='#696969', bins=7)+
    geom_density(lwd=2,linetype = 1,colour = 'blue')+
    labs(x='untransformed residual',y='', color='')+
    scale_color_manual(values=colors)+
    theme_bw()+
    theme(axis.title.y=element_text(colour='black',size=25),
          axis.title.x=element_text(colour='black',size=25),
          axis.text=element_text(colour='black',size=25),
          panel.border=element_blank(),
          axis.line=element_line(colour='black'))


# Specification 4


ordinary_res0 <- y-fitted(fit_BGEV0_SP4)
ordinary_res <- y-fitted(fit_BGEV_SP4)


# BGEV with xi=0
df <- data.frame(x=ordinary_res0)

 ggplot(df,aes(x=x))+
    geom_histogram(aes(y = ..density..),colour='white',fill='#696969', bins=7)+
    geom_density(lwd=2,linetype = 1,colour = 'blue')+
    labs(x='untransformed residual',y='', color='')+
    scale_color_manual(values=colors)+
    theme_bw()+
    theme(axis.title.y=element_text(colour='black',size=25),
          axis.title.x=element_text(colour='black',size=25),
          axis.text=element_text(colour='black',size=25),
          panel.border=element_blank(),
          axis.line=element_line(colour='black'))


# BGEV with xi neq 0

df <- data.frame(x=ordinary_res)

 ggplot(df,aes(x=x))+
    geom_histogram(aes(y = ..density..),colour='white',fill='#696969', bins=7)+
    geom_density(lwd=2,linetype = 1,colour = 'blue')+
    labs(x='untransformed residual',y='', color='')+
    scale_color_manual(values=colors)+
    theme_bw()+
    theme(axis.title.y=element_text(colour='black',size=25),
          axis.title.x=element_text(colour='black',size=25),
          axis.text=element_text(colour='black',size=25),
          panel.border=element_blank(),
          axis.line=element_line(colour='black'))


# **** Ilustration of modelling all parameters -- Median regression  ****#


#---- BGEV Fitted model ----#

fit_BGEV_SP4_S <- gamlss(formula = y~HUM+S+P+WS, 
sigma.formula=~S, nu.formula=~S, tau.formula=~S,
family = BGEV(mu.link = "identity", sigma.link = "log", nu.link = "own", tau.link = "logit"),
method = CG(), control = gamlss.control(n.cyc = 200, trace = FALSE))
summary(fit_BGEV_SP4_S)


#*** Diagnostics ***#

par(mar=c(5.0,5.0,4.0,2.0))
wp(fit_BGEV_SP4_S, cex = 1.5, cex.lab = 2.0)
plot(fit_BGEV_SP4)
