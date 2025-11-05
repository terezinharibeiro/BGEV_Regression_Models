pacman::p_load(fExtremes, tidyverse, bgev, EnvStats, DEoptim, stats)

# Function with correction of bgev package to obtain MLE estimates for BGEV parameters 
bgev.mle.CORRECTED = function(x, itermax = 500, method_envstats = "mle") 
{
  likbgev2 = function(theta) {
    val = -likbgev(x, theta)
    if (val == Inf) 
      return(1e+100)
    else return(val)
  }
  fit_gev <- egevd(x = x, method = method_envstats)
  mu_start = fit_gev$parameters[1]
  sd_start = fit_gev$parameters[2]
  xi_start = -fit_gev$parameters[3]
  starts = c(mu_start, sd_start, xi_start, 0.1)
  xi_min = xi_start/5
  xi_max = xi_start*5
  if(xi_start < 0){
    xi_max = xi_start/5
    xi_min = xi_start*5
  }
  lower = c(mu_start - 2*sd_start, sd_start/5, xi_min,-0.99)
  upper = c(mu_start + 2*sd_start, 5*sd_start, xi_max,5)
  starts.DEoptim = DEoptim(fn = likbgev2, lower, upper, control = DEoptim.control(itermax = itermax, 
                                                                                  trace = FALSE))
  
  esti <- optim(par = starts.DEoptim$optim$bestmem, fn = likbgev2, 
                method = "L-BFGS-B", lower = lower, upper = upper)
  esti
}

DTP <- c(15.6958333, 14.8708333,  9.4666667,  2.8083333,  0.8041667, 12.6916667, 15.3083333, 14.7625000,
         12.8958333,  6.9875000,  5.7875000, 10.3250000, 16.4958333, 15.6333333, 12.8416667,  7.0416667, 4.1166667,
         12.7125000, 14.1958333, 17.1916667, 10.1708333,  4.4500000,  4.3958333,  6.4125000, 13.6916667, 16.5875000,
         14.0833333,  7.0041667,  5.6750000,  6.2791667, 14.8333333, 16.2958333, 12.1916667,  5.4500000,  5.0625000,
         10.4416667, 13.6041667, 15.4500000, 10.2083333,  5.4125000, 2.5375000,  4.7000000, 13.5291667, 15.3166667,
         9.2625000,  8.1833333,  1.4541667, 11.2166667, 13.6000000, 13.4625000, 15.1166667,  4.2375000,  4.2375000,
         3.8333333, 15.9875000, 17.4041667, 12.4958333,  7.6875000,  1.4041667,  3.2875000, 13.0291667, 15.6041667,
         11.9208333,  2.4375000, 0.5291667,  4.3666667, 16.1791667, 14.5458333, 11.8041667,  3.7000000,  3.8500000,
         4.1125000, 12.5166667)


# GEV fit

fitGEV <- gevFit(DTP, type ="mle")
summary(fitGEV)

# BGEV fit

fitBGEV <- bgev.mle.CORRECTED(DTP)
pars <- fitBGEV$par
names(pars) = c("mu", "sigma", "xi", "delta")
print(pars)


# Histogram of DTP with fitted densities GEV (red) and BGEV (blue) distributions

x <- DTP
params <- fitGEV@fit[["par.ests"]]
curve_data <- data.frame(
  x = seq(min(x), max(x), length.out = 500),
  density = fExtremes::dgev(seq(min(x), max(x), length.out = 500), 
                            xi = params[["xi"]], mu = params[["mu"]], beta = params[["beta"]]))

sequencia <- seq(min(DTP), max(DTP), length.out = length(DTP))
curvas <- data.frame(
  x = rep(sequencia, 1),
  densidade = c(
    bgev::dbgev(sequencia, mu = 8.7178821, sigma = 25.0543607, xi = -0.3751647, delta = 0.8609493)
  ),
  distrib = factor(rep(c("μ = 10,0\n\nσ = 1,00\n\nξ = 0,20\n\nδ = 0,50\n"), 
                       each = length(DTP)))
)
(graf_bgevum <- ggplot() +
    geom_histogram(aes(x = DTP, 
                       y = after_stat(density)), 
                   bins=10, fill = "gray60", color = "white") +
    geom_line(data = curvas, aes(x = x, y = densidade, color = distrib),
              color = "#004e89",linewidth = 1)+
    geom_line(data = curve_data, aes(x, density), color = "red", size = 1) +
    labs(x = "min",y = "") +
    theme_bw())
