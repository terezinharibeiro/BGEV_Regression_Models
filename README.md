# R codes and dataset for the paper "A Regression-Type Model for Bimodal Extreme-Valued Data"

This repository contains R codes and dataset used in the application and simulations of the paper "A Regression-Type Model for Bimodal Extreme-Valued Data" by Otiniano, Lisboa, Ribeiro and Fachini-Gomes (2025).

### Application

The directory Application contains the dataset and the R scripts to replicate the results presented in the application section of the paper.

- Application_BGEV.R: R script to replicate the inference and diagnostics results for the DTP dataset. 
- BGEV_GAMLSS.R:  R function that fits BGEV regression models with xi different from zero using the  gamlss package.
- BGEV0_GAMLSS.R: R function that fits BGEV regression models with xi equal to zero using the gamlss package.
- DTPdata.rds: File with DTP dataset.
- IID_fit.R: R script to replicate the parameter estimation for GEV and BGEV distributions in the case of IID data.

### Simulations

The directory Simulations contains folders corresponding to the simulation scenarios presented in the Monte Carlo simulation results section of the paper. For example, folder
Scenario1 contains the following files: 

- Simulation_Scenario1.R: R script to generate the results and files.
- BGEV0_GAMLSS.R:  R function that fits BGEV regression models with xi equal to zero using the gamlss package
