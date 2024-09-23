# INITIALISE TO LOCAL CONTEXT
rm(list=ls())
# set work directory
# setwd("~/OneDrive - Yale University/RSV/CAR/20000 iter/data")
# READ DATA
NYwide<- readRDS("/Volumes/HCUPdata-CC0941-MEDSPH/Rcode_Ke/rjags/rjag_rsv_spatiotemporal/data/NYwide.rds") # monthly time series of observed RSV hospitalizations from each age group in NY
nypopulation <- readRDS("/Volumes/HCUPdata-CC0941-MEDSPH/Rcode_Ke/rjags/rjag_rsv_spatiotemporal/data/nypopulation.rds") #offset term, population for each age group in the area (?matrix)

#reporting_fraction <- readRDS("reporting_fraction.rds") # offset term, recording rates for each age group in NY (?matrix)
reporting_fraction  <- c(0.68 * 1/4 + 0.44 * 1/4 + 0.3 * 1/2, 
                         0.12 * 5/15 + 0.07 * 10/15, 
                         0.04 * 25/45 + 0.11 * 20/45, 
                         0.02 * 20/35 + 0.01 * 15/35)

NYwide <- t(NYwide)
nypopulation <- t(nypopulation)[,1]

names(NYwide) = NULL
names(nypopulation) = NULL

set.seed(1234)
#################################################################
# First stage model 
# NOTE precision = 0.01 --> var = 100. 
# using weakly informative priors 
#################################################################

firststage <- function(y, logpop, agegroup, date, n, m, fraction) {
  library(rjags)
  library(coda)
  
  model_string <- "
model {
  for (i in 1:n.agegroup) {
    for (j in 1:n.date) {
      y[i,j]  ~ dpois(lambda[i,j])
      log(lambda[i,j]) <- logpop[i] + beta0[i] + amp[i]*cos(2*pi*j/12 - phase[i]) + phi[i,j] + log(fraction[i])
      # harmonic Poisson regression model
      phi[i,j] ~ dnorm(0, tau3[i])
    }
    
    phase[i] <- (2*pi*(exp(trans[i]))) / (1 + exp(trans[i])) # support estimates in real line
    beta0[i] ~ dnorm(0, 0.01) # priors
    amp[i] <- exp(beta1[i])
    beta1[i] ~ dnorm(0, 0.01)
    trans[i] ~ dnorm(0, 0.01)
    tau3[i] ~ dgamma(0.01, 0.01)
  }
}"
  
  # Prepare the dataset
  dataset <- list("y" = y, "logpop" = logpop, "n.agegroup" = dim(agegroup)[1], "n.date" = dim(date)[2], pi = pi, "fraction" = fraction)
  
  # Initialize the model
  model <- jags.model(textConnection(model_string),
                      data = dataset,
                      n.chains = 1)
  
  # Burn-in phase
  update(model, n.iter = n)
  
  #trans_variables <- paste0("trans[", 1:4, "]")
  
  # Posterior sampling
  modelresult <- coda.samples(model, variable.names = c("trans"), thin = 10, n.iter = m)
  
  # Convert results to a dataframe
  # result <- as.data.frame(modelresult[[1]])
  
  model.results <<- modelresult
  # Store result in a global variable
  # first.result <<- result
  
}


# run model for NY
firststage(y = NYwide,
           logpop = log(nypopulation),
           agegroup = NYwide,
           date = NYwide,
           n = 10000,
           m = 50000,
           fraction = reporting_fraction)


# extract estimates

result.trans <- as.data.frame(model.results[[1]])
firstny <- first.result
nytransmean<- colMeans(first.result)  # posterior mean of phase estimates
var.ny <- c()
for (i in 1:1745) {
  var.ny[i] <- var(first.result[,i])
  } # calculate variations of phase estimates

 
phase <- (2*pi*(exp(first.result))) / (1+exp(first.result)) # support estimates in real line




saveRDS(firstny,"firstny.rds")  # save posterior estimates for the second stage model
saveRDS(var.ny,"varny.rds")           # save variations for the second stage
saveRDS(first.result,"transny.rds")   # save phase estimates for the second stage
