# INITIALISE TO LOCAL CONTEXT
rm(list=ls())
# set work directory
# setwd("~/OneDrive - Yale University/RSV/CAR/20000 iter/data")
# READ DATA
NYwide<- readRDS("/Volumes/HCUPdata-CC0941-MEDSPH/Rcode_Ke/rjags/rjag_rsv_spatiotemporal/data/NYwide.rds") # monthly time series of observed RSV hospitalizations from each age group in NY
nypopulation <- readRDS("/Volumes/HCUPdata-CC0941-MEDSPH/Rcode_Ke/rjags/rjag_rsv_spatiotemporal/data/nypopulation.rds") #offset term, population for each age group in the area (?matrix)

#reporting_fraction <- readRDS("reporting_fraction.rds") # offset term, recording rates for each age group in NY (?matrix)

NYwide <- t(NYwide)
nypopulation <- t(nypopulation)[,1]

set.seed(1234)
#################################################################
#First stage model
#################################################################
firststage <- function(y, logpop, agegroup, date, n, m) {
  library(rjags)
  library(coda)
  
  model_string <- "
model {
  for (i in 1:n.agegroup) {
    for (j in 1:n.date) {
      y[i,j]  ~ dpois(lambda[i,j])
      log(lambda[i,j]) <- logpop[i] + beta0[i] + amp[i]*cos(2*pi*j/12 - phase[i]) + phi[i,j]
      # harmonic Poisson regression model
      phi[i,j] ~ dnorm(0, tau3[i])
    }
    
    phase[i] <- (2*pi*(exp(trans[i]))) / (1 + exp(trans[i])) # support estimates in real line
    beta0[i] ~ dnorm(0, 5) # priors
    amp[i] <- exp(beta1[i])
    beta1[i] ~ dnorm(0, 10)
    trans[i] ~ dnorm(0, 10)
    tau3[i] ~ dgamma(0.01, 0.01)
  }
}"
  
  # Prepare the dataset
  dataset <- list('y' = y, "logpop" = logpop, "n.agegroup" = length(agegroup), "n.date" = length(date), pi = pi)
  
  # Initialize the model
  model <- jags.model(textConnection(model_string),
                      data = dataset,
                      n.chains = 5)
  
  # Burn-in phase
  update(model, n.iter = n)
  
  trans_variables <- paste0("trans[", 1:4, "]")
  
  # Posterior sampling
  modelresult <- coda.samples(model, variable.names = trans_variables, thin = 10, n.iter = m)
  
  # Convert results to a dataframe
  result <- as.data.frame(modelresult[[5]])
  
  model.results <<- modelresult
  # Store result in a global variable
  first.result <<- result
  
}


#run model for NY
firststage(y=NYwide,logpop=log(nypopulation),agegroup=nrow(NYwide),date=ncol(NYwide),n=10000,m=50000)

firstny <- first.result
nytransmean<- colMeans(first.result)  # posterior mean of phase estimates
var.ny <- c()
for (i in 1:1745) {
  var.ny[i] <- var(first.result[,i])} # calculate variations of phase estimates

saveRDS(firstny,"firstny.rds")  # save posterior estimates for the second stage model
saveRDS(var.ny,"varny.rds")           # save variations for the second stage
saveRDS(first.result,"transny.rds")   # save phase estimates for the second stage

phase4 <- (2*pi*(exp(first.result))) / (1+exp(first.result)) # support estimates in real line

