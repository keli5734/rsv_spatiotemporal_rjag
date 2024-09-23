 


# INITIALISE TO LOCAL CONTEXT
rm(list=ls())
library(rstan)
# set work directory
# setwd("~/OneDrive - Yale University/RSV/CAR/20000 iter/data")
# READ DATA
NYwide<- readRDS("/Volumes/HCUPdata-CC0941-MEDSPH/Rcode_Ke/rjags/rjag_rsv_spatiotemporal/data/NYwide.rds") # monthly time series of observed RSV hospitalizations from each age group in NY
nypopulation <- readRDS("/Volumes/HCUPdata-CC0941-MEDSPH/Rcode_Ke/rjags/rjag_rsv_spatiotemporal/data/nypopulation.rds") #offset term, population for each age group in the area (?matrix)

#reporting_fraction <- readRDS("reporting_fraction.rds") # offset term, recording rates for each age group in NY (?matrix)

NYwide <- t(NYwide)
nypopulation <- t(nypopulation)[,1]


# Example data (replace with your actual data)
data_list <- list(
  n_agegroup = 4,  # Number of age groups
  n_date = 120,     # Number of dates
  y = NYwide,
  logpop = log(nypopulation)     # Example log population sizes (replace with actual values)
)



# Example initial values function for each chain
init_function <- function() {
  list(
    beta0 = rnorm(4, 0, 1),    # 4 values for 4 age groups
    beta1 = rnorm(4, 0, 1),    # Initial values for beta1
    trans = rnorm(4, 0, 1),    # Initial values for trans
    phi = matrix(rnorm(4 * 12, 0, 0.1), 4, 120),  # 4 age groups, 12 time points (dates)
    tau3 = rgamma(4, 0.01, 0.01)  # Initial values for tau3
  )
}



model_fit <- stan(file = "rstan_model.stan",
                      data = data_list,
                      seed = 20240923,  # set random seed for reproducibility
                      pars = c("trans"),
                      iter = 4000,
                      chains = 1,
                      init = list(init_function()),
                      warmup = 2000,
                      control = list(adapt_delta = 0.99, max_treedepth = 25))



 
# Extract results
print(model_fit)