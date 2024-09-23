# INITIALISE TO LOCAL CONTEXT
rm(list=ls())


library(tidyverse)

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
      log(lambda[i,j]) <- logpop[i] + beta0[i] + amp[i]*cos(2*pi*j/12 - phase[i]) + phi[i,j] - log(fraction[i])
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
                      n.chains = 4)
  
  # Burn-in phase
  update(model, n.iter = n)
  
  # Posterior sampling
  modelresult <- coda.samples(model, variable.names = c("trans"), thin = 10, n.iter = m)
  
  model.results <<- modelresult
 
  
}


# run model for NY
firststage(y = NYwide,
           logpop = log(nypopulation),
           agegroup = NYwide,
           date = NYwide,
           n = 10000,
           m = 50000,
           fraction = reporting_fraction)

# check convergence
gelman.diag(model.results)
gelman.plot(model.results)
par(mfrow=c(2,2))
traceplot(model.results)
# extract estimates

result.trans <- mcmc.list(model.results)

# Assuming 'combined_result' is your mcmc.list object
# Convert mcmc.list to a matrix
pooled_matrix <- as.matrix(result.trans)

# Convert the matrix to a data.frame
pooled_df <- as.data.frame(pooled_matrix)

# Now you can access the 'trans' parameters or other variables
result.trans <- pooled_df[, grep("trans", colnames(pooled_df))]

nytransmean<- colMeans(result.trans)  # posterior mean of phase estimates
# var.ny <- c()
# for (i in 1:4) {
#   var.ny[i] <- var(result.trans[,i])
#   } # calculate variations of phase estimates

 
phase <- (2*pi*(exp(result.trans))) / (1+exp(result.trans))  
shift_in_months <- phase/2/pi *12

#saveRDS(shift_in_months, "shift_in_months.rds")

phase_long <- shift_in_months %>%
  pivot_longer(
    cols = starts_with("trans"),
    names_to = "AGE_GROUP",
    values_to = "T"
  ) %>% 
  mutate(AGE_GROUP = ifelse(AGE_GROUP == "trans[1]", "under 5", 
                            ifelse(AGE_GROUP == "trans[2]", "5-19", 
                                   ifelse(AGE_GROUP == "trans[3]", "20-59", "60 above")))) %>% 
  mutate(AGE_GROUP = factor(AGE_GROUP, levels = c("under 5", "5-19", "20-59", "60 above")))

 
# Plot histogram
ggplot(phase_long, aes(x = T, fill = AGE_GROUP)) +
  geom_histogram(binwidth = 0.01, alpha = 0.7) +
  labs(title = "Histogram of Phase Values",
       x = "Phase Values",
       y = "Frequency") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") + 
  facet_wrap(~AGE_GROUP)





  phase_long %>% 
  ggplot(aes(x = AGE_GROUP, y = T, fill = AGE_GROUP)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width=.1, outlier.colour=NA,  aes(group = AGE_GROUP, fill = AGE_GROUP)) +
  #lims(y = c(3,5)) +
  theme_bw() +
  scale_fill_manual(name="AGE GROUP",
                    values=c("under 5" = "#1f78b4",
                             "5-19"  = "#fdbf6f",
                             "20-59" = "#fb9a99",
                             "60 above" = "#a6d96a"),
                    guide = "none") + 
  ylab("Shift in peak timing (months)") +
  xlab("") +
  theme(axis.text.x = element_text(color="black",
                                   size = 15, angle=0),
        axis.text.y = element_text(color="black",
                                   size= 15, angle=0),
        text = element_text(size = 15)) + 
  theme(
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.x = element_text(face="bold"),
    axis.title.y = element_text(face="bold"),
    plot.title = element_text(face="bold", hjust = 0.5),
    legend.position = "NA",
    legend.title = element_text(size = 0))+
  theme(strip.text = element_text(face = "bold"))

# T captures the shift in peak timing (with the 12-month period starting in early July)

  
  


