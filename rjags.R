# INITIALISE TO LOCAL CONTEXT
rm(list=ls())
# set work directory
# setwd("~/OneDrive - Yale University/RSV/CAR/20000 iter/data")
# READ DATA
NYwide<- readRDS("NYwide.rds") # monthly time series of observed RSV hospitalizations from each age group in NY
nyunder5 <- readRDS("nyunder5.rds") #offset term, population for each age group in NY (?matrix)
reporting_fraction <- readRDS("reporting_fraction.rds") # offset term, recording rates for each age group in NY (?matrix)

#################################################################
#First stage model
#################################################################
firststage <- function(y,logpop,zip,date,n,m){
  library(rjags)
  library(coda)
  
  model_string<-"
model{

for(i in 1:n.zip){
for(j in 1:n.date){
y[i,j]  ~ dpois(lambda[i,j])

log(lambda[i,j])<- logpop[i]+ beta0[i] + amp[i]*cos(2*pi*j/12 - phase[i]) + phi[i,j]
# harmonic Poisson regression model

phi[i,j] ~ dnorm (0,tau3[i])
}

phase[i] <- (2*pi*(exp(trans[i]))) / (1+exp(trans[i])) # support estimates in real line
beta0[i] ~ dnorm(0, 0.001) # priors
amp[i]<-exp(beta1[i])
beta1[i] ~ dnorm(0,0.001) 
trans[i] ~ dnorm(0,0.001)
tau3[i] ~ dgamma(0.01, 0.01)} 
}"

dataset <- list('y' = y,"logpop"=logpop, "n.zip"=zip,"n.date"=date,pi=pi)

model<-jags.model(textConnection(model_string),
                  data=dataset,
                  n.chains=1)
update(model, 
       n.iter=n) # burn-in 

modelresult<-coda.samples(model, variable.names=c("trans"),
                          thin=10,
                          n.iter=m) # posterior samples
result <- as.data.frame(modelresult[[1]])
first.result<<-result
}

#run model for NY
firststage(y=NYwide,logpop=nyunder5$log,zip=nrow(NYwide),date=ncol(NYwide),n=10000,m=50000)
firstny <- first.result
nytransmean<- colMeans(first.result)  # posterior mean of phase estimates
var.ny <- c()
for (i in 1:1745) {
  var.ny[i] <- var(first.result[,i])} # calculate variations of phase estimates
saveRDS(firstny,"firstny.rds")  # save posterior estimates for the second stage model
saveRDS(var.ny,"varny.rds")           # save variations for the second stage
saveRDS(first.result,"transny.rds")   # save phase estimates for the second stage
