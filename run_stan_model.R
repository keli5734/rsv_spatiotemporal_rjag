library(rstan)

# Example data (replace with your actual data)
data_list <- list(
  n_agegroup = 4,  # Number of age groups
  n_date = 12,     # Number of dates
  y = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24), 4, 6),  # Example count data
  logpop = rep(8, 4)     # Example log population sizes (replace with actual values)
)

# Compile the Stan model
fit <- stan(model_code = stan_model_string, data = data_list, chains = 4, iter = 2000, warmup = 1000, thin = 1)

# Extract results
print(fit)