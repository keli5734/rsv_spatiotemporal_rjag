
data {
  int<lower=1> n_agegroup;  // Number of age groups
  int<lower=1> n_date;      // Number of time points (dates)
  int<lower=0> y[n_agegroup, n_date];  // Observed count data (must be integers)
  vector[n_agegroup] logpop;     // Log population for each age group
}

parameters {
  vector[n_agegroup] beta0;          // Intercept for each age group
  vector[n_agegroup] beta1;          // Parameter for amplitude (log-scale)
  vector[n_agegroup] trans;          // Transformed phase parameter
  matrix[n_agegroup, n_date] phi;    // Time-varying random effects
  vector<lower=0>[n_agegroup] tau3;  // Standard deviation for phi
}

transformed parameters {
  matrix[n_agegroup, n_date] lambda;  // Poisson rate parameter (intensity)
  vector[n_agegroup] amp;             // Amplitude for each age group
  vector[n_agegroup] phase;           // Phase for each age group

  for (i in 1:n_agegroup) {
    amp[i] = exp(beta1[i]);  // Transform amplitude to positive scale
    phase[i] = (2 * pi() * exp(trans[i])) / (1 + exp(trans[i]));  // Phase transformation

    for (j in 1:n_date) {
      lambda[i, j] = exp(logpop[i] + beta0[i] + amp[i] * cos(2 * pi() * j / 12 - phase[i]) + phi[i, j]);
    }
  }
}

model {
  // Priors
  beta0 ~ normal(0, 5);     // Prior for beta0
  beta1 ~ normal(0, 10);    // Prior for beta1
  trans ~ normal(0, 10);    // Prior for trans
  tau3 ~ gamma(0.01, 0.01); // Prior for tau3 (Gamma distribution)

  // Likelihood
  for (i in 1:n_agegroup) {
    phi[i] ~ normal(0, tau3[i]);  // Prior for phi (random effect)
    for (j in 1:n_date) {
      y[i, j] ~ poisson(lambda[i, j]);  // Poisson likelihood
    }
  }
}
