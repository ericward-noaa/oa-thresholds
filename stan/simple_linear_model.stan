data {
  int<lower=1> n_obs; // number of observations
  vector[n_obs] logRS; // real observations
  vector[n_obs] S; // spawners
  vector[n_obs] env; // environmental covariate
}
parameters {
  real<lower=0> sigma_resid;
  real<upper=0> b_s; // density dep effect has to be negative
  real<lower=0> b_0; // constrain intercept to be positive
  real slope;
}
transformed parameters {
  vector[n_obs] pred;
  for(i in 1:n_obs) {
    pred[i] = b_0 + b_s*S[i] + slope*env[i];
  }
}
model {
  b_0 ~ normal(0,1);
  b_s ~ normal(0,1);
  slope ~ normal(0,1);
  sigma_resid ~ student_t(3,0,2);
  logRS ~ normal(pred, sigma_resid);
}
generated quantities {
  // log_lik stores log liklihood for looic
  vector[n_obs] log_lik;
  for(i in 1:n_obs) {
    log_lik[i] = normal_lpdf(logRS[i] | pred[i], sigma_resid);
  }
}
