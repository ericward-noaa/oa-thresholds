data {
  int<lower=1> n_obs; // number of observations
  vector[n_obs] logRS; // real observations
  vector[n_obs] S; // spawners
  vector[n_obs] env; // environmental covariate
  int n_breakpoints; // 1 or 2 for now
}
parameters {
  real<lower=0> sigma_resid;
  ordered[n_breakpoints] cutpoints; // on the scale of z-scored env data
  //vector[n_breakpoints+1] slopes; // estimated slopes between segments
  real<lower=0> slope_1;
  real<upper=0> slope_2;
  real<upper=0> b_s; // density dep effect has to be negative
  real<lower=0> b_0; // constrain intercept to be positive
  real b0; // constrain intercept to be positive
}
transformed parameters {
  vector[n_obs] effect;
  vector[n_obs] pred;
  vector[2] slopes;
  slopes[1] = slope_1;
  slopes[2] = slope_2;
  
  for(i in 1:n_obs) {
    if(env[i] < cutpoints[1]) {
      effect[i] = b0 + slopes[1]*env[i];
    } else {
      // 2 or 3 segments
      if(n_breakpoints==1) {
        effect[i] = (b0 + slopes[1]*cutpoints[1]) + slopes[2]*env[i];
      } else {
        if(env[i] < cutpoints[2]) {
          //effect[i] = (b0 + slopes[1]*cutpoints[1]) + slopes[2]*env[i];
          effect[i] = (b0 + slopes[1]*cutpoints[1]);
        } else {
          //effect[i] = (slopes[2]*cutpoints[2]) + slopes[3]*env[i];
          effect[i] = (b0 + slopes[1]*cutpoints[1]) + slopes[2]*env[i];
        }
      }
    }
    pred[i] = b_0 + b_s*S[i] + effect[i];
  }
  
}
model {
  b_0 ~ normal(0,1);
  b_s ~ normal(0,1);
  //slopes ~ normal(0,1);
  slope_1 ~ normal(0,10);
  slope_2 ~ normal(0,10);
  sigma_resid ~ student_t(3,0,2);
  if(n_breakpoints==1) {
    cutpoints[1] ~ normal(0,3);
  } else {
    cutpoints[1] ~ normal(-1,3);
    cutpoints[2] ~ normal(1,3);
  }
  logRS ~ normal(pred, sigma_resid);
}
generated quantities {
  // log_lik stores log liklihood for looic
  vector[n_obs] log_lik;
  for(i in 1:n_obs) {
    log_lik[i] = normal_lpdf(logRS[i] | pred[i], sigma_resid);
  }
}