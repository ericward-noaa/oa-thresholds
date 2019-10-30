
library(rstan)

# generate random data
set.seed(123)
N = 40
env = cumsum(rnorm(N, 0, 1)) # env covariate
env = scale(env)
S = exp(rnorm(N, 5)) # spawners
effect = 200 - 20*env*env# quadratic effect of env variable
plot(env, effect) 

logRS = rnorm(N, mean = 3 + effect - 0.1*S, sd = 0.01)
plot(S, logRS)

# try to fit 1-breakpoint model in Stan
data_list = list(n_obs=N, logRS=logRS, 
  S=S, env=as.numeric(env), 
  n_breakpoints=1)
fit = stan(file="model.stan", data=data_list, 
  chains=3, iter=1000, 
  pars = c("b0", "slopes","b_0","b_s","cutpoints","pred","log_lik","effect"))
pars = rstan::extract(fit)
plot(apply(pars$pred,2,mean), logRS)
plot(env, apply(pars$effect,2,median))

# try to fit 2-breakpoint model in Stan
data_list = list(n_obs=N, logRS=logRS, 
  S=S, env=as.numeric(env), 
  n_breakpoints=2)
fit = stan(file="model.stan", data=data_list, 
  chains=3, iter=1000, 
  pars = c("b0", "slopes","b_0","b_s","cutpoints","pred","log_lik","effect"))
pars = rstan::extract(fit)
plot(apply(pars$pred,2,mean), logRS)
plot(env, apply(pars$effect,2,median))


