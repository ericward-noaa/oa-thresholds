###########################################################################
##   Fit simulated stock-recruit data with non-linear covariate effect   ##
###########################################################################
rm(list=ls()) 
library(rstan)
options(mc.cores=parallel::detectCores()) ## make stan() use multiple cores
rstan_options(auto_write=T) ## avoid unnecessary recompilation 

###################################################### generate random data
set.seed(13) 
N=50 ## data points
env=cumsum(rnorm(N,0,1)) # env covariate
env=as.numeric(scale(env))
S=round(exp(rnorm(N,5,1))) ## spawners 
alpha<-10 ## maximum productivity 
beta<-0.002 ## density dependence parameter
err<-rnorm(N,0,0.1) ## residual error (normal on log-scale)
# effect<-0.5*env 				## linear environmental effect 
effect<- -env^2 				 ## OR quadratic effect
# effect<- 2*env-env^2 			## OR linear plus quadratic effect
effect<-as.numeric(scale(effect,scale=2)) ## use 'scale' for effect size
effect<-effect-max(effect) ## set largest effect to zero
logRS<-log(alpha)+effect-beta*S+err
R<-exp(logRS)*S

######################################################## fit models in Stan
niter<-1e4 ## number of iterations 
start.time<-Sys.time() 

################################## simple linear model
data_list=list(n_obs=N,logRS=logRS,S=S,env=env)
params<-c("slope","b_0","b_s","pred","log_lik")#"b0","effect",
fit0=stan(file="stan/simple_linear_model.stan",data=data_list,chains=3,iter=niter,pars=params)
saveRDS(fit0,"results/fit0.rds")

################################## 1-breakpoint model 
data_list=list(n_obs=N,logRS=logRS,S=S,env=env,n_breakpoints=1)
params<-c("slopes","b_0","b_s","cutpoints","pred","log_lik")#"b0","effect",
fit1=stan(file="stan/model.stan",data=data_list,chains=3,iter=niter,pars=params)
saveRDS(fit1,"results/fit1.rds")

################################## 2-breakpoint model 
data_list=list(n_obs=N,logRS=logRS,S=S,env=env,n_breakpoints=2)
params<-c("slopes","b_0","b_s","cutpoints","pred","log_lik")#"b0","effect",
fit2=stan(file="stan/model.stan",data=data_list,chains=3,iter=niter,pars=params)
saveRDS(fit2,"results/fit2.rds")

###########################################################################
elapsed<-Sys.time()-start.time
print(elapsed)
###########################################################################
