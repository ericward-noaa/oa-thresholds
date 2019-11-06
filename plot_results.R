###########################################################################
##                   Plot results for Stan model output                  ##
###########################################################################
# fit0<-readRDS("results/fit0.rds") 
# fit1<-readRDS("results/fit1.rds") 
# fit2<-readRDS("results/fit2.rds") 

################################################################# plot data
pdf("results/Data_simulated.pdf",height=5,width=10)
layout(matrix(1:8,ncol=4,nrow=2,byrow=T))
par(mar=c(5,5,1,1),pch=16)
plot(env,type="o")
plot(env,effect)
lines(loess.smooth(env,effect,degree=2))
abline(v=env[effect==max(effect)],lty=3)
plot(env,log(alpha)+effect+err) ## density-independent effects
plot(env,logRS);abline(h=0,lty=3) ##logRS with replacement line
plot(R,type="l",ylim=c(0,max(R)));abline(h=min(R),lty=3)
plot(S,type="l",ylim=c(0,max(S)));abline(h=min(S),lty=3)
plot(S,logRS);lines(S,predict(lm(logRS~S)))
plot(S,R,xlim=c(0,max(S)),ylim=c(0,max(R)))
##----------------------------------## fit stock-recruitment curve
mod<-lm(logRS~S)
inter<-summary(mod)$coefficients[1,1] ## intercept
slope<-summary(mod)$coefficients[2,1] ## slope
a.raw<-exp(inter) ## productivity parameter (alpha)
b.raw<--slope ## density-dependence parameter (beta)
sig<-summary(mod)$sigma ## residual sigma to correct for log-normal error
a<-a.raw+0.5*sig^2 ## corrected alpha value 
b<-b.raw*(a/a.raw) ## corrected beta value 
spns<-seq(0,max(S),length=1e3)
Ricker<-function(S,a,b) { a*S*exp(-b*S) } 
lines(spns,Ricker(S=spns,a=a,b=b),lwd=2) ## Ricker fit
abline(1,1,lty=3) ## replacement line
abline(0,a,lty=3,col=3) ## estimated productivity
##----------------------------------## S_msy (Hilborn & Walters 1992)
S_msy<-(log(a)/b)*(0.5-0.07*log(a)) 
abline(v=S_msy,lty=2,lwd=2)
##----------------------------------## or using surplus production
recs<-a*spns*exp(-b*spns)
surplus<-recs-spns
index<-which(surplus==max(surplus)) 
E_msy<-spns[index]
abline(v=E_msy,lty=1,lwd=0.5,col=2)
dev.off()

############################################## model comparison using looic
## leave-one-out cross validation
##-------------------------------## 0 breakpoints
loo_stats0<-loo(fit0,pars="log_lik")
loo0<-loo_stats0$estimates[1,1]
##-------------------------------## 1 breakpoints
loo_stats1<-loo(fit1,pars="log_lik")
loo1<-loo_stats1$estimates[1,1]
##-------------------------------## 2 breakpoints
loo_stats2<-loo(fit2,pars="log_lik")
loo2<-loo_stats2$estimates[1,1]
##-------------------------------## looic table
looic<-data.frame(cbind(breaks=c(0,1,2),looic=c(loo0,loo1,loo2)))
looic$delta_looic=looic$looic-min(looic$looic,na.rm=T)
looic$weights<-round(loo::loo_model_weights(list(loo_stats0,loo_stats1, loo_stats2),method="stacking"),6)
looic
write.csv(looic,file="results/looic-table.csv")

########################################################### get parameters

#################################### simple linear (0-breakpoint) model
pars0=rstan::extract(fit0)
logRS_pred0<-apply(pars0$pred,2,median) ## median for logRS
R_pred0<-exp(logRS_pred0)*S
alpha0<-median(exp(pars0$b_0)) ## median for alpha
beta1<-median(pars0$b_s) ## median for beta
slope0<-median(pars0$slope)

#################################### fit 1-breakpoint model 
pars1=rstan::extract(fit1)
logRS_pred1<-apply(pars1$pred,2,median) ## median for logRS
R_pred1<-exp(logRS_pred1)*S
alpha1<-median(exp(pars1$b_0)) ## median for alpha
beta1<-median(pars1$b_s) ## median for beta
cuts1<-median(pars1$cutpoints)
slopes1<-apply(pars1$slopes,2,function(x) median(x))

#################################### fit 2-breakpoint model 
pars2=rstan::extract(fit2)
logRS_pred2<-apply(pars2$pred,2,median) ## median for logRS
R_pred2<-exp(logRS_pred2)*S
alpha2<-median(exp(pars2$b_0)) ## median for alpha
beta2<-median(pars2$b_s,prob=probs) ## median for beta
cuts2<-apply(pars2$cutpoints,2,function(x) median(x))
slopes2<-apply(pars2$slopes,2,function(x) median(x))

############################################################ results plots

######################################## priors and posteriors of cutpoints
pdf("results/Prior_posterior_breakpoints.pdf",height=3,width=9)
layout(matrix(1:3,ncol=3,nrow=1,byrow=T))
par(mar=c(5,5,1,1))
xlim<-c(-2,2)
##-------------------------------## 1 breakpoints
plot(density(rnorm(1e5,0,1)),xlab="",main="",lwd=2,ylim=c(0,5),xlim=xlim)
lines(density(pars1$cutpoints[,1]),lwd=2,col=2)
##-------------------------------## 2 breakpoints
plot(density(rnorm(1e5,-1,1)),xlab="",main="",lwd=2,ylim=c(0,5),xlim=xlim)
lines(density(pars2$cutpoints[,1]),lwd=2,col=2)
plot(density(rnorm(1e5,1,1)),xlab="",main="",lwd=2,ylim=c(0,5),xlim=xlim)
lines(density(pars2$cutpoints[,2]),lwd=2,col=2)
dev.off()

############################################ observed vs predicted log(R/S)
pdf("results/Model_comparison_obs_vs_pred_logRS.pdf",height=3,width=9)
layout(matrix(1:3,ncol=3,nrow=1,byrow=T))
par(mar=c(5,5,1,1))
##-------------------------------## 0 breakpoints
plot(logRS_pred0,logRS)
lmfit<-lm(logRS~logRS_pred0)
lines(logRS_pred0,predict(lmfit))
r2<-signif(summary(lmfit)$adj.r.squared,2)
mtext(paste0(" r2=",r2),side=3,line=-1,cex=0.6,adj=0)
##-------------------------------## 1 breakpoints
plot(logRS_pred1,logRS)
lmfit<-lm(logRS~logRS_pred1)
lines(logRS_pred1,predict(lmfit))
r2<-signif(summary(lmfit)$adj.r.squared,2)
mtext(paste0(" r2=",r2),side=3,line=-1,cex=0.6,adj=0)
##-------------------------------## 2 breakpoints
plot(logRS_pred2,logRS)
lmfit<-lm(logRS~logRS_pred2)
lines(logRS_pred2,predict(lmfit))
r2<-signif(summary(lmfit)$adj.r.squared,2)
mtext(paste0(" r2=",r2),side=3,line=-1,cex=0.6,adj=0)
dev.off()

############################################ predicted environmental effect
pdf("results/Model_comparison_env_effect.pdf",height=3,width=9)
layout(matrix(1:3,ncol=3,nrow=1,byrow=T))
par(mar=c(4.5,4.5,1.5,0.5),oma=c(0.5,0.5,0.5,0.5),cex.lab=1.5,pch=16)
nsamp<-dim(pars2$pred)[1]
probs<-c(0,0.05,0.5,0.95,1)
cols<-c("goldenrod1","chocolate2","firebrick2","midnightblue")
nsub<-500 ## number of samples to use from posterior
subsamp<-sample(seq(nsamp),nsub,replace=F) ## random subsample
tcols<-alpha(cols,0.1);lwd<-1
xlim<-c(min(env),max(env))
ymax<-max(logRS);ylim<-c(min(effect),ymax)
##=========================================================## 0 breakpoints
plot(NA,NA,xlab="Covariate value",ylab="Effect",xlim=xlim,ylim=ylim)
abline(v=env[effect==max(effect)],lty=3)
abline(h=0,lty=3)
points(env,effect,cex=1,pch=16,col=1) ## env effect
points(env,logRS,cex=1,pch=16,col="darkgray") ## logRS
##----------------------------------------## slope
for(i in subsamp) {
x0<-min(env)
y0<-pars0$slope[i]*min(env)
x1<-max(env)
y1<-pars0$slope[i]*max(env)
segments(x0,y0,x1,y1,lwd=lwd,col=tcols[1])
}
##----------------------------------------## medians
x0<-min(env)
y0<-median(pars0$slope)*min(env)
x1<-max(env)
y1<-median(pars0$slope)*max(env)
segments(x0,y0,x1,y1,lwd=1,col=1)
##=========================================================## 1 breakpoints
plot(NA,NA,xlab="Covariate value",ylab="Effect",xlim=xlim,ylim=ylim)
abline(v=env[effect==max(effect)],lty=3)
abline(h=0,lty=3)
points(env,effect,cex=1,pch=16,col=1) ## env effect
points(env,logRS,cex=1,pch=16,col="darkgray") ## logRS
##----------------------------------------## cutspoints
cuts_1<-quantile(pars1$cutpoints,prob=probs)
points(cuts_1[3],ymax,cex=1.5,col=cols[2],xpd=T)
segments(cuts_1[1],ymax,cuts_1[5],ymax,lwd=0.5,col=cols[2],xpd=T)
segments(cuts_1[2],ymax,cuts_1[4],ymax,lwd=2,col=cols[2],xpd=T)
##----------------------------------------## slopes
for(i in subsamp) {
x0<-min(env)
y0<-pars1$slopes[i,1]*min(env)
x1<-pars1$cutpoints[i,1]
y1<-pars1$slopes[i,1]*pars1$cutpoints[i,1]
segments(x0,y0,x1,y1,lwd=lwd,col=tcols[2])
x2<-max(env)
y2<-pars1$slopes[i,1]*pars1$cutpoints[i,1]+pars1$slopes[i,2]*max(env)
segments(x1,y1,x2,y2,lwd=lwd,col=tcols[2])
}
##----------------------------------------## medians
x0<-min(env)
y0<-median(pars1$slopes[,1])*min(env)
x1<-median(pars1$cutpoints[,1])
y1<-median(pars1$slopes[,1])*median(pars1$cutpoints[,1])
segments(x0,y0,x1,y1,lwd=2,col=1)
x2<-max(env)
y2<-median(pars1$slopes[,1])*median(pars1$cutpoints[,1])+median(pars1$slopes[i,2])*max(env)
segments(x1,y1,x2,y2,lwd=2,col=1)
##=========================================================## 2 breakpoints
plot(NA,NA,xlab="Covariate value",ylab="Effect",xlim=xlim,ylim=ylim)
abline(v=env[effect==max(effect)],lty=3)
abline(h=0,lty=3)
points(env,effect,cex=1,pch=16,col=1) ## env effect
points(env,logRS,cex=1,pch=16,col="darkgray") ## logRS
##----------------------------------------## cutpoints
cuts_2<-apply(pars2$cutpoints,2,function(x) quantile(x,prob=probs))
points(cuts_2[3,1],ymax,cex=1.5,col=cols[3],xpd=T)
segments(cuts_2[1,1],ymax,cuts_2[5,1],ymax,lwd=0.5,col=cols[3],xpd=T)
segments(cuts_2[2,1],ymax,cuts_2[4,1],ymax,lwd=2,col=cols[3],xpd=T)
ymax<-ymax+0.1;points(cuts_2[3,2],ymax,cex=1.5,col=cols[3],xpd=T)
segments(cuts_2[1,2],ymax,cuts_2[5,2],ymax,lwd=0.5,col=cols[3],xpd=T)
segments(cuts_2[2,2],ymax,cuts_2[4,2],ymax,lwd=2,col=cols[3],xpd=T)
##----------------------------------------## slopes
for(i in subsamp) {
x0<-min(env)
y0<-pars2$slopes[i,1]*min(env)
x1<-pars2$cutpoints[i,1]
y1<-pars2$slopes[i,1]*pars2$cutpoints[i,1]
segments(x0,y0,x1,y1,lwd=lwd,col=tcols[3])
x2<-pars2$cutpoints[i,2]
y2<-pars2$slopes[i,1]*pars2$cutpoints[i,1]
segments(x1,y1,x2,y2,lwd=lwd,col=tcols[3])
x3<-max(env)
y3<-pars2$slopes[i,1]*pars2$cutpoints[i,1]+pars2$slopes[i,2]*max(env)
segments(x2,y2,x3,y3,lwd=lwd,col=tcols[3])
}
##----------------------------------------## medians
x0<-min(env)
y0<-median(pars2$slopes[,1])*min(env)
x1<-median(pars2$cutpoints[,1])
y1<-median(pars2$slopes[,1])*median(pars2$cutpoints[,1])
segments(x0,y0,x1,y1,lwd=2,col=1)
x2<-median(pars2$cutpoints[,2])
y2<-median(pars2$slopes[,1])*median(pars2$cutpoints[,1])
segments(x1,y1,x2,y2,lwd=2,col=1)
x3<-max(env)
y3<-median(pars2$slopes[,1])*median(pars2$cutpoints[,1])+median(pars2$slopes[,2])*max(env)
segments(x2,y2,x3,y3,lwd=2,col=1)
dev.off()

####################################### median and CIs of Ricker parameters
pdf("results/Model_comparison_alpha_beta.pdf",height=3,width=6)
layout(matrix(1:2,ncol=2,nrow=2,byrow=T))
par(mar=c(5,5,1,1),pch=16)
probs<-c(0.05,0.5,0.95)
max_breaks<-2
##====================================## log alphas
plot(NA,NA,xlab="# of breakpoints",xaxt="n",ylab="log(alpha)", xlim=c(-0.5,max_breaks+0.5),ylim=c(0,6))
axis(1,labels=T,at=seq(0,max_breaks,1))
abline(h=log(alpha),col=2)
##------------------------------------## 0 breakpoints
alpha0<-quantile(pars0$b_0,prob=probs)
points(0,alpha0[2])
segments(0,alpha0[1],0,alpha0[3])
##------------------------------------## 1 breakpoints
alpha1<-quantile(pars1$b_0,prob=probs)
points(1,alpha1[2])
segments(1,alpha1[1],1,alpha1[3])
##------------------------------------## 2 breakpoints
alpha2<-quantile(pars2$b_0,prob=probs)
points(2,alpha2[2])
segments(2,alpha2[1],2,alpha2[3])
##====================================## betas
plot(NA,NA,xlab="# of breakpoints",xaxt="n",ylab="beta", xlim=c(-0.5,max_breaks+0.5),ylim=c(1.5*min(beta1,beta2),0))
axis(1,labels=T,at=seq(0,max_breaks,1))
abline(h=-beta,col=2)
##------------------------------------## 0 breakpoints
beta0<-quantile(pars0$b_s,prob=probs) 
points(0,beta0[2])
segments(0,beta0[1],0,beta0[3])
##------------------------------------## 1 breakpoints
beta1<-quantile(pars1$b_s,prob=probs) 
points(1,beta1[2])
segments(1,beta1[1],1,beta1[3])
##------------------------------------## 2 breakpoints
beta2<-quantile(pars2$b_s,prob=probs)
points(2,beta2[2])
segments(2,beta2[1],2,beta2[3])
dev.off()

############################################## parameter correlation matrix
# library(PerformanceAnalytics)
# dev.new(width=6,height=6);par(pch=16)
# plotdat0<-cbind(b_0=pars0$b_0,b_s=pars0$b_s,b0=pars0$b0,slope=pars0$slope)
# chart.Correlation(plotdat0,histogram=T)
# quartz.save("results/Parameter_correlations_breaks=0.png",type="png")
# dev.off()
# 
# dev.new(width=6,height=6);par(pch=16)
# plotdat1<-cbind(b_0=pars1$b_0,b_s=pars1$b_s,slope1=pars1$slopes[,1], slope2=pars1$slopes[,2],cutpoint=pars1$cutpoints[,1])
# chart.Correlation(plotdat1,histogram=T)
# quartz.save("results/Parameter_correlations_breaks=1.png",type="png")
# dev.off()
# 
# dev.new(width=6,height=6);par(pch=16)
# plotdat2<-cbind(b_0=pars2$b_0,b_s=pars2$b_s,b0=pars2$b0, slope1=pars2$slopes[,1],slope2=pars1$slopes[,2], cutpoint1=pars2$cutpoints[,1],cutpoint2=pars2$cutpoints[,2])
# chart.Correlation(plotdat2,histogram=T)
# quartz.save("results/Parameter_correlations_breaks=2.png",type="png")
# dev.off()

###########################################################################