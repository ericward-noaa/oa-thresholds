###########################################################################
##                                                                       ##
##                  PWS wild pink salmon data analysis                   ##
##                                                                       ##
###########################################################################
# rm(list=ls()) 
pkgs<-c("readxl","tidyverse","Hmisc","RColorBrewer", "nimble","mcmcplots","coda","bayesplot","hexbin","LearnBayes")
if(length(setdiff(pkgs,rownames(installed.packages())))>0) { install.packages(setdiff(pkgs,rownames(installed.packages())), dependencies=TRUE) }
lapply(pkgs,library,character.only=T)

###########################################################################
###################################################################### data
###########################################################################

########################################### PWS wild pink salmon abundances
pws_pink_data<-data.frame(read_excel("code_NIMBLE/data/Rich/2018_PWS_Wild_Pink_SR.xlsx"))
nY<-dim(pws_pink_data)[1]
fac<-1e-6 ## run and esc in millions of fish
pws_pink_data$Harvest<-round(pws_pink_data$Harvest*fac,6) 
pws_pink_data$Run.Size<-round(pws_pink_data$Run.Size*fac,6) 
pws_pink_data$Escapement<-round(pws_pink_data$Escapement*fac,6)
pws_pink_data$Spawners<-c(NA,NA,pws_pink_data$Escapement[-c(nY-1,nY)])
pws_pink_data$logRS<-log(pws_pink_data$RperS) ## log(recruits/spawner)

################################################# PWS hatchery release data
release_data<-data.frame(read_excel("code_NIMBLE/data/Rich/HR011558.xlsx"))
##-----------------------------------------------------------## filter data
releases_all<-dplyr::filter(release_data,Rearing.Code=="H", Species=="PINK",!is.na(Total.Released),!is.na(Year.Released), Location..Facility.or.Wild.Stock.!="TUTKA BAY")
##--------------------------------------------------------## select columns
releases<-dplyr::select(releases_all,Year.Brood,Stock,Stage,Weight,Length, Year.Released, Release.Site,Date.Last.Released,Total.Released)
##--------------------------------------------------## sum releases by year
pws_pink_releases<-aggregate(list(Total=releases$Total.Released), list(Year.Released=releases$Year.Released),sum)
pws_pink_releases$Total<-pws_pink_releases$Total*1e-6 #millions

###################################### NE Pacific salmon competitor biomass
competitor_bio<-data.frame(read_excel("code_NIMBLE/data/Ruggerone_Irvine_2018/Total-Salmon-Biomass.xlsx")) ## in thousand metric tons
competitor_bio<-dplyr::select(competitor_bio,Return.Year,Chum,Sockeye)
competitor_bio$competitors<-rowSums(competitor_bio[,-1])
competitor_bio<-dplyr::select(competitor_bio,Return.Year,competitors)
names(competitor_bio)<-c("Year","competitors")
competitor_bio$competitors<-competitor_bio$competitors*1e-3 ## million tons

########################################### regional chum/sockeye abundance
regions<-c("Kod","CI","PWS","SEAK","NBC","SBC","WA") ## without western AK
##------------------------------------------------------------------## chum
chum_count<-read_excel("code_NIMBLE/data/Ruggerone_Irvine_2018/Chum-Salmon-Total-Abundance.xlsx")
chum_count$regional<-rowSums(chum_count[,names(chum_count) %in% regions])
chum_count<-dplyr::select(chum_count,Year,regional) %>% data.frame() 
##---------------------------------------------------------------## sockeye
sock_count<-read_excel("code_NIMBLE/data/Ruggerone_Irvine_2018/Sockeye-Salmon-Total-Abundance.xlsx")
sock_count$regional<-rowSums(sock_count[,names(sock_count) %in% regions])
sock_count<-dplyr::select(sock_count,Year,regional) %>% data.frame()  
##-----------------------------------------------------------------## merge
regional_comp<-merge(chum_count,sock_count,by="Year")
regional_comp$comp_regional<-rowSums(regional_comp[,-1])
regional_comp<-dplyr::select(regional_comp,Year,comp_regional)

#################################################### NEW environmental data
env_wide<-data.frame(read.csv("code_NIMBLE/data/Mary/updated_goa_climate_data.csv"))

##################################################################### pCO2
## 1) Prince William Sound in spring or summer of BY+1 (outmigration)
## 2) Northern GoA east of Kodiak in the fall of BY+1 (outmigration)
## 3) Western GoA south of Aleutians in winter of BY+2 (return year)
##-------------------------------------------------------------## get data
pCO2_PWS<-read.csv("code_NIMBLE/data/Claudine/Mean_pCO2_by_year_month_PWS.csv")
pCO2_nGoA<-read.csv("code_NIMBLE/data/Claudine/Mean_pCO2_by_year_month_nGoA.csv")
pCO2_wGoA<-read.csv("code_NIMBLE/data/Claudine/Mean_pCO2_by_year_month_wGoA.csv")
names(pCO2_PWS)[1]<-names(pCO2_nGoA)[1]<-names(pCO2_wGoA)[1]<-"Year"
##------------------------------------------------------------------## PWS
pCO2_PWS$pCO2_PWS<-rowMeans(pCO2_PWS[,names(pCO2_PWS) %in% c("April","May","June")])
pCO2_PWS<-dplyr::select(pCO2_PWS,Year,pCO2_PWS) 
##-----------------------------------------------------------------## nGoA
pCO2_nGoA$pCO2_nGoA<-rowMeans(pCO2_nGoA[,names(pCO2_nGoA) %in% c("October","November","December")])
pCO2_nGoA<-dplyr::select(pCO2_nGoA,Year,pCO2_nGoA)
##-----------------------------------------------------------------## wGoA
pCO2_wGoA$pCO2_wGoA<-rowMeans(pCO2_wGoA[,names(pCO2_wGoA) %in% c("January","November","December")])
pCO2_wGoA<-dplyr::select(pCO2_wGoA,Year,pCO2_wGoA)

###########################################################################
###################################################### model data with lags
###########################################################################
data<-pws_pink_data
data$Year<-data$Run.Year ## express everything in run year
is.even<-function(x) { x %% 2 == 0 }
data$BL<-factor(is.even(data$Year),labels=c("odd","even"))
data<-dplyr::select(data,-Harvest,-Escapement,-RperS,-Brood.Line,-Run.Year)
nY<-dim(data)[1]
##-----------------------------------------------## previous years run size
## inter-brood-line competition
data$prevYearRunSize<-c(NA,data[-nY,]$Run.Size) 
##-----------------------## PWS pink releases in year prior to outmigration
## interbroodline competition 
pws_release_ts<-pws_pink_releases
names(pws_release_ts)<-paste0("PWS_released_prev_",names(pws_release_ts))
pws_release_ts$Year<-pws_release_ts$PWS_released_prev_Year.Released+2 
data<-merge(data,pws_release_ts[,-1],by="Year",all=T)
data$PWS_released_prev_Total[is.na(data$PWS_released_prev_Total)]<-0
## assumes no data means no hatchery releases
##-----------------------------------------## regional competitor abundance
same_year_region_comp<-regional_comp
names(same_year_region_comp)<-paste0("sameYr",names(same_year_region_comp))
data<-merge(data, same_year_region_comp,by.x="Year",by.y="sameYrYear",all=T)
##-----------------------------------------------------## ocean temperature
env_ts<-dplyr::select(env_wide,year,egoa.spr.sst)
env_ts$year<-env_ts$year+1 ## lag-1 (previous year) 
data<-merge(data,env_ts,by.x="Year",by.y="year",all=T)
##-------------------------------------------------------------------## pCO2
# pCO2_PWS$Year<-pCO2_PWS$Year+1 ## lag-1: previous year in PWS in summer
# data<-merge(data,pCO2_PWS,by="Year",all=T)
# pCO2_nGoA$Year<-pCO2_nGoA$Year+1 ## lag-1: previous year in nGoA in fall
# data<-merge(data,pCO2_nGoA,by="Year",all=T)
# pCO2_wGoA$Year<-pCO2_wGoA$Year ## lag-0: return year in wGoA in winter
# data<-merge(data,pCO2_wGoA,by="Year",all=T)
##-----------------------------------------------------## add regime shifts
## variance in Aleutian Low shifted in 1989
data$regime<-NA
for(i in 1:nY) { 
  if(data$Year[i]<=1989) { data$regime[i]<-"upto1988" }
  if(data$Year[i]>=1990) { data$regime[i]<-"since1989" }
}
data$regime<-factor(data$regime)

############################################################## select years
# years<-seq(1962,2015,1) ## years used in analysis
# data<-dplyr::filter(data,Year %in% years)
##-------------------------------------------------## or use complete cases
data<-data[complete.cases(data),]
years<-as.numeric(data$Year)

########################################## continuous variable correlations
test<-dplyr::select(data,-Year,-BL,-regime)
test<-apply(test,2,function(x) as.numeric(as.character(x)))
out<-rcorr(test)$r
out_mat<-matrix(out,nrow=dim(out)[1],ncol=dim(out)[2])
diag(out_mat)<-NA
colnames(out_mat)<-rownames(out_mat)<-colnames(out)
round(out_mat,4) ## still lots of covariation

###########################################################################
#################################################################### NIMBLE
###########################################################################
y<-data$logRS
x1<-data$Spawners # -mean(data$Spawners)
x2<-data$PWS_released_prev_Total # -mean(data$PWS_released_prev_Total)
x3<-data$egoa.spr.sst # -mean(data$egoa.spr.sst)
x4<-data$sameYrcomp_regional # -mean(data$sameYrcomp_regional)
x5<-as.numeric(data$BL) ## 1=odd, 2=even
x6<-as.numeric(data$regime) 
n<-dim(data)[1]

##============================================================## full model
code<-nimbleCode({
  sigma~dunif(0,10) ## prior on variance 
  beta0~dnorm(0,sd=10) ## prior on intercept
  for(i in 1:nbetas) {
    beta[i]~dnorm(0,sd=10) ## priors on slopes
  }  
  for(i in 1:n) {
    y[i]~dnorm(beta0 
               +beta[1]*x1[i]			## Spawners
               +beta[2]*x2[i]			## Releases
               +beta[3]*x2[i]^2		## Releases^2
               +beta[4]*x3[i]			## SST
               +beta[5]*x3[i]^2		## SST^2
               +beta[6]*x4[i]			## Competitors
               +beta[7]*x4[i]^2		## Competitors^2
               +beta[8]*x5[i]			## Brood line
               +beta[9]*x6[i]			## Regime
               +beta[10]*x6[i]*x3[i]	## Regime:SST
               +beta[11]*x6[i]*x1[i],	## Regime:Spawners
               sd=sigma) 
  }
})
##--------------------------------------------------## data/constants/inits
nbetas<-11
mydata<-list(y=y,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6)
constants<-list(n=n,nbetas=nbetas)
inits<-list(sigma=1,beta0=0,beta=rnorm(constants$nbetas))
##--------------------------------------## build model and get mcmc samples
model<-nimbleModel(code,constants=constants,data=mydata,inits=inits) 
mymcmc<-nimbleMCMC(model,niter=55000,nburnin=5000,thin=50,nchains=3, summary=TRUE,WAIC=TRUE) 
sample_chains<-mymcmc$samples
summary<-mymcmc$summary
WAIC_full<-mymcmc$WAIC
samples<-rbind(sample_chains[[1]],sample_chains[[2]],sample_chains[[3]])

##==================================================## WAIC model selection
model_selection<-"no"  
if(model_selection=="yes") {
  code<-nimbleCode({
    sigma~dunif(0,10) ## prior for variance 
    beta0~dnorm(0,sd=10) ## prior on intercept
    for(i in 1:nbetas) {
      beta[i]~dnorm(0,sd=10)	## priors on slopes
      zbeta[i]<-z[i]*beta[i]	## indicator variables for covar selection
    }  
    for(i in 1:n) {
      y[i]~dnorm(beta0 
                 +zbeta[1]*x1[i]			## Spawners
                 +zbeta[2]*x2[i]			## Releases
                 +zbeta[3]*x2[i]^2		## Releases^2
                 +zbeta[4]*x3[i]			## SST
                 +zbeta[5]*x3[i]^2		## SST^2
                 +zbeta[6]*x4[i]			## Competitors
                 +zbeta[7]*x4[i]^2		## Competitors^2
                 +zbeta[8]*x5[i]			## Brood line
                 +zbeta[9]*x6[i]			## Regime
                 +zbeta[10]*x6[i]*x3[i]	## Regime:SST
                 +zbeta[11]*x6[i]*x1[i],	## Regime:Spawners
                 sd=sigma) 
    }
  })
##----------------------------------------------------## all betas in model
  allbetas<-c("Spawners","Releases","Releases^2","SST","SST^2","Competitors", "Competitors^2","BroodLine", "Regime","Regime:SST","Regime:Spawners")
  nbetas<-length(allbetas)
##--------------------------------------------## list of indicator vectors
  nM<-1e6 ## number of randomly selected models
  ind_array<-array(dim=c(nM,nbetas))
  for(j in 1:nM){ ind_array[j,]<-sample(c(0,1),nbetas,replace=T) }
  ind_array<-ind_array[!duplicated(ind_array),]
  nmod<-dim(ind_array)[1] ## max 2048 models
##------------------------------------------## loop through indicator list
  start<-Sys.time()
  WAICs<-NA
  for(j in 1:nmod){
    print(paste0("model being fit = ",j))
    z<-ind_array[j,]
    mydata<-list(y=y,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,z=z)
    constants<-list(n=n,nbetas=nbetas)
    inits<-list(sigma=1,beta0=0,beta=rnorm(constants$nbetas))
    model<-nimbleModel(code,constants=constants,data=mydata,inits=inits) 
    mymcmc<-nimbleMCMC(model,niter=5000,nburnin=1000,thin=10,nchains=3, summary=TRUE,WAIC=TRUE) ## NOTE: based on relatively few samples 
    WAICs[j]<-mymcmc$WAIC
  }
  end<-Sys.time() 
  print(end-start) ## (~0.5min/run)
  ind_array<-data.frame(ind_array)
##----------------------------------------------## make model results table
  full_mod<-rep(1,nbetas)
  selected_mod<-rep(0,nbetas)
  selected_mod[keep_betas]<-1
  all_mods<-data.frame(rbind(selected_mod,full_mod,ind_array))
  names(all_mods)<-allbetas
  rownames(all_mods)<-c("selected","full",paste0("model",seq(nmod)))
  all_mods$WAIC<-c(WAIC_selected,WAIC_full,WAICs)
  all_mods$WAIC<-round(all_mods$WAIC,3)
  all_mods$deltaWAIC<-all_mods$WAIC-min(all_mods$WAIC)
##-----------------------------------------## save table with model results
  all_mods<-all_mods[order(all_mods$deltaWAIC),]
  write.csv(all_mods,"plots_nimble/model_selection_table.csv")
##------------------------------------------------------## end if statement
}

##====================================================## fit selected model
code<-nimbleCode({
  sigma~dunif(0,10) ## prior for variance 
  beta0~dnorm(0,sd=10) ## prior on intercept
  for(i in 1:nbetas) {
    beta[i]~dnorm(0,sd=10) ## priors on slopes
  }  
  for(i in 1:n) { 
    y[i]~dnorm(beta0 		## Intercept
               +beta[1]*x1[i] 			## Spawners
               +beta[2]*x2[i]			## Releases
               +beta[3]*x3[i] 			## SST
               +beta[4]*x4[i]			## Competitors
               +beta[5]*x4[i]^2 		## Competitors^2
               +beta[6]*x5[i] 			## Brood line
               +beta[7]*x6[i] 			## Regime
               +beta[8]*x6[i]*x3[i],	## Regime:SST
               sd=sigma) 				## Error term
  }  
})
##--------------------------------------------------## data/constants/inits
nbetas<-8
mydata<-list(y=y,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6)
constants<-list(n=n,nbetas=nbetas)
inits<-list(sigma=1,beta0=0,beta=rnorm(constants$nbetas))
##-----------------------------------------------------------## build model
model<-nimbleModel(code,constants=constants,data=mydata,inits=inits) 
##----------------------------------------------------------## mcmc samples
mymcmc<-nimbleMCMC(model,niter=55000,nburnin=5000,thin=50,nchains=3, summary=TRUE,WAIC=TRUE) 
sample_chains<-mymcmc$samples
summary<-mymcmc$summary
WAIC_selected<-mymcmc$WAIC
##---------------------------## combine chanins and re-name sample columns
samples<-rbind(sample_chains[[1]],sample_chains[[2]],sample_chains[[3]])
varnames<-c("Spawners","Releases","SST","Competitors","Competitors^2","BL", "Regime", "Interaction","Intercept","Sigma")
nvar<-length(varnames)
colnames(samples)<-varnames

##=====================================================## model diagnostics
pdf("code_NIMBLE/plots/traceplots_by_chain.pdf")
color_scheme_set("mix-blue-red")
mcmc_trace(sample_chains,size=0.1)
dev.off()

pdf("code_NIMBLE/plots/posterior_by_chain.pdf")
mcmc_hist_by_chain(sample_chains)
dev.off()

pdf("code_NIMBLE/plots/posteriors.pdf")
# mcmc_hist(samples) ## histograms
mcmc_dens(samples) ## density plots
dev.off()

pdf("code_NIMBLE/plots/posteriors_by_chain.pdf")
mcmc_dens_overlay(sample_chains)
dev.off()

pdf("code_NIMBLE/plots/posterior_intervals.pdf")
mcmc_intervals(samples[,1:nbetas])
dev.off()

pdf("code_NIMBLE/plots/corrrelations.pdf")
# pairs(samples,pch=16,cex=0.5)
pairs(samples[,c(1,2,3,4,6,7)],pch=16,cex=0.5)
# mcmc_pairs(sample_chains,diag_fun="dens",off_diag_fun="scatter", np_style=pairs_style_np(div_size=0.1,td_size=0.1)) ## or use bayesplot fn
dev.off()

pdf("code_NIMBLE/plots/corrrelation_releases_regime.pdf")
mcmc_scatter(samples,pars=c("Releases","Regime"),alpha=0.5,size=1)
dev.off()

pdf("code_NIMBLE/plots/corrrelation_hex_releases_regime.pdf",height=4,width=4)
col_fun<-colorRampPalette(brewer.pal(9,'Blues'))
releases<-samples[,colnames(samples)=="Releases"]
regime<-samples[,colnames(samples)=="Regime"]
hexbinplot(regime~releases,xlab="Regime",ylab="Releases",colramp=col_fun, colorkey=F)
# mcmc_hex(samples,pars=c("Releases","Regime")) ## or use bayesplot fn
dev.off()

###########################################################################
###########################################################################
###########################################################################