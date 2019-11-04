# This model code is used for analyses in manuscript Yang W, Kehm RD, and Terry MB 2019 (Survival model methods for analyses of cancer incidence trends in young adults)
# CAUTION: for use in other research, please read the code carefully and modify it to suit your own need.

# Functions for cohort specific survival analysis
# including 1) different survivial models
# 2) analyzing cohort trend using regression
# 3) make predictions for later cohorts

fn_cntNA=function(x) sum(!is.na(x))
# Function to plot model fits and save in a pdf file
fn_plotFits=function(best.fit,dat,fig.title='Survival model-fit',add.axis=T,add.legend=T,add.ylab=T,add.xlab=T){ # do not save as pdf so we can combine figures
  cc=palette(brewer.pal(n=9,name='OrRd')); colfunc=colorRampPalette(cc);cols=colfunc(ncol(best.fit)+2)[-c(1,2)]
  cc=palette(brewer.pal(n=9,name='OrRd')); colfunc=colorRampPalette(cc);cols=colfunc(ncol(best.fit)+2)[-c(1,2)]
  ymax=max(c(best.fit,unlist(dat[,-1])))*1.1
  
  matplot(best.fit,ylim=c(0,ymax),xaxt='n',xlab='',ylab='',col=cols,type='l',lty=1)
  matpoints(dat[,-1],col=cols,pch=1:ncol(best.fit))
  mtext(fig.title,side=3,line=-1.1,adj=.02,outer = F,cex=.8)
  
  if (add.axis==T) axis(1,at=1:nrow(dat),labels = dat$Age)
  if (add.ylab==T) mtext('Rate per 100,000',side=2,outer = F,line=.95,cex=.75,xpd=T)
  if (add.xlab==T) mtext('Age',side=1,outer = F,line=.95,cex=.75,xpd=T)
  if(add.legend==T) legend(1,ymax*.92,legend = colnames(dat)[-1],col=cols,lty=1,pch=1:ncol(best.fit),cex=.8,bty='n')
  if(add.note==T) mtext(fig.note,side=3,outer = T,adj = .02,line=.5,cex=.75)
}
fn_plotTrend=function(res.surv,fitted.cohort.trend,est.cohort.eff,fig.title='Estimated cohort trend',add.axis=T,add.legend=T,add.ylab=T,add.xlab=T,printRR=T){
  ymin=0; ymax=max(res.surv$hazard.rate*1e6)*1.57
  plot(res.surv$cohort.start,res.surv$hazard.rate*1e6,ylim=c(ymin,ymax), xaxt='n',xlab='',ylab='',pch=20)  # ifelse(surv.model=="fn_surv.reg.acc.fixed.hr0",'Change in hazard rate','Mean hazard rate')
  lines(res.surv$cohort.start,predict(fitted.cohort.trend)*1e6)
  mtext(fig.title,side=3,line=-1.1,adj=.02,outer = F,cex=.8)
  if (add.axis==T) axis(1,at=res.surv$cohort.start,labels = res.surv$cohort.start)
  if (add.ylab==T) mtext(expression(paste('Mean hazard rate (x10'^-6,')')),side=2,outer = F,line=.95,cex=.75,xpd=T)
  if (add.xlab==T) mtext('Cohort',side=1,outer = F,line=.95,cex=.75,xpd=T)
  if (printRR==T){
    mtext(paste0('Detected Cohort Effect: P=',format(est.cohort.eff$p.value,digits = 2)),col=ifelse(est.cohort.eff$p.value<.05,'red','black'),side=1,line=-2.3,adj=.05,outer = F,cex=.75)
    mtext(bquote(Delta*'RR ='~.(round(est.cohort.eff$deltaRR,2))*'% (95% CI: '~.(round(est.cohort.eff$deltaRR.lwr,2))*'%, '~.(round(est.cohort.eff$deltaRR.upr,2))*'%)'),col=ifelse(est.cohort.eff$p.value<.05,'red','black'),side=1,line=-1.2,adj=.05,outer = F,cex=.75)
  }
}

##########################################################################################
# FUNCTION FOR SURVIVIAL MODELING
##########################################################################################
## The 'constant' survival model: log(surv) ~ -cumulative hazard = -h*(t-t0) [assume: before t0, h=0]
fn_surv.reg.const=function(dat,age.min,age.max,N=1e5,plotit=F,return.bestfit=F){
  # input data structure - dat
  # column 1: age; columns 2 to n: cohort-specific incidence rate per 100,000 person-year
  # column names: age, and cohort birth years
  # survival model: log(surv) ~ -cumulative hazard = -h*(a-a0) [assume: before a0, h=0]
  
  # find a0 for the cancer site (same for all cohorts)
  for(a0 in 0:age.min){ # set the upper bound to age.min to make sure the intercept is non-negative
    res.surv=NULL; fits=NULL;
    
    for(ic in colnames(dat)[-1]){ # loop through cohorts
      tmp=dat[,c('Age',ic),with=F]
      tmp$Age=pmax(0,tmp$Age-a0) # need to adjust the original!
      tmp$Surv=1-tmp[,2]/N # survival rate = 1 - cancer rate
      fit=lm(log(Surv)~0+Age,data=tmp) # fit to the survival model
      
      y.fitted=pmax(0,(1-exp(fit$fitted.values))*N) # should >=0
      fits=cbind(fits,y.fitted) # save model fits for that cohort
      
      adj.r.squared=summary(fit)$adj.r.squared # get the adjusted r-squared
      f=summary(fit)$fstatistic
      p=pf(f[1],f[2],f[3],lower.tail=F) # get the P value
      # compute AIC, BIC 
      aic=AIC(fit); bic=BIC(fit)
      res.surv=rbind(res.surv,c(as.numeric(unlist(str_split(ic,'-'))),-fit$coefficients,summary(fit)$coefficients[,'Std. Error']*sqrt(nrow(tmp)),a0,adj.r.squared,p,aic,bic))
    }
    res.surv=data.frame(res.surv); colnames(res.surv)=c('cohort.start','cohort.end','hazard.rate','sd.hr','age0','adj.r.squared','p.value','aic','bic')
    if(a0==0) { # record the best fit
      best.surv=res.surv; best.fit=fits;
    } else {
      # select the best a0 based on AIC/BIC
      if (mean(res.surv$bic,na.rm = T)<mean(best.surv$bic,na.rm = T)) {best.surv=res.surv; best.fit=fits}
    }
  } # a0
  if(plotit==T){ # plot the fits
    fn_plotFits(best.fit,dat)
  }
  # compute the combined metric for comparison
  best.surv=rbind(best.surv,c(best.surv$cohort.start[1],tail(best.surv$cohort.end,1),colMeans(best.surv[,3:(ncol(best.surv)-3)]),colSums(best.surv[,ncol(best.surv)-2:0])))
  if(return.bestfit==F){ 
    return(best.surv)
  } else {
    return(list(best.surv=best.surv,best.fit=best.fit))
  }
  
}

# The 'acceleration' survival model: log(surv) ~ -cumulative hazard = -(h0*(a-a0)+1/2*delta*(a-a0)^2) [assume: before t0, h=0, and afterwards linear increase]
fn_surv.reg.acc=function(dat,age.min,age.max,N=1e5,plotit=F,return.bestfit=F){
  # input data structure - dat
  # column 1: age; columns 2 to n: cohort-specific incidence rate per 100,000 person-year
  # column names: age, and cohort birth years
  # model: log(surv) ~ -cumulative hazard = -(h0*(a-a0)+1/2*delta*(a-a0)^2) [assume: before t0, h=0, and afterwards linear increase]
  
  # find a0 for the cancer site (same for all cohorts)
  for(a0 in 0:age.min){ # set the upper bound to age.min to make sure the intercept is non-negative
    res.surv=NULL; fits=NULL;

    for(ic in colnames(dat)[-1]){ # loop through cohorts
      tmp=dat[,c('Age',ic),with=F]
      tmp$Age=pmax(0,tmp$Age-a0) # need to adjust the original!
      tmp$Age2=tmp$Age^2
      tmp$Surv=1-tmp[,2]/N
      fit=lm(log(Surv)~0+Age+Age2,data=tmp)
      if(-fit$coefficients['Age']<0) {
        fit=lm(log(Surv)~0+Age2,data=tmp); hr0=0;
      } else {
        hr0='nonneg'
      }
      # save model fits
      y.fitted=pmax(0,(1-exp(fit$fitted.values))*N) # should >=0
      fits=cbind(fits,y.fitted)
      
      adj.r.squared=summary(fit)$adj.r.squared # get the adjusted r-squared
      f=summary(fit)$fstatistic
      p=pf(f[1],f[2],f[3],lower.tail=F) # get the P value
      # compute AIC, BIC 
      aic=AIC(fit); bic=BIC(fit)
      if (hr0=='nonneg'){
        hrs=-fit$coefficients; hr.sds=summary(fit)$coefficients[,'Std. Error']*sqrt(nrow(tmp)); 
      } else {
        hrs=c(0,-fit$coefficients); hr.sds=c(0,summary(fit)$coefficients[,'Std. Error'])*sqrt(nrow(tmp)); # hazard.rate.mean=sum(hrs*c(1,max(tmp$Age)))
      }
      hazard.rate.mean.from.a0=sum(hrs*c(1,max(tmp$Age))); # v.ave=v0+1/2*a*t from a0 to a.max
      # hr1=hrs[1]+2*hrs[2]*(age.min-a0); hr2=hrs[1]+2*hrs[2]*(age.max-a0); => hr.mean=1/2*(hr1+hr2)
      hazard.rate.mean.from.amin=hrs[1]+hrs[2]*(age.max+age.min-2*a0)
      res.surv=rbind(res.surv,c(as.numeric(unlist(str_split(ic,'-'))),hrs,hr.sds,hazard.rate.mean.from.a0,hazard.rate.mean.from.amin,a0,adj.r.squared,p,aic,bic))
    }
    
    res.surv=data.frame(res.surv); colnames(res.surv)=c('cohort.start','cohort.end','hazard.rate0','half.hr.accel.rate','sd.hr0','sd.accel.hr','hazard.rate.mean.from.a0','hazard.rate.mean.from.amin','age0','adj.r.squared','p.value','aic','bic')
    if(a0==0) { # record the best fit
      best.surv=res.surv; best.fit=fits
    } else {
      # select the best a0 based on AIC/BIC
      if (mean(res.surv$bic,na.rm = T)<mean(best.surv$bic,na.rm = T)) {best.surv=res.surv; best.fit=fits}
    }
  } # a0

  if(plotit==T){ # plot the fits
    fn_plotFits(best.fit,dat)
  }
  # compute the combined metric for comparison
  best.surv=rbind(best.surv,c(best.surv$cohort.start[1],tail(best.surv$cohort.end,1),colMeans(best.surv[,3:(ncol(best.surv)-3)]),colSums(best.surv[,ncol(best.surv)-2:0])))
  if(return.bestfit==F){ 
    return(best.surv)
  } else {
    return(list(best.surv=best.surv,best.fit=best.fit))
  }
}

## The "2-steps" survival model: log(surv) ~ -cumulative hazard = -(h1*(a1-a0)+h2*(a2-a1)) # 2 steps
fn_surv.reg.2steps=function(dat,age.min,age.max,N=1e5,plotit=F,return.bestfit=F){
  # input data structure - dat
  # column 1: age; columns 2 to n: cohort-specific incidence rate per 100,000 person-year
  # column names: age, and cohort birth years
  # model: log(surv) ~ -cumulative hazard = -(h1*(a1-a0)+h2*(a2-a1)) # 2 steps
  
  # find a0 for the cancer site (same for all cohorts)
  age.mid=mean(c(age.min,age.max))
  for(a0 in 0:age.min){ # set the upper bound to age.min to make sure the intercept is non-negative
    res.surv=NULL; fits=NULL;

    for(ic in colnames(dat)[-1]){
      tmp=dat[,c('Age',ic),with=F]
      
      tmp$Age1=pmax(0,tmp$Age-a0) # need to adjust the original
      tmp$Age2=pmax(0,tmp$Age-age.mid) # those >age.mid
      tmp$Age1[which(tmp$Age2!=0)]=age.mid-a0   # those <=age.mid
      tmp$Surv=1-tmp[,2]/N
      fit=lm(log(Surv)~0+Age1+Age2,data=tmp)
      # save model fits
      y.fitted=pmax(0,(1-exp(fit$fitted.values))*N) # should >=0
      fits=cbind(fits,y.fitted)
      
      adj.r.squared=summary(fit)$adj.r.squared # get the adjusted r-squared
      f=summary(fit)$fstatistic
      p=pf(f[1],f[2],f[3],lower.tail=F) # get the p value
      # compute AIC, BIC 
      aic=AIC(fit); bic=BIC(fit)
      res.surv=rbind(res.surv,c(as.numeric(unlist(str_split(ic,'-'))),-fit$coefficients,summary(fit)$coefficients[,'Std. Error']*sqrt(nrow(tmp)),a0,age.mid,adj.r.squared,p,aic,bic))
    }
    res.surv=data.frame(res.surv); colnames(res.surv)=c('cohort.start','cohort.end','hazard.rate0','hazard.rate1','sd.hr0','sd.hr1','age0','age1','adj.r.squared','p.value','aic','bic')
    if(a0==0) { # record the best fit
      best.surv=res.surv; best.fit=fits
    } else {
      # select the best a0 based on AIC/BIC
      if (mean(res.surv$bic,na.rm = T)<mean(best.surv$bic,na.rm = T)) {best.surv=res.surv; best.fit=fits}
    }
  } # a0
  if(plotit==T){ # plot the fits
    fn_plotFits(best.fit,dat)
  }
  # compute the combined metric for comparison
  best.surv=rbind(best.surv,c(best.surv$cohort.start[1],tail(best.surv$cohort.end,1),colMeans(best.surv[,3:(ncol(best.surv)-3)]),colSums(best.surv[,ncol(best.surv)-2:0])))
  if(return.bestfit==F){ 
    return(best.surv)
  } else {
    return(list(best.surv=best.surv,best.fit=best.fit))
  }
}

##########################################################################################
# FUNCTION FOR ESTIMATING COHORT EFFECT
##########################################################################################
fn_get.cohort.trend=function(res.surv,dat,surv.model="fn_surv.reg.acc",return.res=T){
  # input: res.surv=results from the survival model
  #        dat=original incidence data
  
  dat=dat[Age<=age.max & Age>=age.min]
  idx=which(unlist(lapply(dat,fn_cntNA))>=(age.max-age.min))
  dat=dat[,idx,with=F]
  if(sum(dat[,-1])<1e-2) return(NULL)
  ir.mean=mean(as.numeric(unlist(dat[,-1]))) # mean risk across ages and cohorts
  ir.sd=sd(as.numeric(unlist(dat[,-1])))
  
  if(res.surv$cohort.start[1]==res.surv$cohort.start[nrow(res.surv)]) res.surv=res.surv[-nrow(res.surv),] # if the last row is the mean across cohorts, delete it
  fitted.cohort.trend=lm(hazard.rate~cohort.start,data = res.surv) # linear regression against the cohort birth year
  beta=unname(fitted.cohort.trend$coefficients['cohort.start']) # get the slope: it represents the changes in risk per 5-year birth cohort. 
  p=summary(fitted.cohort.trend)$coefficients['cohort.start',4]  # get significant level of beta
  sign.beta=sign(beta)
  
  sd.beta=summary(fitted.cohort.trend)$coefficients['cohort.start','Std. Error']*sqrt(nrow(res.surv)); # compute the standard devidation of beta
  beta.lwr=unname(confint(fitted.cohort.trend)['cohort.start',][1]); beta.upr=unname(confint(fitted.cohort.trend)['cohort.start',][2]); # get the confidence intervals
  # compute relative changes in risk (survival rate) per cohort compared to the mean risk
  ref.risk=ir.mean # use the mean incidence rate as reference 
  deltaRR=sign(beta)*(1-exp(-abs(beta)*(age.max-age.min+1)))*N/ref.risk*100; # in %
  deltaRR.lwr=sign(beta.lwr)*(1-exp(-abs(beta.lwr)*(age.max-age.min+1)))*N/ref.risk*100; 
  deltaRR.upr=sign(beta.upr)*(1-exp(-abs(beta.upr)*(age.max-age.min+1)))*N/ref.risk*100;
  
  est.cohort.eff=data.frame(cohort.start.min=min(res.surv$cohort.start),cohort.start.max=max(res.surv$cohort.start),
                            ir.mean=ir.mean,ir.sd=ir.sd,beta=beta,sd.beta=sd.beta,beta.lwr=beta.lwr,beta.upr=beta.lwr,
                            deltaRR=deltaRR,deltaRR.lwr=deltaRR.lwr,deltaRR.upr=deltaRR.upr,
                            sign=sign.beta,p.value=p)
  if(return.res==T)
    return(list(fitted.cohort.trend=fitted.cohort.trend,est.cohort.eff=est.cohort.eff))
}

# do both survival modeling and fitting the cohort trend together
fn_fit.trend=function(dat,surv.model="fn_surv.reg.acc",N=1e5,fig.title1='Survival model-fit',fig.title2='Estimated cohort trend',
                      add.axis=add.axis,add.legend=add.legend,add.ylab=add.ylab,add.xlab=add.xlab,printRR=T,return.res=T){
  # input data structure - dat
  # column 1: age; columns 2 to n: cohort-specific incidence rate per 100,000 person-year
  # column names: age, and cohort birth years
  
  fn_surv.reg=get(surv.model)
  dat=dat[Age<=age.max & Age>=age.min]
  idx=which(unlist(lapply(dat,fn_cntNA))>=(age.max-age.min))
  dat=dat[,idx,with=F]
  if(sum(dat[,-1])<1e-2) return(NULL)
  ir.mean=mean(as.numeric(unlist(dat[,-1]))) # mean risk across ages and cohorts
  ir.sd=sd(as.numeric(unlist(dat[,-1])))
  
  # fit to the survival model
  res=fn_surv.reg(dat,age.min,age.max,plotit=F,return.bestfit = T) 
  res.surv=res$best.surv; surv.fits=res$best.fit; age0=res.surv$age0[1]
  
  # reshape model output from the survival model before linear regresssion 
  if(surv.model=="fn_surv.reg.acc"){ 
    res.surv=res.surv[-nrow(res.surv),]
    # NOTE: HERE IT'S the mean HAZARD RATE
    res.surv$hazard.rate=res.surv$hazard.rate.mean.from.amin
  } else if (grepl('steps',surv.model)){ # take the mean of the muti segments
    res.surv=res.surv[-nrow(res.surv),]
    res.surv$hazard.rate=rowMeans(res.surv[,grepl('hazard.rate',colnames(res.surv))])
  }
  
  # get cohort trend
  res=fn_get.cohort.trend(res.surv,dat) 
  fitted.cohort.trend=res$fitted.cohort.trend; est.cohort.eff=res$est.cohort.eff
  
  fn_plotFits(best.fit=surv.fits,dat,fig.title=fig.title1,add.axis=add.axis,add.legend=add.legend,add.ylab=add.ylab,add.xlab=add.xlab)
  fn_plotTrend(res.surv,fitted.cohort.trend,est.cohort.eff,fig.title=fig.title2,add.axis=add.axis,add.legend=add.legend,add.ylab=add.ylab,add.xlab=add.xlab,printRR=T)
  
  if(return.res==T)
    return(list(est.cohort.eff=est.cohort.eff, fitted.inci.rates=surv.fits, fitted.cohort.trend=fitted.cohort.trend))
}


##########################################################################################
# FUNCTION FOR GENERATING PREDICTIONS
##########################################################################################
# prediction using the 'acceleration' model
fn_get.pred.acc=function(res.surv,surv.model='fn_surv.reg.acc', N=1e5,age.min=25,age.max=49,num.pred=2,cohort.interval=5){
  # input: res.surv=output from the survival model
  
  age0=res.surv$age0[1];
  cohort.start=res.surv$cohort.start; 
  pred.start=tail(res.surv$cohort.end,1)+1;
  
  pred.cohorts=NULL  # get the cohorts for which predictions are made
  for(i in 1:num.pred) pred.cohorts=c(pred.cohorts,paste0(cohort.interval*(i-1)+pred.start,'-',pred.start+cohort.interval*i-1));
  
  if(res.surv$cohort.start[1]==res.surv$cohort.start[nrow(res.surv)]) res.surv=res.surv[-nrow(res.surv),] # if the last row is the mean across cohorts, delete it
  
  # 1. predict hr0 for later cohort
  fitted.h0=lm(hazard.rate0~cohort.start,data = res.surv) 
  pred.hr0=predict(fitted.h0,new = data.frame(cohort.start=seq(pred.start,length.out = num.pred,by=5)),se.fit = TRUE, interval = "confidence", level = 0.95)$fit
  
  fitted.delta=lm(half.hr.accel.rate~cohort.start,data = res.surv) # accel hr
  
  # make prediction based on the trend
  pred.delta=predict(fitted.delta,new = data.frame(cohort.start=seq(pred.start,length.out = num.pred,by=5)),se.fit = TRUE, interval = "confidence", level = 0.95)$fit
  ages=age.min:age.max
  dt=pmax(0,(ages-age0)); 
  
  pred.ir=pred.ir.lwr=pred.ir.upr=matrix(0,length(dt),num.pred)
  for(i in 1:num.pred){
    pred.ir[,i]=pmax(0,(1-exp(-pred.hr0[i,'fit']*dt -pred.delta[i,'fit']* dt^2))*N)
    pred.ir.lwr[,i]=pmax(0,(1-exp(-pred.hr0[i,'lwr']*dt -pred.delta[i,'lwr']* dt^2))*N)
    pred.ir.upr[,i]=pmax(0,(1-exp(-pred.hr0[i,'upr']*dt -pred.delta[i,'upr']* dt^2))*N)
  }
  colnames(pred.ir)=paste0('pred.ir',pred.cohorts);
  colnames(pred.ir.lwr)=paste0('pred.ir.lwr',pred.cohorts);
  colnames(pred.ir.upr)=paste0('pred.ir.upr',pred.cohorts);
  res.pred=data.frame(age=ages,pred.ir,pred.ir.lwr,pred.ir.upr )
  return(res.pred)
}
# prediction using the 'n-steps' model
fn_get.pred.nsteps=function(res.surv,surv.model='fn_surv.reg.2steps',N=1e5,age.min=25,age.max=49,num.pred=2){
  # input: res.surv=output from the survival model
  
  if(res.surv$cohort.start[1]==res.surv$cohort.start[nrow(res.surv)]) res.surv=res.surv[-nrow(res.surv),] # if the last row is the mean across cohorts, delete it
  age0=res.surv$age0[1];
  cohort.start=res.surv$cohort.start; 
  pred.start=tail(res.surv$cohort.end,1)+1;
  
  pred.cohorts=NULL  # get the cohorts for which predictions are made
  for(i in 1:num.pred) pred.cohorts=c(pred.cohorts,paste0(5*(i-1)+pred.start,'-',pred.start+5*i-1));
  
  hrs=res.surv[,grepl('hazard.rate',colnames(res.surv))]
  # make prediction for each age segment
  ages=age.min:age.max
  n.step=as.numeric(substr(surv.model,regexpr('steps',surv.model)[1]-1,regexpr('steps',surv.model)[1]-1))
  age.mids=seq(age.min,age.max,length.out = n.step+1)
  age.mids[1]=age0
  dt=matrix(0,age.max-age.min+1,n.step+1)
  dt[,1]=age.min:age.max
  for(j in ncol(dt):2){
    dt[,j]=pmax(dt[,1]-age.mids[j-1],0)
  }
  for(j in (ncol(dt)-1):2){
    dt[which(dt[,j+1]!=0),j]=age.mids[j]-age.mids[j-1]
  } 
  
  pred.hr.mean=pred.hr.lwr=pred.hr.upr=matrix(0,n.step,num.pred)
  for(j in 1:ncol(hrs)){
    y=unlist(hrs[j])
    tmp.fit=lm(y~cohort.start)
    tmp.pred=predict(tmp.fit,new = data.frame(cohort.start=seq(pred.start,length.out = num.pred,by=5)),se.fit = TRUE, interval = "confidence", level = 0.95)$fit
    pred.hr.mean[j,]= tmp.pred[,'fit']
    pred.hr.lwr[j,]=tmp.pred[,'lwr']
    pred.hr.upr[j,]=tmp.pred[,'upr']
  }
  
  pred.ir=pred.ir.lwr=pred.ir.upr=matrix(0,nrow(dt),num.pred)
  for(i in 1:num.pred){
    pred.ir[,i]=pmax(0,(1-exp(-(dt[,-1]%*%c(pred.hr.mean[,i]))))*N)
    pred.ir.lwr[,i]=pmax(0,(1-exp(-(dt[,-1]%*%c(pred.hr.lwr[,i]))))*N)
    pred.ir.upr[,i]=pmax(0,(1-exp(-(dt[,-1]%*%c(pred.hr.upr[,i]))))*N)
  }
  colnames(pred.ir)=paste0('pred.ir',pred.cohorts);
  colnames(pred.ir.lwr)=paste0('pred.ir.lwr',pred.cohorts);
  colnames(pred.ir.upr)=paste0('pred.ir.upr',pred.cohorts);
  res.pred=data.frame(age=ages,pred.ir,pred.ir.lwr,pred.ir.upr )
  return(res.pred)
}
