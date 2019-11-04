# This is to use simulation to test the ability of survival modeling to evaluate cohort-effect
# The method is tested on 5 scenarios:
# 1) Cohort effect alone, for which both the baseline risk and aging-related risk are increased by 10% per later cohort. 
# 2) Period effect due to changes in exposure alone, for which, to mimic the folate intake and colorectal cancer risk example, the aging-related risk is reduced by 10% for all ages from 1998 and partially reduced for all ages between 1973 and 1997 based on the prevalence of folate supplement use. 
# 3) Period effect due to changes in detection alone, for which, to mimic the CT use and thyroid cancer diagnosis example, the diagnosis rate and thus incidence are increased based on the prevalence of CT use over time. 
# 4) Cohort effect combining period effect due to changes in exposure, i.e., combining 1 and 2. 
# 5) Cohort effect combining period effect due to changes in detection, i.e., combining 1 and 3.  

# test 5 different cohorts
# add period-specific "noise" to represent period effect 
# age-effect: survial model

library(RColorBrewer); library(data.table); library(magrittr); library(stringr);

source('./Functions_to_analyze_cancer_trends_by_cohort.R') # load the functions to analyze cohort trends


surv.model='fn_surv.reg.acc' # use the acceleration model
fn_surv.reg=get(surv.model) # get the model

cohort.width=5; # cohort width: 5 year
age.min=25; age.max=49 # set the age range
cohorts=seq(1945,1965,by=cohort.width) # set cohorts, those born from 1945 to 1969
ages=0:age.max

# changes for the cohort effect
ch.incr.rate=.1/5*cohort.width # per cohort % increase [10% for 5-yr]

# Simulaiton of Folic Acid Fortification and Colorectal Cancer Risk
# changes for the first period effect (based on folate intake)
# From: Keum N, Giovannucci EL. Folic acid fortification and colorectal cancer risk. Am J Prev Med. 2014;46(3 Suppl 1):S65-72. PubMed PMID: 24512932.
# In 1998, the U.S. Food and Drug Administration required that folic acid be added to enriched grain products (such as bread, pasta, rice, and cereal). This is called fortification.
# According to the National Health and Nutrition Examination Survey I (1971–1975), the prevalence of supplement use among the adults aged 20–74 years was 27.5% in men and 37.9% in women, and whites were the major users
period.rate=-.1 # period effect (10% reduction in hazard rate if with 100% intake)
eff.year1=1973 # start year
eff.year2=2100 # end year
affected.ages=0:age.max # all ages are affected
# prevalence of fortificaiton: 0% before 1970, 25% with 1% increase per year form 1973 to 1997; 100% from 1998 onwards
prev=data.frame(year=1920:2020,prev=c(rep(0,1972-1920+1),.25+.01*(1973:1997-1973),rep(1,2020-1998+1)))

# Simulation of widespread CT use and increased cancer detection (e.g. thyroid cancer)
# changes for the 2nd period effect (based on CT use)
# according to Mettler et al. 2010, in the US, frequency of diagnostic radiologic exam has increased almost 10-fold from 1950 to 2006
# Brenner 2010 JAMA: "In 1980 fewer than 3 million CT scans were performed, but the annual number now approaches 80 million and is increasing by approximately 10% per year."
# based on Hess et al. 2014 CT scan increased by 10% from 2000 to 2010 -> 1% increase per year
eff.year3=1980  # if it's after 1980, increase the rate
eff.year4=2005  # peak in 2005; reduce afterwards or detection rate levels off
incr.rate=.1; # % increase per year
dr0=.5 # basic detection rate


# Use the acceleration model as an example
# acceleration model, hazard at age a: h(a)=h0 + delta * (a-a0)
# -> cumulative hazard: H(a)=h0 * (a-a0) + 1/2 * delta * (a-a0)^2
# -> Incidence rate per N=100,000: (1 - S(a)) * N = (1 - exp(-H(a)))*N
N=1e5; # per 100,000 person-year
a0=15; # age with nonzero hazard
h0=2e-6; # initial hazard rate
delta0=3e-7 # increase in harzard rate per 1-year age

# Generate test data
############################################
# 1. with cohort effect only
h0s=h0*seq(1,length.out = length(cohorts),by=ch.incr.rate) # initial hazard rate for each cohort
deltas=delta0*seq(1,length.out = length(cohorts),by=ch.incr.rate) # hazard increase rate for each cohort
IR=matrix(0,length(ages),length(cohorts))
for(ic in 1:length(cohorts)){
  
  a=pmax(ages,a0)
  H=h0s[ic] * (a-a0) + 1/2 * deltas[ic] * (a-a0)^2  # compute the cumulative hazard
  ir = (1 - exp(-H))*N  # compute the true total cancer rate including undetected
  
  IR[,ic]=ir * dr0; # assume 50% detection rate, for comparison with Scenario 3 and 5
}

colnames(IR)=apply(cbind(cohorts,cohorts+cohort.width-1),1,paste0,collapse='-')
dat=cbind(Age=ages,IR) %>% data.table()
dat=dat[Age<=age.max & Age>=age.min]
# save it
dat1=dat

# plot data
cc=palette(brewer.pal(n=9,name='OrRd')); colfunc=colorRampPalette(cc);cols=colfunc(length(cohorts)+2)[-c(1,2)]
matplot(IR,type='l',lty=1,col=cols,xlim=c(age.min,age.max))


############################################
# 2. Period effect due to changing exposure
IR=matrix(0,length(ages),length(cohorts))
for(ic in 1:length(cohorts)){
  cohort=cohorts[ic]
  cohort.mid=cohort+cohort.width/2 # mid.year of the cohhort
  
  a=pmax(ages,a0)
  H=h0 * (a-a0) + 1/2 * delta0 * (a-a0)^2 # cumulative harzard for each age
  ir.a = (1 - exp(-H))*N # due to aging 
  
  # add harzard due to exposure in certain period
  c.years=ages+cohort.mid  # the corresponding years of detection
  period.eff=data.table(year=floor(c.years),age=ages,p.eff=0)
  period.eff=merge(period.eff,prev,by='year',all.x=T) # include the prevalence of exposure
  # adjust the period effect based on calendar year
  period.eff[year %in% eff.year1:eff.year2 & (age %in% affected.ages) & age>a0]$p.eff=period.rate*period.eff[year %in% eff.year1:eff.year2 & (age %in% affected.ages) & age>a0]$prev  
  p.effs=period.eff$p.eff
  
  H.period = h0 * p.effs * pmax(0,c.years-eff.year1) # get the period effect, it's a change related to the base-rate hazard (h0)
  
  ir = (1 - exp(-H-H.period))*N # add period effect
  
  IR[,ic]=ir
}
colnames(IR)=apply(cbind(cohorts,cohorts+cohort.width-1),1,paste0,collapse='-')
dat=cbind(Age=ages,IR) %>% data.table()
dat=dat[Age<=age.max & Age>=age.min]
# save it
dat2=dat

# plot
cc=palette(brewer.pal(n=9,name='OrRd')); colfunc=colorRampPalette(cc);cols=colfunc(length(cohorts)+2)[-c(1,2)]
matplot(IR,type='l',lty=1,col=cols,xlim=c(age.min,age.max))

############################################
# 3. Period effect due to changing detection rate 
IR=matrix(0,length(ages),length(cohorts))
for(ic in 1:length(cohorts)){
  cohort=cohorts[ic]
  cohort.mid=cohort+cohort.width/2 # mid.year of the cohhort
  
  a=pmax(ages,a0)
  H=h0 * (a-a0) + 1/2 * delta0 * (a-a0)^2
  ir.a = (1 - exp(-H))*N # due to aging 

  c.years=ages+cohort.mid  # the corresponding years of detection

  period.eff=data.table(year=floor(c.years),num_CT=1,p.eff=1)
  period.eff[year %in% eff.year3:eff.year4]$num_CT=(1+incr.rate)^(period.eff[floor(year) %in% eff.year3:eff.year4]$year-eff.year3) # 10% increase from 1980 to 2005
  period.eff[floor(year) > eff.year4]$num_CT= period.eff[floor(year) == eff.year4]$num_CT*(1+incr.rate/10)^(period.eff[floor(year) > eff.year4]$year-eff.year4) # slower increase after 2005
  # get probability of diagnosis, assuming a 50% positivity rate per CT scan
  period.eff$p.eff = 1 - (1-dr0)^period.eff$num_CT
  
  p.effs=period.eff$p.eff
  
  ir = ir.a * p.effs
  IR[,ic]=ir
}
colnames(IR)=apply(cbind(cohorts,cohorts+cohort.width-1),1,paste0,collapse='-')
dat=cbind(Age=ages,IR) %>% data.table()
dat=dat[Age<=age.max & Age>=age.min]

# save it
dat3=dat

# plot and see
cc=palette(brewer.pal(n=9,name='OrRd')); colfunc=colorRampPalette(cc);cols=colfunc(length(cohorts)+2)[-c(1,2)]
matplot(IR,type='l',lty=1,col=cols,xlim=c(age.min,age.max))


##########################
# 4. Exposure period effect + cohort effect
h0s=h0*seq(1,length.out = length(cohorts),by=ch.incr.rate) # initial hazard rate for each cohort
deltas=delta0*seq(1,length.out = length(cohorts),by=ch.incr.rate) # hazard increase rate for each cohort
IR=matrix(0,length(ages),length(cohorts))
for(ic in 1:length(cohorts)){
  cohort=cohorts[ic]
  cohort.mid=cohort+cohort.width/2 # mid.year of the cohhort
  
  a=pmax(ages,a0)
  H=h0s[ic] * (a-a0) + 1/2 * deltas[ic] * (a-a0)^2
  ir.a = (1 - exp(-H))*N  # due to aging 
  
  # add harzard due to exposure in certain period
  c.years=ages+cohort.mid  # the corresponding years of detection
  period.eff=data.table(year=floor(c.years),age=ages,p.eff=0)
  period.eff=merge(period.eff,prev,by='year',all.x=T) # include the prevalence of exposure
  # adjust the period effect based on calendar year
  period.eff[year %in% eff.year1:eff.year2 & (age %in% affected.ages) & age>a0]$p.eff=period.rate*period.eff[year %in% eff.year1:eff.year2 & (age %in% affected.ages) & age>a0]$prev  
  p.effs=period.eff$p.eff
  
  H.period = h0 * p.effs * pmax(0,c.years-eff.year1) # get the period effect, it's a change related to the base-rate hazard (h0)
  
  ir = (1 - exp(-H -H.period))*N
  
  IR[,ic]=ir
}
colnames(IR)=apply(cbind(cohorts,cohorts+cohort.width-1),1,paste0,collapse='-')
dat=cbind(Age=ages,IR) %>% data.table()
dat=dat[Age<=age.max & Age>=age.min]
# save it
dat4=dat

# plot and see
cc=palette(brewer.pal(n=9,name='OrRd')); colfunc=colorRampPalette(cc);cols=colfunc(length(cohorts)+2)[-c(1,2)]
matplot(IR,type='l',lty=1,col=cols,xlim=c(age.min,age.max))

##########################
## 5. multiplicative: cohort effect & period effect
h0s=h0*seq(1,length.out = length(cohorts),by=ch.incr.rate) # initial hazard rate for each cohort
deltas=delta0*seq(1,length.out = length(cohorts),by=ch.incr.rate) # hazard increase rate for each cohort
IR=matrix(0,length(ages),length(cohorts))
for(ic in 1:length(cohorts)){
  cohort=cohorts[ic]
  cohort.mid=cohort+cohort.width/2 # mid.year of the cohhort
  
  a=pmax(ages,a0)
  H=h0s[ic] * (a-a0) + 1/2 * deltas[ic] * (a-a0)^2  # cumulative hazard for that cohort
  ir.a = (1 - exp(-H))*N  # due to aging for that cohort
  
  c.years=ages+cohort.mid  # the corresponding years of detection
  
  period.eff=data.table(year=floor(c.years),num_CT=1,p.eff=1)
  period.eff[year %in% eff.year3:eff.year4]$num_CT=(1+incr.rate)^(period.eff[floor(year) %in% eff.year3:eff.year4]$year-eff.year3) # 10% increase from 1980 to 2005
  period.eff[floor(year) > eff.year4]$num_CT= period.eff[floor(year) == eff.year4]$num_CT*(1+incr.rate/10)^(period.eff[floor(year) > eff.year4]$year-eff.year4) # slower increase after 2005
  # get probability of diagnosis, assuming a 50% positivity rate per CT scan
  period.eff$p.eff = 1 - (1-dr0)^period.eff$num_CT
  
  p.effs=period.eff$p.eff
  
  ir = ir.a * p.effs
  
  IR[,ic]=ir
}
colnames(IR)=apply(cbind(cohorts,cohorts+cohort.width-1),1,paste0,collapse='-')
dat=cbind(Age=ages,IR) %>% data.table()
dat=dat[Age<=age.max & Age>=age.min]
# save it
dat5=dat

# plot and see
cc=palette(brewer.pal(n=9,name='OrRd')); colfunc=colorRampPalette(cc);cols=colfunc(length(cohorts)+2)[-c(1,2)]
matplot(IR,type='l',lty=1,col=cols,xlim=c(age.min,age.max))

################################################################
scenarios=c('(A) Cohort effect alone','(B) Period effect: exposure','(C) Period effect: detection','(D) Period exposure + cohort eff.','(E) Period detection + cohort eff.')
scenarios=c('Cohort effect alone',"Period effect (exposure)",'Period effect (detection)','Period/exposure + cohort','Period/detection + cohort')

num.pred=2; cols=c('blue','red'); pchs=c('x','+')
n.row=5;n.col=2

pdf(paste0('./Fig_sim_fit.trend_by.',surv.model,'_a',age.min,'-',age.max,'.pdf'),width=7,height = 1.6*n.row)
par(mfrow=c(n.row,n.col),mar=c(0,2.2,0,.5),oma=c(2.5,.5,1.5,.5),cex=.8,cex.lab=.9,cex.axis=.8,mgp=c(1,.2,0),tck=-.02)
for(i in 1:5){

  fig.title1=paste0('(',LETTERS[2*(i-1)+1],') ',scenarios[i],': Survival model-fit'); # paste0(site,': ',sex.names[sex])
  fig.title2=paste0('(',LETTERS[2*(i-1)+2],') ',scenarios[i],': Cohort trend'); 
  dat=get(paste0('dat',i)); 
  add.ylab=T; add.note=F
  add.axis=ifelse(i==5,T,F); add.xlab=ifelse(i==5,T,F); 
  add.legend=ifelse(i==1,T,F);
  fn_fit.trend(dat,surv.model=surv.model,fig.title1=fig.title1,fig.title2=fig.title2,add.axis=add.axis,add.legend=add.legend,add.ylab=add.ylab,add.xlab=add.xlab,printRR=T,return.res=F)
  
}
dev.off()

