#Lew Coggins 11/2013
rm(list=ls(all=TRUE)) 
graphics.off()
windows(record=T)
setwd ("\\\\CCFHR-S-1720365\\Restrict3\\LCoggins\\2013 Work\\GCSDM")

#prepare the historical data for comparison with the predictions
histdat=read.csv(file='rbt_hist_sum.csv')
histdat=subset(histdat,reach!=" RM < 0")
histdat=subset(histdat,reach!=" RM 65.7-76.3") 
histdat=cbind(histdat,rep(1:7,length(unique(histdat$year))));names(histdat)[9]='reachnum'
histdat=subset(histdat,select=c(year,reachnum,Nsamples,EF,catch_all))
#histdat=cbind(histdat,histdat[,4],histdat[,5]);names(histdat)[6]='effort';names(histdat)[7]='catch'
#histdat=subset(histdat,select=c(year,reachnum,Nsamples,CPUE))
names(histdat)[5]='catch'
histdat=subset(histdat,year>1989)
histdat=subset(histdat,year<2010)
histdat=cbind(histdat,rep(NA,dim(histdat)[1]));names(histdat)[6]='predCatch'

#####################################
# Define Functions
#####################################

expit=function(x){1/(1+exp(-x))}

##################################################################################
####################This function attempts to find the best fit to 2000-2009 data
##################################################################################
tune=function(pars){
Q=pars[1]
#Q=Qinit
sigma=pars[2]
#sigma=siginit
#mdm=mdminit
mdm=pars[3]
logitzeroinf=-20
zeroinf=expit(logitzeroinf)
#print(paste('Q is ',Q))
#print(paste('sigma is ',sigma))
#print(paste('mdm is ',mdm))
#print(paste('zeroinf is ',zeroinf))
Q=Q/1000000
mdm=abs(mdm)
# sigma=sig
rbtN=lf2lcrFixedPars$InitVec
RbtByYear=matrix(NA,nrow=lf2lcrFixedPars$duration+1,ncol=7)
remreachN=0  
MRTRIG=c(rep(0,13),rep(1,4),rep(0,2),1)
MechVec=lf2lcrFixedPars$MechVec
S=lf2lcrFixedPars$S
MRpass=lf2lcrFixedPars$MRpass
MRpcap=lf2lcrFixedPars$MRpcap
UppEndMR=lf2lcrFixedPars$UppEndMR
LowEndMR=lf2lcrFixedPars$LowEndMR
AnnRec=lf2lcrFixedPars$AnnRec #this is the default if no outmigration pattern is specifie
rivm=lf2lcrFixedPars$rivm

  
rbtvec= lf2lcrFixedPars$InitVec=rep(50,151)   
  #Set up a matrix to keep the Annual Abundance by Reach
rbtmat=matrix(0,ncol=7,nrow=12)
  #Set up the vector to get the RBT Abundance in the LCR Reach
RbtIn4b=rep(NA,12)
  
    
    #Now get the Movement Matrix
    sqdim=length(lf2lcrFixedPars$rivm)
    totreach=lf2lcrFixedPars$rivm
    toofar=lf2lcrFixedPars$toofar
      movemat=matrix(NA,nrow=sqdim,ncol=sqdim)
      for(curloca in 1:sqdim){
        #get Basic Normal Movement from current location
        #moveprob=dnorm(totreach-curloca,mdm,sd=sigma)
        moveprob=dcauchy(totreach-curloca,location=mdm,scale=sigma)
        #Now add increased fidelity with the zero inflation
        moveprob[curloca]=moveprob[curloca]*(1+zeroinf)  #version1
        #moveprob[curloca]=moveprob[curloca+mdm]*(1+zeroinf) #version2
        #Eliminate movement if greater than defined maximum movement
        moveprob[abs(curloca-1:sqdim)>toofar]=0
        # Now renormalize
        moveprob=moveprob/sum(moveprob)
        movemat[curloca,]=moveprob
      }
    MMat=movemat
    RbtIn4b=numeric()
    #matplot(t(MMat),type='l')


for(yr in lf2lcrFixedPars$startyr:(lf2lcrFixedPars$startyr+lf2lcrFixedPars$duration)){
      AnnRec= obsOutmig[yr-(lf2lcrFixedPars$startyr-1)] 
      #start the Loop by Month
      for(mon in 1:12){
        
        #Add Monthly Recruitment to River Mile 1 
        rbtvec[1]=rbtvec[1]+AnnRec/12  
        
    #Now construct the survival matrix including Mechanical Removal if Authorized
    #SMat=getSurvMatrix(month=mon,MRtrig=MRtrig)
    
    MRtrig=MRTRIG[yr-(lf2lcrFixedPars$startyr-1)]
    MechVec=lf2lcrFixedPars$MechVec
    MRpcap=lf2lcrFixedPars$MRpcap
    MRpass=lf2lcrFixedPars$MRpass
    UppEndMR=lf2lcrFixedPars$UppEndMR
    LowEndMR=lf2lcrFixedPars$LowEndMR
    sqdim=length(lf2lcrFixedPars$rivm)
    S=lf2lcrFixedPars$S
      SurvMR=1
      #check if MR authorized
      if(MRtrig==1&MechVec[mon]==1){SurvMR=(1-MRpcap)^MRpass}
      Surv=rep(1,sqdim)
      Surv[UppEndMR:LowEndMR]=Surv[UppEndMR:LowEndMR]*SurvMR
      #Now Add natural Mortality to the survival vector
      Surv=S*Surv;#print(Surv)
      #Now create and return the survival matrix
      SurvMatrix=matrix(0,nrow=sqdim,ncol=sqdim)
      diag(SurvMatrix)=Surv
    SMat=SurvMatrix    
    #matplot(t(SMat),type='l')
    
    #print(diag(SMat))
    #print(SMat)
    #print(dim(SMat))
    
    #Now Survive and Move the fish
    #print(rbtvec);print(dim(rbtvec))
    rbtvec=rbtvec%*%SMat#;print(rbtvec)
    rbtvec=round(rbtvec%*%MMat,0)#;print(rbtvec)
    
    rbtmat[mon,]=c(sum(rbtvec[1:10]),
                   sum(rbtvec[11:20]),
                   sum(rbtvec[21:30]),
                   sum(rbtvec[31:40]),
                   sum(rbtvec[41:50]),
                   sum(rbtvec[51:61],rbtvec[62]*0.5),
                   sum(rbtvec[62]*0.5,rbtvec[63:65],rbtvec[66]*0.7))
    RbtIn4b<-c(RbtIn4b,sum(rbtvec[64]*0.6,rbtvec[65]*.8))
        
    #print(rbtvec)
    #print(rbtmat)
    #plot(rbtvec,ylim=c(0,1000))
  }
  

rbtN=rbtvec #set next january 1 to december 31
RbtByYear[(yr-(lf2lcrFixedPars$startyr-1)),]=rbtmat[6,]#;print(RbtByYear)
#print(yr)  
}
#matplot(t(MMat),type='l')

rbt.jcm<<-RbtIn4b
Nbyreach<<-cbind(histdat$year,as.vector(t(RbtByYear)))
pred<<-cbind(histdat$year,as.vector(t(RbtByYear))*Q*histdat[,4])

obs<<-cbind(histdat$year,histdat[,5])
predsub=subset(pred,pred[,1]>1999&is.na(pred[,2])==FALSE);obssub=subset(obs,obs[,1]>1999&is.na(pred[,2])==FALSE)
predsub=predsub[predsub[,1]<2010,];obssub=obssub[obssub[,1]<2010,]
#ssq=na.omit((predsub-obssub)^2)
#ssq=sum(ssq)*100
#return(ssq=ssq)
lik<-dpois(obssub,predsub,log=T)
llik=sum(-1*lik)
return(llik)

}

#End Function Definitions
##########################################################
##########################################################






###########################################################
############################################################
#START THE Main PROGRAM
############################################################
#specify all of the parameters and the initial conditions


lf2lcrFixedPars=list(InitVec=rep(50,151),#1000*(15:1), #This is the Initial RBT Abundance
                     MechVec=rep(c(rep(1,3),rep(0,3)),2), #this is the bit that specifies which months CAN have mechanical removal
                     S=0.96,
#                     S=0.92,
                     MRpass=5,    #Number of MechRem passes
                     MRpcap=0.10,  #capture probability of each MR pass
                     UppEndMR=56, #upper end of the removal reach (must be in whole miles)
                     LowEndMR=66, #lower end of the removal reach (must be in whole miles)
                     AnnRec=5000, #this is the default if no outmigration pattern is specified
                     rivm=1:151,  #These are the ending river miles for the 1 mile reaches that we will model 
                     sigma=6.0,    #diffusion
                     zeroinf=0.0, #used to adjust "fidelity" to the average amount of movement
                     toofar=100.0,   #used to specify the max movement in a month
                     mdm=0.0, #used to specify an average amount of movement in miles
                     startyr=1990,
                     duration=19) 

obsOutmig=c(7500.948, 2965.82, 10838.41, 24784.04, 58564.31,
            37788.56, 43309.8, 88600.62, 67744.24, 70546.81,
            44704.52, 12792.74, 18986.56, 22107.93, 14111.24,
            1517.380, 10834.17, 99020.74, 131703.6, 134571.0,NA)


Qinit=3.4
siginit=3.4
mdminit=0
zeroinfinit=-20
parms=c(Qinit,siginit,mdminit)
#parms=c(Qinit,siginit)
#parms=c(Qinit)
#parms=c(siginit,mdminit)
#parms=c(siginit,zeroinfinit)
#parms=c(siginit)

res=optim(parms,tune,control=list(trace=TRUE))
#res=optim(parms,tune,method='BFGS',control=list(trace=TRUE))
#optim(parms,tune,method='Brent',lower=1,upper=30,control=list(trace=TRUE))

#histdat[,4]=round(obs[,2],0)
#histdat[,5]=round(pred[,2],0)
histdat[,5]=obs[,2]
histdat[,6]=pred[,2]
histdat

histdattmp=subset(histdat,year>1999)
histdattmp=subset(histdattmp,year<2010)

windows(record=T)
op=par(mfrow=c(2,2))
ptcex=.6
noxlab=c(2000,2001,2004,2005)
for(i in 2000:2009){
  x=histdat[histdat[,1]==i,2];y2=histdat[histdat[,1]==i,5];y=histdat[histdat[,1]==i,6]
  test=!any(i==noxlab);print(test)
  if(test){
    plot(x,y,type='b',main=i,ylim=c(0,max(c(y,y2),na.rm = TRUE)),axes=F,ylab='Annual Catch',xlab='Reach')
  }else{
    plot(x,y,type='b',main=i,ylim=c(0,max(c(y,y2),na.rm = TRUE)),axes=F,ylab='Annual Catch',xlab=' ')
  }
  axis(2,las=1)
  axis(1,at=1:7,
       labels=c('RM 0-10','RM 10-20','RM 20-30','RM 30-40','RM 40-50','RM 50-61.5','RM 61.5-65.7'),
       cex.axis=.8,las=2,tick=F,line=-.9)
  box()
  points(x,y2,pch=3,cex=ptcex)
  if(i==2008 | i==2007){
    legend('topright',legend=c('predicted catch','observed catch'),pch=c(1,3),lty=c(1,0),pt.cex=c(1,ptcex))
  }else if (i==2009){
    legend('bottom',legend=c('predicted catch','observed catch'),pch=c(1,3),lty=c(1,0),pt.cex=c(1,ptcex))
  }else{
    legend('topleft',legend=c('predicted catch','observed catch'),pch=c(1,3),lty=c(1,0),pt.cex=c(1,ptcex))
    
  }  
}
par=op
#res
