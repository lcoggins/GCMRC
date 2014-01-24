#Lew Coggins 11/2013
rm(list=ls(all=TRUE)) 
graphics.off()
windows(record=T)
setwd ("\\\\CCFHR-S-1720365\\Restrict3\\LCoggins\\Git\\GCMRC")
require(memisc)
#####################################
# Define Functions
#####################################

expit=function(x){1/(1+exp(-x))}

#this function figures out the number of MR trips based on the september abundance (x)
# getnumtrips=function(x){
#   y=lf2lcrFixedPars$MrTrigVec
#   nmtrips=cases(
#     x<y[1]->0,
#     x>=y[2] &x<y[3]->2,
#     x>=y[3] &x<y[4]->3,
#     x>=y[4] &x<y[5]->4,
#     x>=y[5] &x<y[6]->5,
#     x>=y[6] ->6
#   )
#   return(nmtrips)
# }
getnumtrips=function(x){
  y=lf2lcrFixedPars$MrTrigVec
  
  
  if(x<y[1]){nmtrips=0
    }else{
      if(x>=y[2] &x<y[3]){nmtrips=2
      }else{
        if(x>=y[3] &x<y[4]){nmtrips=3
        }else{                  
          if(x>=y[4] &x<y[5]){nmtrips=4
          }else{
            if(x>=y[5] &x<y[6]){nmtrips=5
            }else{
              if(x>=y[6]){nmtrips=6}}}}}}
   return(nmtrips)
}

#this function uses the number of trips (x) to assign the monthly mechanical removal vector
getMecvec=function(x){
     switch(x+1,
     c(0,0,0,0,0,0,0,0,0,0,0,0),
     c(0,0,0,0,0,0,1,0,0,0,0,0),
     c(0,0,0,0,0,1,1,0,0,0,0,0),
     c(0,0,0,0,1,1,1,0,0,0,0,0),
     c(0,0,0,1,1,1,1,0,0,0,0,0),
     c(0,0,1,1,1,1,1,0,0,0,0,0),
     c(0,1,1,1,1,1,1,0,0,0,0,0)
  )
}

#This function makes the movement matrix conditional on the variance of the normal kernal,
#the mean distance moved each month, zero inflation to have a higher probability of moving 
#the mean distance and a limit of the maximum monthly distance moved
getMoveMatrix=function(sqdim=length(lf2lcrFixedPars$rivm),
                       sigma=lf2lcrFixedPars$sigma,
                       logitzeroinf=lf2lcrFixedPars$logitzeroinf,
                       totreach=lf2lcrFixedPars$rivm,
                       toofar=lf2lcrFixedPars$toofar,
                       meandistmove=lf2lcrFixedPars$meandistmove){
  movemat=matrix(NA,nrow=sqdim,ncol=sqdim)
  zeroinf=expit(logitzeroinf)
  for(curloca in 1:sqdim){
    #get Basic Normal Movement from current location
    #moveprob=dnorm(totreach-curloca,meandistmove,sd=sigma)
    moveprob=dcauchy(totreach-curloca,meandistmove,scale=sigma)
    #Now add increased fidelity with the zero inflation
    moveprob[curloca]=moveprob[curloca]*(1+zeroinf)  #version1
    #moveprob[curloca]=moveprob[curloca+meandistmove]*(1+zeroinf) #version2
    #Eliminate movement if greater than defined maximum movement
    moveprob[abs(curloca-1:sqdim)>toofar]=0
    # Now renormalize
    moveprob=moveprob/sum(moveprob)
    movemat[curloca,]=moveprob
  }
  return(movemat)
}

##This Function creates the survival Matrix conditional on Mechanical Removal trigger
getSurvMatrix=function(MRtrig=0,
                       month,
                       MechVec=lf2lcrFixedPars$MechVec,
                       MRpcap=lf2lcrFixedPars$MRpcap,
                       MRpass=lf2lcrFixedPars$MRpass,
                       UppEndMR=lf2lcrFixedPars$UppEndMR,
                       LowEndMR=lf2lcrFixedPars$LowEndMR,
                       sqdim=length(lf2lcrFixedPars$rivm),
                       S=lf2lcrFixedPars$S){
  SurvMR=1
  #check if MR authorized
  if(MRtrig==1&MechVec[month]==1){SurvMR=(1-MRpcap)^MRpass}
  Surv=rep(1,sqdim)
  Surv[UppEndMR:LowEndMR]=Surv[UppEndMR:LowEndMR]*SurvMR
  #Now Add natural Mortality to the survival vector
  Surv=S*Surv;#print(Surv)
  #Now create and return the survival matrix
  SurvMatrix=matrix(0,nrow=sqdim,ncol=sqdim)
  diag(SurvMatrix)=Surv
  return(SurvMatrix)
}


#This function Does the  moving and surviving and keeping track of the important output
#it will need to be supplied with the current year, status of the trout trigger for that year,
#and the rbt abundance by river mile the last day of the previous year 
MoveTrout=function(year,
                   MRtrig=0,
                   rbtvec,
                   MechVec=lf2lcrFixedPars$MechVec,
                   S=lf2lcrFixedPars$S,
                   MRpass=lf2lcrFixedPars$MRpass,
                   MRpcap=lf2lcrFixedPars$MRpcap,
                   UppEndMR=lf2lcrFixedPars$UppEndMR,
                   LowEndMR=lf2lcrFixedPars$LowEndMR,
                   AnnRec=lf2lcrFixedPars$AnnRec, #this is the default if no outmigration pattern is specified
                   rivm=lf2lcrFixedPars$rivm,
                   meandistmove=lf2lcrFixedPars$meandistmove){
  
#Add the Annual Recruitment
rbtvec[1]=rbtvec[1]+AnnRec/12
#Set up a matrix to keep the Annual Abundance by Reach
rbtmat=matrix(0,ncol=7,nrow=12)
#Set up the vector to get the RBT Abundance in the LCR Reach
RbtIn4b=rep(NA,12)


#Now get the Movement Matrix
MMat=getMoveMatrix()
#matplot(t(MMat),type='l',ylab="Proportion following Movement", xlab='River Mile')

#start the Loop by Month
for(mon in 1:12){

#Add Monthly Recruitment to River Mile 1 
rbtvec[1]=rbtvec[1]+AnnRec/12  

#Now construct the survival matrix including Mechanical Removal if Authorized
SMat=getSurvMatrix(month=mon,MRtrig=MRtrig)
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
#RbtIn4b[mon]=rbtmat[mon,7]
RbtIn4b[mon]=sum(rbtvec[64]*.6,rbtvec[65]*0.8)


#print(rbtvec)
#print(rbtmat)
#plot(rbtvec,ylim=c(0,1000))
}

#matplot(t(MMat),type='l')
#matplot(t(rbtmat),type='l')
return(list(rbtvec=rbtvec,rbtmat=rbtmat,RbtIn4b=RbtIn4b))
}

#End Function Definitions
##########################################################
##########################################################




###########################################################
############################################################
#START THE Main PROGRAM Segment
############################################################
#specify all of the parameters and the initial conditions


lf2lcrFixedPars=list(InitVec=rep(10,151), #This is the Initial RBT Abundance
#lf2lcrFixedPars=list(InitVec=5*(151:1), #This is the Initial RBT Abundance
                     MechVec=rep(c(rep(1,3),rep(0,3)),2), #this is the bit that specifies which months CAN have mechanical removal
                     S=0.96,
                     MRpass=5,    #Number of MechRem passes
                     MRpcap=0.1,  #capture probability of each MR pass
                     UppEndMR=56, #upper end of the removal reach (must be in whole miles)
                     LowEndMR=66, #lower end of the removal reach (must be in whole miles)
                     AnnRec=5000, #this is the default if no outmigration pattern is specified
                     rivm=1:151,  #These are the ending river miles for the 1 mile reaches that we will model 
                     sigma=3.38,  # use for cauchy diffusion
                     #sigma=9.91,  #use for normal diffusion
                     logitzeroinf=-20.0, #used to adjust "fidelity" to the average amount of movement
                     toofar=100,  #used to specify the max movement in a month
                     meandistmove=0, #used to specify an average amount of movement in miles
                     startyr=1990,
                     duration=19,
                     MrTrigVec=c(760,760,926,1649,2855,4864)
                     ) 

obsOutmig=c(7500.948, 2965.82, 10838.41, 24784.04, 58564.31,
            37788.56, 43309.8, 88600.62, 67744.24, 70546.81,
            44704.52, 12792.74, 18986.56, 22107.93, 14111.24,
            1517.380, 10834.17, 99020.74, 131703.6, 134571.0,NA)

# ParVar=0  #allow variation in parameter values for some funcitonal relationships among trials (0=no, >0=yes (# of trials))
# IaVar=0  #interannual variation in predictions from some functional relationships (0-no,1=yes)
# if(IaVar==0){
#   Pstay_sd=0.000001
#   
# } else {
#   Pstay_sd=lf2crFixedPars$Pstay_sd
# }  


#Define the initial rbt abundance by mile
rbtN=lf2lcrFixedPars$InitVec

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Make a vector to keep track of the midyear abundance by 10 mile section
#This can be used to "tune" the movement parameters to fit the historical patterns
RbtByYear=matrix(NA,nrow=lf2lcrFixedPars$duration+1,ncol=7)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Now Loop over years.  
#This is meant to be akin to how these functions will be incorporated into the full simulation
for(yr in lf2lcrFixedPars$startyr:(lf2lcrFixedPars$startyr+lf2lcrFixedPars$duration)){
remreachN=0
MRTRIG=0
if(yr>lf2lcrFixedPars$startyr){remreachN=sum(AnnResults$rbtvec[63:64])+0.5*AnnResults$rbtvec[65]}
AduChubN=runif(1,6000,6000) #this line is just a placeholder for using predicted HBC abundance
if(AduChubN<7000&remreachN>760){MRTRIG=1}
numtrips=getnumtrips(remreachN)
Mecvec=getMecvec(numtrips)*MRTRIG
print(remreachN)
print(numtrips)
print(Mecvec)



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#AnnResults contains the list of model output
AnnResults=MoveTrout(year=yr,rbtvec=rbtN,MechVec=Mecvec,AnnRec=obsOutmig[yr-(lf2lcrFixedPars$startyr-1)],MRtrig=MRTRIG)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rbtN=AnnResults$rbtvec
#this is the matrix that can be used to tune the model
RbtByYear[(yr-(lf2lcrFixedPars$startyr-1)),]=AnnResults$rbtmat[6,]

#Use the lines below to see some of the output
plot(1:length(lf2lcrFixedPars$rivm),AnnResults$rbtvec,type='l',ylim=c(0,7500),
     main=paste('this is the abundance at the end of December ',yr),xlab='rivermile')
#print(AnnResults$RbtIn4b)
#print(AnnResults$rbtmat[6,])
#plot(1:7,AnnResults$rbtmat[6,],type='l',ylim=c(0,20000),
     #main=paste('these are the 10mile reach abundances at the end of July ',yr),xlab='10mile Reach')
#Add the midyear abundances to the tuning matrix

}

#END THE MAIN PROGRAM SEGMENT HERE
##########################################################
##########################################################

cbind(lf2lcrFixedPars$startyr:(lf2lcrFixedPars$startyr+lf2lcrFixedPars$duration),RbtByYear)
plot(lf2lcrFixedPars$startyr:(lf2lcrFixedPars$startyr+lf2lcrFixedPars$duration),RbtByYear[,7],pch=19)


