mu=0
location=0
scale=3.38
sigma=9.91
x=seq(-20,20,by=.1)
dcauch=dcauchy(x,location,scale)
#dcauch=dcauch/sum(dcauch)
plot(x,dcauch,type='l',xlab='Movement Distance (miles)',ylab='Probability',lwd=2)
#mode is:
x[dcauch==max(dcauch)]
lines(x,dnorm(x,mu,sigma),lty=3)
legend('topright',legend=c('Cauchy Movement Distribution','Normal Movement Distribution'),
       lty=c(1,3),lwd=c(2,1),bty='n')
rcauch=rcauchy(1000000,location,scale)
#plot(density(rcauch),xlim=range(x))
median(rcauch)
mean(rcauch)
var(rcauch)
