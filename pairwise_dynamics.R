setwd("/Users/eleanorbrush/Dropbox/signaling_network")
source("gillespie.R")
library(fields)
library(lattice)

e=5
bs=2
bf=1
rf=1
dominance=.5
leak=.2
ramp=.5
maxtime=100

threshvals=seq(0,20,by=.5)
Lt=length(threshvals)
domvals=seq(.5,1,by=0.1)
Ld=length(domvals)

exp.waits=array(,c(Lt,Ld))
prob.weaker.sigs=array(,c(Lt,Ld))
trans.probs=array(,c(Lt,Ld,4,4))
trans.times=array(,c(Lt,Ld,4,4))
for(i in 1:Lt){
	for(j in 1:Ld){
		thresh=threshvals[i]
		dominance=domvals[j]
		simulation=sims(e,maxtime,ramp,thresh,thresh,rf,dominance,leak,bs,bf)
		s=stats(simulation)
		exp.waits[i,j]=s$exp.wait
		prob.weaker.sigs[i,j]=s$prob.weaker.sigs
		trans.probs[i,j,,]=s$trans.probs
		trans.times[i,j,,]=s$trans.times
		
	}
}
exp.waits[is.nan(exp.waits)]=NA
prob.weaker.sigs[is.nan(prob.weaker.sigs)]=NA

L=list(threshvals=threshvals,domvals=domvals,exp.waits,prob.weaker.sigs=prob.weaker.sigs,trans.probs=trans.probs,trans.times=trans.times)
save(file="/Users/eleanorbrush/Desktop/signetsims.rdata",L)

domvals=c(1-domvals[length(domvals):2],domvals)

waits=L$exp.waits
waits[is.na(waits)]=maxtime
waits=cbind(waits[,dim(waits)[2]:2],waits)
rownames(waits)=paste("thresh",threshvals)
colnames(waits)=paste("dom",domvals)

focalreceives=L$prob.weaker.sigs
focalreceives[is.na(focalreceives)]=0.5
focalreceives=cbind(1-focalreceives[,dim(focalreceives)[2]:2],focalreceives)
rownames(right)=paste("thresh",threshvals)
colnames(right)=paste("dom",domvals)

layout(matrix(1:2,nrow=1))
image.plot(threshvals,domvals,waits,xlab="Threshold",ylab="Probability that focal animal wins a fight")
image.plot(threshvals,domvals,focalreceives,xlab="Threshold",ylab="Probability that focal animal wins a fight")

weights=rep(1/length(domvals),length(domvals))
weights=exp(-10*abs(domvals-.5))
weights=dnorm(seq(0,1,by=0.1),mean=.5,log=TRUE)
weights=c(1,1,10,rep(1,8))
weights=c(rep(0,10),1)
weights=weights/sum(weights)
expected.wait=apply(waits,1,weighted.mean,w=weights)
expected.receive=apply(focalreceives,1,weighted.mean,weights=weights)

layout(matrix(1:2,nrow=1))
plot(threshvals,expected.waits,xlab="Threshold",ylab="Expected time to first signal",type='l')
plot(threshvals,expected.receive,xlab="Threshold",ylab="Probability that the focal animal signals first",type='l')
