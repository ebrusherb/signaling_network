setwd("/Users/eleanorbrush/Dropbox/signaling_network")
source("gillespie.R")

tau=1
e=5
bs=2
bf=1
rf=1
dominance=.75
leak=0
ramp=.5
thresh=1

xmax=5
possible=seq(-xmax,xmax,by=1)
s=length(possible)

states=array(,c(s^2,2))
colnames(states)=c("1 state","2 state")
rates=array(,c(s^2,4))
colnames(rates)=c("1 signals","2 signals","1 wins","2 wins")
transitions=array(0,c(s^2,s^2))

v<-possible
states[,1]=as.vector(sapply(v, function (x) rep(x,s)))
states[,2]=rep(v,s)
rates[,1]=apply(as.matrix(states[,1],ncol=1),1,sig_rate,ramp=ramp,thresh=thresh)
rates[,2]=apply(as.matrix(states[,2],ncol=1),1,sig_rate,ramp=ramp,thresh=thresh)
rates[,3]=rf*dominance*(1-rates[,1])*(1-rates[,2])
rates[,4]=rf*(1-dominance)*(1-rates[,1])*(1-rates[,2])
rates=rates

for(i in 1:s^2){
	state=states[i,]
	leakedstate=(1-leak*tau)*state
	newi=which(apply(states==leakedstate,1,sum)==2)
	probs=cbind(dpois(0:s,rates[i,1]),dpois(0:s,rates[i,2]),dpois(0:s,rates[i,3]),dpois(0:s,rates[i,4]))
	for(n1 in 1:s){
		for(n2 in 1:s){
			for(n3 in 1:s){
				for(n4 in 1:s){
					p=probs[n1+1,1]*probs[n2+1,2]*probs[n3+1,3]*probs[n4+1,4]
					change=c((n2)*bs+(n3)*bf-(n4)*bf,(n1)*bs-(n3)*bf+(n4)*bf)
					if(round(p,3)>0){
					print(c(round(p,3),n1,n2,n3,n4,change))}
					changedstate=state+change
					if(changedstate[1]>xmax){changedstate[1]=xmax}
					if(changedstate[2]>xmax){changedstate[2]=xmax}
					if(changedstate[1]<(-xmax)){changedstate[1]=-xmax}
					if(changedstate[2]<(-xmax)){changedstate[2]=-xmax}
					newnewi=which(apply(states==changedstate,1,sum)==2)
						transitions[i,newnewi]=transitions[i,newnewi]+p
				}
			}
		}
	}
}

values=eigen(t(transitions))$values
vectors=eigen(t(transitions))$vectors
stationary=vectors[,1]/sum(vectors[,1])
stationary=Re(stationary)

a1=factor(states[,1])
b1=split(stationary,a1)
p1=sapply(b1,sum)

a2=factor(states[,2])
b2=split(stationary,a2)
p2=sapply(b2,sum)

expectedrates=c(sum(stationary*rates[,1]),sum(stationary*rates[,2]),sum(stationary*rates[,3]),sum(stationary*rates[,4]))



