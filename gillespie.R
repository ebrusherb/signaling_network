sig_rate<-function(x,ramp=1,thresh=1){
	1/(1+exp(4*ramp*(x+thresh)))
}

sims<-function(e,maxtime,ramp,thresh1,thresh2,rf,dominance,leak,bs,bf){

	its=100
	boost=matrix(c(0,bs,bs,0,bf,-bf,-bf,bf,0,0),nrow=2)
	X=list()
	TIMES=list()
	RXNS=list()

	for(j in 1:its){
		x=array(0,c(2,1))
		x.now=x
		times=0
		time.now=times
		rxns=array(0,5)
		while(time.now<maxtime){
			rxn.now=array(0,5)
			a0=e
			s1=sig_rate(x.now[1],ramp=ramp,thresh=thresh1)
			s2=sig_rate(x.now[2],ramp=ramp,thresh=thresh2)
			s=s1+s2-s1*s2
			rates=c(s1,s2,(1-s)*rf*dominance,(1-s)*rf*(1-dominance))
			r=runif(2,0,1)
			tau=1/a0*log(1/r[1])
			rxn.now[min(which(r[2]<=c(cumsum(rates),1)))]=1
			deltax=boost%*%rxn.now
			x.now=(1-leak*tau)*x.now+deltax
			time.now=time.now+tau
			x=cbind(x,x.now)
			times=c(times,time.now)
			rxns=cbind(rxns,rxn.now)
		}
	X[[j]]=cbind(x,x.now)
	TIMES[[j]]=c(times,maxtime)
	RXNS[[j]]=cbind(rxns,rep(1,5))
	dimnames(RXNS[[j]])[[2]]=NULL
	}	
return(list(estimates=X,times=TIMES,rxns=RXNS))
}

stats<-function(simulation){
	estimates=simulation$estimates
	times=simulation$times
	rxns=simulation$rxns
	L=length(times) #number of simulations 
	num.rxns=array(,L)
	waits=array(,L)
	firstsignaler=array(0,c(2,L))
	switches.count=array(0,L)
	switches.times=list()
	transitions=list()
	regimes=NULL
	for(i in 1:L){
		num.rxns[i]=length(times[[i]])-2
		
		signaling=rxns[[i]][1:2,]
		signaling.happened=which(apply(signaling,2,sum)>0)
		signaling=signaling[,signaling.happened]
		signaling=matrix(signaling,nrow=2)
		
		t=times[[i]][signaling.happened]
		waits[i]=t[1]
		firstsignaler[,i]=signaling[,1]
		
		state=array(0,c(2,length(times[[i]])))
		state[1,][which(estimates[[i]][1,]<(-thresh1))]=1
		state[2,][which(estimates[[i]][2,]<(-thresh2))]=2
		state=apply(state,2,sum)
		d=diff(state)
		regimes=rbind(regimes,cbind(c(state[d!=0],state[length(state)]),diff(c(0,times[[i]][d!=0],times[[i]][length(times[[i]])]))),c(4,0))
		
		# if(length(signaling.happened)>1){
			# d=apply(signaling,1,diff) #find changes in signaling 
			# d=matrix(d,ncol=2)
			# numswitches=sum(apply(d,1,prod)==-1) #number of changes in signaling
			# switches.count[i]=numswitches
			# if(numswitches>0){
				# switches.times[[i]]=t[which(apply(d,1,prod)==-1)+1]
			# }
		# }
	}
	weaker.signaled=intersect(which(firstsignaler[2,]==1),which(firstsignaler[1,]==0))
	stronger.signaled=intersect(which(firstsignaler[2,]==0),which(firstsignaler[1,]==1))
	signaling.happened=union(weaker.signaled,stronger.signaled)
	exp.wait=mean(waits[signaling.happened])
	prob.weaker.sigs=length(weaker.signaled)/length(signaling.happened)
	
	m=matrix(0,nrow=4,ncol=4)
	M=matrix(0,nrow=4,ncol=4)
		for(j in 1:4){
			w=which(regimes[,1]==j-1)
			for(k in 1:4){
				m[j,k]=length(which(regimes[,1][w+1]==(k-1)))/length(which(regimes[,1][w+1]!=4))
				M[j,k]=mean(regimes[,2][w][which(regimes[,1][w+1]==(k-1))])
			}
		}
		transition.probs=m
		transition.times=M
	return(list(num.rxns=num.rxns,num.fights=num.fights,num.sig1=num.sig1,num.sig2=num.sig2,waits=waits,firstsignaler=firstsignaler,exp.wait=exp.wait,prob.weaker.sigs=prob.weaker.sigs,trans.probs=transition.probs,trans.times=transition.times))
}