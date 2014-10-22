re1<-1 #rate of error in animal 1
re2<-re1 #rate of error in animal 2
me1<-0 #mean error in animal 1
me2<-0 #mean error in animal 2
se1<-100 #variance of error in animal 1
se2<-100 #variance of error in animal 2

rf<-5 #rate of fighting
mf<-5 #mean fight outcome
sf<-0 #variance of fight outcome

bs1=0 #boost to animal 2 when animal 1 signals
bs2=0 #boost to animal 1 when animal 2 signals

rs<-10

T<-25 #signaling threshold

erf<-function(x){
	return(2*pnorm(x*sqrt(2))-1)
	}
f<-function(x){
	return(rs*(erf(-(x-2*(-T)))+1))
	}

timesteps=1000 #how many reactions to have
bins=seq(0,timesteps,by=0.5) #bin time so we can compare evolution in time between runs

B<-function(re1,re2,me1,me2,se1,se2,rf,mf,sf,bs1,bs2){ #Gillespie algorithm

runs=500 #how many separate runs to do
fight_rates=array(,runs) #time spent fighting divided by total time each run
error1_rates=array(,runs) #rates of errors for each run
error2_rates=array(,runs)
signal1_rates=array(,runs)
signal2_rates=array(,runs)
dom1=array(,runs) #sum of intervals in which animal 1 is dominant
dom2=array(,runs)
hitting2=array(,runs)



where1=list() #each element of the list will be an array of dominance values for animal 1 at a point in time
where2=list()
for(i in 1:length(bins)){
	where1[[i]]=array()
	where2[[i]]=array()
	}

for(k in 1:runs){
time=array(,timesteps+1) #keep track of timepoints
time[1]=0 #initialize
X=array(,c(2,timesteps+1))
X[,1]=c(0,0)
rxns=array(,timesteps) #keep track of which reactions occur

for(i in 1:timesteps+1){
	X[,i]=X[,i-1]
	rates=c(re1,re2,rf,f(X[1,i-1]),f(X[2,i-1]))
	R=sum(rates)
	r1=runif(1,0,1) #random variable to increment time
	r2=runif(1,0,1)*R #random variable to decide reaction
	tau=1/R*log(1/r1) #exponentially distributed time step
	time[i]=tau+time[i-1] #add to total time
	C=cumsum(rates) 
	j=min(which(C>r2))
	rxns[i-1]=j
	if(j==1){ #animal 1 makes an error
		E=rnorm(1,me1,se1)
		X[1,i]=X[1,i]+E
		}
	if(j==2){ #animal 2 makes an error
		E=rnorm(1,me2,se2)
		X[2,i]=X[2,i]+E
		}
	if(j==3){ #fight occurs
		F=rnorm(1,mf,sf)
		X[1,i]=X[1,i]+F
		X[2,i]=X[2,i]-F
		}
	if(j==4){ #animal 2 signals
		X[2,i]=X[2,i]+bs1
		}
	if(j==5){ #animal 1 signals
		X[1,i]=X[1,i]+bs2
		}	

	}
	
for(i in 1:length(bins)){ #for each bin concatenate the dominance values at the beginning of that time bin
	if(bins[i]<=max(time)){
		M=max(which(time<=bins[i]))
		where1[[i]]=c(where1[[i]],X[1,M])
		where2[[i]]=c(where2[[i]],X[2,M])
		}

	}
	
	
fight_rates[k]=length(which(rxns==3))/time[timesteps+1]
error1_rates[k]=length(which(rxns==1))/time[timesteps+1]
error2_rates[k]=length(which(rxns==2))/time[timesteps+1]
signal1_rates[k]=length(which(rxns==4))/time[timesteps+1]
signal2_rates[k]=length(which(rxns==5))/time[timesteps+1]
dom1[k]=sum(diff(time)[which(X[1,2:timesteps+1]>0)])/time[timesteps+1]
dom2[k]=sum(diff(time)[which(X[2,2:timesteps+1]>0)])/time[timesteps+1]
if(sum(X[1,]>-T)==timesteps+1){hitting2[k]=time[timesteps+1]} else hitting2[k]=time[min(which(X[1,]<=-T))]
	
plot(time,X[1,],type="l",ylim=range(X)) #plot dominance estimates over time
lines(time,X[2,],col="red")
abline(h=-T)

} #end of k loop

for(i in 1:length(bins)){ #get rid of first null value of each list
	where1[[i]]=where1[[i]][-1]
	where2[[i]]=where2[[i]][-1]
	}


layout(matrix(1:2,nrow=1))
boxplot(list(fight_rates,error1_rates,error2_rates,signal1_rates,signal2_rates),names=c("Fights","Animal 1 Error","Animal 2 Error","Animal 1 Signals","Animal 2 Signals"),main="Rates")
boxplot(list(dom1,dom2),names=c(1,2),main="Time Estimate is Positive")

to_return=list(where1,where2,hitting2)

return(to_return)
}

#for(i in 1:timesteps){
#	if(rxns[i]==3){abline(v=i)}
#	}


L=B(re1,re2,me1,me2,se1,se2,rf,mf,sf,bs1,bs2)

layout(1)

where1=L[[1]]
where2=L[[2]]
M=50
pvals1<-array(,M)
for(k in 1:M){
	v=seq(min(where1[[k]]),max(where1[[k]]),by=.5)
	v=as.matrix(v,nrow=1)
	t=apply(v,2,dnorm,mean=(me2-mf)*bins[k],sd=(se1+sf)*sqrt(bins[k])) #generate 
	hist(where2[[k]],freq=FALSE,main=paste(k,bins[k]))
	lines(v,t)
	pvals1[k]=ks.test(where2[[k]],"pnorm",(me1+mf)*bins[k],se1*sqrt(bins[k]))$p.value
	}
hist(pvals1)

hit<-L[[3]]
tvec=seq(min(hit),max(hit),by=0.02)
tvec=matrix(tvec,nrow=1)
g<-function(t,C){
	return(abs(C)/((se1+sf)*sqrt(2*pi*t^3))*exp(-(C-(me1-mf)*t)^2/(2*(se1+sf)^2*t)))
	}
gvec=apply(tvec,2,g,-T)
