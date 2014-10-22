re1<-1 #rate of error in animal 1
re2<-re1 #rate of error in animal 2
me1<-0 #mean error in animal 1
me2<-me1 #mean error in animal 2
se1<-1 #variance of error in animal 1
se2<-se1 #variance of error in animal 2

rf<-1 #rate of fighting
mf<-0.1 #mean fight outcome
sf<-0.0 #variance of fight outcome

bs1=0 #boost to animal 2 when animal 1 signals
bs2=0 #boost to animal 1 when animal 2 signals

rs<-1 #maximum signaling rate

T<-2.5  #threshold

erf<-function(x){
	return(2*pnorm(x*sqrt(2))-1)
	}

f<-function(x){
	return(rs*(-erf((x-(-T)))+1)/2)
	}

#for drift diffusion model (size of error does not depend on state)	
error_rate<-function(re1,re2,me1,me2,se1,se2,rf,mf){
	z=T/(me1+mf)
	a=(me1+mf)^2/(se1^2+sf^2)
	x=1/(1+exp(2*z*a)) #from psych.pdf
#	x=1-(1-exp(-2*(me1+mf)*(se1^2+sf^2)))/(1-exp(-4*(me1+mf)*(se1^2+sf^2))) #from stoch calc class
	return(x)
	}
	
dec_time<-function(re1,re2,me1,me2,se1,se2,rf,mf){
	z=T/(me1+mf)
	a=(me1+mf)^2/(se1^2+sf^2)
	x=z*tanh(z*a)
	return(x)
	}
	
#for ornstein-uhlenbeck model (size of error DOES depend on state)
#error_rate<-function(re1,re2,me1,me2,se1,se2,rf,mf){
#	a=(me1+mf)^2/(se1^2+sf^2)
	z=T/(me1+mf)
#	}

timesteps=100 #how many reactions to have
runs=500 #how many separate runs to do
bins=seq(0,timesteps,by=0.5) #bin time so we can compare evolution in time between runs

B<-function(re1,re2,me1,me2,se1,se2,rf,mf,sf,bs1,bs2){ #Gillespie algorithm


#fight_rates=array(,runs) #time spent fighting divided by total time each run
#error1_rates=array(,runs) #rates of errors for each run
#error2_rates=array(,runs)
#signal1_rates=array(,runs)
#signal2_rates=array(,runs)
#dom1=array(,runs) #sum of intervals in which animal 1 is dominant
#dom2=array(,runs)
#hitting1=NULL
#hitting2=NULL
decision1=NULL
decision2=NULL
hitting=NULL
decision=NULL
switches=NULL

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
	#rates=c(re1,0,rf,0,0) #the only processes that happen are animal 1 makes errors and fights--- to compare to 1d stochastic process where there's only one variable being changed
	R=sum(rates)
	r1=runif(1,0,1) #random variable to increment time
	r2=runif(1,0,1)*R #random variable to decide reaction
	tau=1/R*log(1/r1) #exponentially distributed time step
	time[i]=tau+time[i-1] #add to total time
	C=cumsum(rates) 
	j=min(which(C>r2))
	rxns[i-1]=j
	if(j==1){ #animal 1 makes an error
		E=rnorm(1,me1,se1) #ddm variant
		#E=rnorm(1,me1*X[1,i-1],se1) #ornstein-uhlenbeck variant
		X[1,i]=X[1,i]+E
		}
	if(j==2){ #animal 2 makes an error
		E=rnorm(1,me2,se2) #ddm variant
		#E=rnorm(1,me2*X[2,i-1],se2) #ornstein-uhlenbeck variant
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

	
#fight_rates[k]=length(which(rxns==3))/time[timesteps+1]
#error1_rates[k]=length(which(rxns==1))/time[timesteps+1]
#error2_rates[k]=length(which(rxns==2))/time[timesteps+1]
#signal1_rates[k]=length(which(rxns==4))/time[timesteps+1]
#signal2_rates[k]=length(which(rxns==5))/time[timesteps+1]
#dom1[k]=sum(diff(time)[which(X[1,2:timesteps+1]>0)])/time[timesteps+1]
#dom2[k]=sum(diff(time)[which(X[2,2:timesteps+1]>0)])/time[timesteps+1]
#if(sum(X[1,]>-T)==timesteps+1){hitting1[k]=time[timesteps+1]} else hitting1[k]=time[min(which(X[1,]<=-T))]
#if(sum(X[2,]>-T)==timesteps+1){hitting2[k]=time[timesteps+1]} else hitting2[k]=time[min(which(X[2,]<=-T))]
#if(sum(X[1,]>-T)<timesteps+1){hitting1=c(hitting1,time[min(which(X[1,]<=-T))])}
#if(sum(X[2,]>-T)<timesteps+1){hitting2=c(hitting2,time[min(which(X[2,]<=-T))])}
if(sum(abs(X[1,])<T)<timesteps+1){
dt=min(which(abs(X[1,])>=T))
#hitting1=c(hitting1,time[dt])
decision1=c(decision1,X[1,dt])
	}
if(sum(abs(X[2,])<T)<timesteps+1){
dt=min(which(abs(X[2,])>=T))
#hitting2=c(hitting2,time[dt])
decision2=c(decision2,(X[2,dt]))
}
v=apply(X,2,min)
if(sum(v>(-T))<length(v)){
dt=min(which(v<=(-T)))
hitting=c(hitting,time[dt])
decision=c(decision,X[which.min(X[,dt]),dt])} 
v=apply(X,2,prod)
w=which(v<0)
count=0
s=sign(X[1,w[1]])
if(length(w)>1){
for(i in 2:length(w)){
	if(sign(X[1,w[i]])==-sign(X[1,w[i-1]])){
		count=count+1
		}
	s=sign(X[1,w[i]])
	}
	}
switches=c(switches,count)

#layout(1)	
#plot(time,X[1,],type="l",ylim=range(X)) #plot dominance estimates over time
#lines(time,X[2,],col="red")
#abline(h=-T)
#cols=rainbow(5)
#for(j in 5:1){
#	w=which(rxns==j)
#	points(time[w+1],X[1,w+1],col=cols[j],pch=20)
#	points(time[w+1],X[2,w+1],col=cols[j],pch=20)
#	}
#legend(10,max(X[1,])-500,legend=c("Animal 1 Error","Animal 2 Error","Fight","Animal 1 Signal","Animal 2 Signal"),col=cols,lty=rep(1,5),bty="n")

} #end of k loop

for(i in 1:length(bins)){ #get rid of first null value of each list
	where1[[i]]=where1[[i]][-1]
	where2[[i]]=where2[[i]][-1]
	}	

#layout(matrix(1:2,nrow=1))
#boxplot(list(fight_rates,error1_rates,error2_rates,signal1_rates,signal2_rates),names=c("Fights","Animal 1 Error","Animal 2 Error","Animal 1 Signals","Animal 2 Signals"),main="Rates")
#boxplot(list(dom1,dom2),names=c(1,2),main="Time Estimate is Positive")

#to_return=list(where1,where2,hitting1,hitting2,decision1,decision2)
to_return=list(where1,where2,hitting,decision1,decision2)

##compare observed statistics from simulations to what's predicted:
#print(paste("Predicted error rate is ",error_rate(re1,re2,me1,me2,se1,se2,rf,mf)))
print(paste("Observed error rate is ",1-length(which(decision1>0))/length(decision1)))
#print(paste("Expected decision time is ",dec_time(re1,re2,me1,me2,se1,se2,rf,mf)))
print(paste("Mean time until one starts signaling is ",mean(hitting)))
print(paste("On average ",mean(switches)," switches occurred"))

return(to_return)
}

T<-5  #threshold

me1<--2 #mean error in animal 1
me2<-me1 #mean error in animal 2
se1<-3 #variance of error in animal 1
se2<-se1 #variance of error in animal 2

mf<-1 #mean fight outcome
sf<-.5 #variance of fight outcome

bs1=10 #boost to animal 2 when animal 1 signals
bs2=bs1 #boost to animal 1 when animal 2 signals

runs=100

L=B(re1,re2,me1,me2,se1,se2,rf,mf,sf,bs1,bs2)

#l=which(bins==ceiling(max(bins)/3))
l=which(bins==30)
layout(matrix(1:2,nrow=1))
boxplot(L[[1]][1:l],axes=FALSE)
axis(1,at=1:l,bins[1:l])
axis(2)
boxplot(L[[2]][1:l],axes=FALSE)
axis(1,at=1:l,bins[1:l])
axis(2)	
