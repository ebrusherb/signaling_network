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
error_rate<-function(T,me1,se1,mf,sf){
	z=T/(me1+mf)
	a=(me1+mf)^2/(se1^2+sf^2)
	x=1/(1+exp(2*z*a)) #from psych.pdf
#	x=1-(1-exp(-2*(me1+mf)*(se1^2+sf^2)))/(1-exp(-4*(me1+mf)*(se1^2+sf^2))) #from stoch calc class
	return(x)
	}
	
dec_time<-function(T,me1,se1,mf,sf){
	z=T/(me1+mf)
	a=(me1+mf)^2/(se1^2+sf^2)
	x=z*tanh(z*a)
	return(x)
	}
	

timesteps=100 #how many reactions to have
runs=500 #how many separate runs to do
bins=seq(0,timesteps,by=0.5) #bin time so we can compare evolution in time between runs

B<-function(re1,re2,me1,me2,se1,se2,rf,mf,sf,bs1,bs2){ #Gillespie algorithm


hitting=NULL
decision=NULL

where1=list() #each element of the list will be an array of dominance values for animal 1 at a point in time
for(i in 1:length(bins)){
	where1[[i]]=array()
	}


for(k in 1:runs){
#time=array(,timesteps+1) #keep track of timepoints
#time[1]=0 #initialize
#X=array(,c(timesteps+1))
#X[1]=0
#rxns=array(,timesteps) #keep track of which reactions occur
#
#for(i in 1:timesteps+1){
#	X[i]=X[i-1]
#	rates=c(re1,rf)
#	R=sum(rates)
#	r1=runif(1,0,1) #random variable to increment time
#	r2=runif(1,0,1)*R #random variable to decide reaction
#	tau=1/R*log(1/r1) #exponentially distributed time step
#	time[i]=tau+time[i-1] #add to total time
#	C=cumsum(rates) 
#	j=min(which(C>r2))
#	rxns[i-1]=j
#	if(j==1){ #animal 1 makes an error
#		E=rnorm(1,me1,se1) #ddm variant
#		#E=rnorm(1,me1*X[i-1],se1) #ornstein-uhlenbeck variant
#		X[i]=X[i]+E
#		}
#	if(j==2){ #fight occurs
#		F=rnorm(1,mf,sf)
#		X[i]=X[i]+F
#		}
#	}

time=array(,timesteps+1) #keep track of timepoints
time[1]=0 #initialize
X=array(,c(timesteps+1))
X[1]=0

deltat=0.001
for(i in 1:timesteps+1){
	z1=rnorm(1,0,1)
	z2=rnorm(1,0,1)
	X[i]=X[i-1]+(me1+mf)*deltat+se1*sqrt(deltat)*z1+sf*sqrt(deltat)*z2
	time[i]=time[i-1]+deltat
	}

	
for(i in 1:length(bins)){ #for each bin concatenate the dominance values at the beginning of that time bin
	if(bins[i]<=max(time)){
		M=max(which(time<=bins[i]))
		where1[[i]]=c(where1[[i]],X[M])
		}

	}

	

if(sum(abs(X)<(T))<length(X)){
dt=min(which(abs(X)>(T)))
hitting=c(hitting,time[dt])
decision=c(decision,sign(X[dt]))} 

layout(1)	
plot(time,X[],type="l",ylim=range(X)) #plot dominance estimates over time
abline(h=-T)
abline(h=T)
cols=rainbow(2)
#for(j in 2:1){
#	w=which(rxns==j)
#	points(time[w+1],X[w+1],col=cols[j],pch=20)
#	}
abline(v=time[dt])
#legend(10,max(X[1,])-500,legend=c("Animal 1 Error","Animal 2 Error","Fight","Animal 1 Signal","Animal 2 Signal"),col=cols,lty=rep(1,5),bty="n")

} #end of k loop

for(i in 1:length(bins)){ #get rid of first null value of each list
	where1[[i]]=where1[[i]][-1]
	}	

#layout(matrix(1:2,nrow=1))
#boxplot(list(fight_rates,error1_rates,error2_rates,signal1_rates,signal2_rates),names=c("Fights","Animal 1 Error","Animal 2 Error","Animal 1 Signals","Animal 2 Signals"),main="Rates")
#boxplot(list(dom1,dom2),names=c(1,2),main="Time Estimate is Positive")

#to_return=list(where1,where2,hitting1,hitting2,decision1,decision2)
to_return=list(where1,decision,hitting)

##compare observed statistics from simulations to what's predicted:
#print(paste("Predicted error rate is ",error_rate(re1,re2,me1,me2,se1,se2,rf,mf)))
print(paste("Observed error rate is ",length(which(decision==-1))/length(decision)))
#print(paste("Expected decision time is ",dec_time(re1,re2,me1,me2,se1,se2,rf,mf)))
print(paste("Mean time until one starts signaling is ",mean(hitting)))

return(to_return)
}

T<-2  #threshold

me1<--.5 #mean error in animal 1
me2<-me1 #mean error in animal 2
se1<-1 #variance of error in animal 1
se2<-se1 #variance of error in animal 2

mf<-1 #mean fight outcome
sf<-2 #variance of fight outcome

bs1=0 #boost to animal 2 when animal 1 signals
bs2=bs1 #boost to animal 1 when animal 2 signals

timesteps=2000
bins=seq(0,timesteps,by=0.5) #bin time so we can compare evolution in time between runs
runs=1000

L=B(re1,re2,me1,me2,se1,se2,rf,mf,sf,0,0)
error_rate(T,me1,se1,mf,sf)
dec_time(T,me1,se1,mf,sf)


#l=which(bins==ceiling(max(bins)/3))
l=which(bins==30)
boxplot(L[[1]][1:l],axes=FALSE)
axis(1,at=1:l,bins[1:l])
axis(2)
abline(h=-mf/me1,col="red")

#k=401
#hist(L[[1]][[k]],col=rgb(1,0,0,alph=.25),freq=FALSE,xlim=range(c(L[[1]][[k]],Lb[[1]][[k]])),main=paste("Animal 1's Estimate at t=",bins[k]),xlab="Estimate",breaks=20)
#hist(Lb[[1]][[k]],col=rgb(0,0,1,alph=.25),freq=FALSE,add=T,breaks=20)
#legend(40,0.08,legend=c("w/o Feedback","w/ Feedback"),fill=c(rgb(1,0,0,.25),rgb(0,0,1,.25)),bty="n")

