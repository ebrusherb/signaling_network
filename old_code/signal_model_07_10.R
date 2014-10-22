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

timesteps=100 #how many reactions to have
runs=500 #how many separate runs to do
deltat=0.01
bins=seq(0,timesteps*deltat,by=deltat) #bin time so we can compare evolution in time between runs

B<-function(re1,re2,me1,me2,se1,se2,rf,mf,sf,bs1,bs2){ #Gillespie algorithm


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

for(i in 1:timesteps+1){
	z=rnorm(5,0,1)
	X[1,i]=X[1,i-1]+(me1*X[1,i-1]+mf+bs2*f(X[2,i-1]))*deltat+sf*sqrt(deltat)*z[1]+se1*sqrt(deltat)*z[2]+bs2*sqrt(f(X[2,i-1]))*z[3]
	X[2,i]=X[2,i-1]+(me1*X[2,i-1]-mf+bs1*f(X[1,i-1]))*deltat+sf*sqrt(deltat)*z[1]+se2*sqrt(deltat)*z[4]+bs1*sqrt(f(X[1,i-1]))*z[5]
	time[i]=time[i-1]+deltat
	}
	
for(i in 1:length(bins)){ #for each bin concatenate the dominance values at the beginning of that time bin
		where1[[i]]=c(where1[[i]],X[1,i])
		where2[[i]]=c(where2[[i]],X[2,i])
	}

v=apply(X,2,min)
if(sum(v>(-T))<length(v)){
dt=min(which(v<=(-T)))
hitting=c(hitting,time[dt])
decision=c(decision,which.min(X[,dt]))} 

u=apply(X,2,max)
w=intersect(which(v<(-T)),which(u>T))
count=0
s=which.max(X[,w[1]])
if(length(w)>1){
for(i in 2:length(w)){
	if(which.max(X[,w[i]])!=which.max(X[,w[i-1]])){
		count=count+1
		}
	s=sign(X[1,w[i]])
	}
	}
switches=c(switches,count)

layout(1)	
plot(time,X[1,],type="l",ylim=range(X)) #plot dominance estimates over time
lines(time,X[2,],col="red")
abline(h=-T)

} #end of k loop

for(i in 1:length(bins)){ #get rid of first null value of each list
	where1[[i]]=where1[[i]][-1]
	where2[[i]]=where2[[i]][-1]
	}	

to_return=list(where1,where2,decision,hitting,switches)

##compare observed statistics from simulations to what's predicted:
#print(paste("Predicted error rate is ",error_rate(re1,re2,me1,me2,se1,se2,rf,mf)))
print(paste("Observed error rate is ",length(which(decision==1))/length(decision)))
#print(paste("Expected decision time is ",dec_time(re1,re2,me1,me2,se1,se2,rf,mf)))
print(paste("Mean time until one starts signaling is ",mean(hitting)))
print(paste("On average ",mean(switches)," switches occurred"))

return(to_return)
}

T<-.5  #threshold

me1<--.5 #mean error in animal 1
me2<-me1 #mean error in animal 2
se1<-1 #variance of error in animal 1
se2<-se1 #variance of error in animal 2

mf<-5 #mean fight outcome
sf<-.5 #variance of fight outcome

bs1=0 #boost to animal 2 when animal 1 signals
bs2=bs1 #boost to animal 1 when animal 2 signals

timesteps=1000 #how many reactions to have
runs=1000 #how many separate runs to do
deltat=0.1
bins=seq(0,timesteps*deltat,by=deltat) #bin time so we can compare evolution in time between runs

L=B(re1,re2,me1,me2,se1,se2,rf,mf,sf,0,0)
Lb=B(re1,re2,me1,me2,se1,se2,rf,mf,sf,bs1,bs2)

#l=which(bins==ceiling(max(bins)/3))
l=which(bins==100)
layout(matrix(1:2,nrow=1))
boxplot(L[[1]][1:l],axes=FALSE)
axis(1,at=1:l,bins[1:l])
axis(2)
abline(h=-mf/me1,col="red")
boxplot(L[[2]][1:l],axes=FALSE)
axis(1,at=1:l,bins[1:l])
axis(2)	
abline(h=mf/me1,col="red")

#to confirm that without signaling, the estimates reach the stationary distributions they should:
k=951
layout(matrix(1:2,nrow=1))
hist(L[[1]][[k]],freq=FALSE,xlab="Estimate",ylab="Probability",main=paste("Animal 1's Estimate at t=",bins[k]))
v=seq(min(L[[1]][[k]]),max(L[[1]][[k]]),by=0.01)
v=matrix(v,nrow=1)
v2=apply(v,2,dnorm,mean=-mf/me1,sd=sqrt(-(sf^2+se1^2)/(2*me1)))
lines(v,v2,col="red")
hist(L[[2]][[k]],freq=FALSE,xlab="Estimate",ylab="Probability",main=paste("Animal 1's Estimate at t=",bins[k]))
v=seq(min(L[[2]][[k]]),max(L[[2]][[k]]),by=0.01)
v=matrix(v,nrow=1)
v2=apply(v,2,dnorm,mean=mf/me1,sd=sqrt(-(sf^2+se1^2)/(2*me1)))
lines(v,v2,col="red")

ks.test(L[[1]][[k]],pnorm,mean=-mf/me1,sd=sqrt(-(sf^2+se1^2)/(2*me1)))
ks.test(L[[2]][[k]],pnorm,mean=mf/me1,sd=sqrt(-(sf^2+se1^2)/(2*me1)))

#to compare the stationary distributions with and without signaling feedback:
k=1001
layout(matrix(1:2,nrow=1))
hist(L[[1]][[k]],col=rgb(0,0,1,alph=1),freq=FALSE,xlim=range(c(L[[1]][[k]],Lb[[1]][[k]]+5)),main=paste("Animal 1's Estimate at t=",bins[k]),xlab="Estimate",breaks=10)
hist(Lb[[1]][[k]],col=rgb(0,0,1,alph=.25),freq=FALSE,add=T,breaks=10)
abline(v=-mf/me1,col="red")
legend(10,0.2,legend=c("w/o Feedback","w/ Feedback"),fill=c(rgb(0,0,1,1),rgb(0,0,1,.25)),bty="n")

hist(L[[2]][[k]],col=rgb(0,0,1,alph=1),freq=FALSE,xlim=range(c(L[[2]][[k]],Lb[[2]][[k]]+5)),main=paste("Animal 2's Estimate at t=",bins[k]),xlab="Estimate",breaks=10)
hist(Lb[[2]][[k]],col=rgb(0,0,1,alph=.25),freq=FALSE,add=T,breaks=10)
abline(v=mf/me1,col="red")

#to determine the joint effects of mean fight outcome and signaling feedback on the number of switches in dominance regime:
fight_diff<-seq(0,10,by=1)
signal_feedback<-seq(0,10,by=1)
A=array(,c(length(fight_diff),length(signal_feedback)))
for(i in 1:length(fight_diff)){
	for(j in 1:length(signal_feedback)){
		L=B(re1,re2,me1,me2,se1,se2,rf,fight_diff[i],sf,signal_feedback[j],signal_feedback[j])
		A[i,j]=mean(L[[5]])
		}
	}
image(A,axes=FALSE,xlab="Difference in Fighting Ability",ylab="Signaling Feedback",main=paste("Numbers of switches, me=",me1,"se=",se1,"sf=",sf))
axis(1,at=seq(0,1,length.out=6),seq(min(fight_diff),max(fight_diff),length.out=6))
axis(2,at=seq(0,1,length.out=6),seq(min(signal_feedback),max(signal_feedback),length.out=6))