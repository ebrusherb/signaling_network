#needed to run:
erf<-function(x){
	return(2*pnorm(x*sqrt(2))-1)
	}

f<-function(x){
	rs=1
	return(rs*(-erf(10*(x-(-T)))+1)/2)
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


################
################
################
B<-function(re1,re2,me1,me2,se1,se2,rf,mf,sf,bs1,bs2){ 

rs=1 # max signaling rate

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
	X[1,i]=X[1,i-1]+(me1*X[1,i-1]+mf+bs2*f(X[2,i-1]))*deltat+sf*sqrt(deltat)*z[1]+se1*sqrt(deltat)*z[2]#+bs2*sqrt(f(X[2,i-1]))*z[3]
	X[2,i]=X[2,i-1]+(me2*X[2,i-1]-mf+bs1*f(X[1,i-1]))*deltat-sf*sqrt(deltat)*z[1]+se2*sqrt(deltat)*z[4]#+bs1*sqrt(f(X[1,i-1]))*z[5]
	time[i]=time[i-1]+deltat
#	print(c(sf*sqrt(deltat)*z[1],X[1,i-1],X[1,i]))
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
if(length(w)>1){
for(i in 2:length(w)){
	if(which.max(X[,w[i]])!=which.max(X[,w[i-1]])){
		count=count+1
		}
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
#################################
#################################
#################################

T<-.5  #threshold

me1<--.5 #mean error in animal 1
me2<-me1 #mean error in animal 2
se1<-.5 #variance of error in animal 1
se2<-se1 #variance of error in animal 2

mf<-0 #mean fight outcome
sf<-.5 #variance of fight outcome

bs1=0 #boost to animal 2 when animal 1 signals
bs2=bs1 #boost to animal 1 when animal 2 signals

timesteps=1000 #how many reactions to have
runs=100 #how many separate runs to do
deltat=0.05
bins=seq(0,timesteps*deltat,by=deltat) #bin time so we can compare evolution in time between runs

L=B(re1,re2,me1,me2,se1,se2,rf,mf,sf,0,0)
Lb=B(re1,re2,me1,me2,se1,se2,rf,mf,sf,3,3)

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
#abline(v=-mf/me1,col="red")
legend(20,0.2,legend=c("w/o Feedback","w/ Feedback"),fill=c(rgb(0,0,1,1),rgb(0,0,1,.25)),bty="n")

hist(L[[2]][[k]],col=rgb(0,0,1,alph=1),freq=FALSE,xlim=range(c(L[[2]][[k]],Lb[[2]][[k]]+5)),main=paste("Animal 2's Estimate at t=",bins[k]),xlab="Estimate",breaks=10)
hist(Lb[[2]][[k]],col=rgb(0,0,1,alph=.25),freq=FALSE,add=T,breaks=10)
#abline(v=mf/me1,col="red")

#to determine the joint effects of mean fight outcome and signaling feedback on the number of switches in dominance regime:

T<-.5  #threshold

me1<--.5 #mean error in animal 1
me2<-me1 #mean error in animal 2
se1<-.5 #variance of error in animal 1
se2<-se1 #variance of error in animal 2

mf<-0 #mean fight outcome
sf<-.5 #variance of fight outcome

bs1=0 #boost to animal 2 when animal 1 signals
bs2=bs1 #boost to animal 1 when animal 2 signals

timesteps=1000 #how many reactions to have
runs=100 #how many separate runs to do
deltat=0.05
bins=seq(0,timesteps*deltat,by=deltat) #bin time so we can compare evolution in time between runs

fight_diff<-seq(0,5,by=.5)
fight_var<-seq(0,10,by=1)
signal_feedback<-seq(0,20,by=1)
switch_mat=array(,c(length(fight_diff),length(fight_var),length(signal_feedback)))
error_mat=array(,c(length(fight_diff),length(fight_var),length(signal_feedback)))
time_mat=array(,c(length(fight_diff),length(fight_var),length(signal_feedback)))
for(i in 1:length(fight_diff)){
	for(k in 1:length(fight_var)){
	for(j in 1:length(signal_feedback)){
	#for(j in 1:1){
		L=B(re1,re2,me1,me2,se1,se2,rf,fight_diff[i],fight_var[k],signal_feedback[j],signal_feedback[j])
		switch_mat[i,k,j]=mean(L[[5]])
		error_mat[i,k,j]=length(which(L[[3]]==1))/length(L[[3]])
		time_mat[i,k,j]=mean(L[[4]])
		}
		}
	}
#image(A,axes=FALSE,xlab="Difference in Fighting Ability",ylab="Signaling Feedback",main=paste("Numbers of switches, me=",me1,"se=",se1,"sf=",sf))
#axis(1,at=seq(0,1,length.out=6),seq(min(fight_diff),max(fight_diff),length.out=6))
#axis(2,at=seq(0,1,length.out=6),seq(min(signal_feedback),max(signal_feedback),length.out=6))
library(fields)
layout(matrix(1:3,nrow=1))
k=5
image.plot(fight_diff,signal_feedback,switch_mat[,k,],xlab="Difference in Fighting Ability",ylab="Signaling Feedback",main=paste("Switches, me=",me1,"se=",se1,"sf=",fight_var[k]))
image.plot(fight_diff,signal_feedback,error_mat[,k,],xlab="Difference in Fighting Ability",ylab="Signaling Feedback",main=paste("Incorrect Signaling, me=",me1,"se=",se1,"sf=",fight_var[k]))
image.plot(fight_diff,signal_feedback,time_mat[,k,],xlab="Difference in Fighting Ability",ylab="Signaling Feedback",main=paste("Time to Signaling, me=",me1,"se=",se1,"sf=",fight_var[k]))

k=2
image.plot(fight_diff,fight_var,switch_mat[,,k],xlab="Difference in Fighting Ability",ylab="Variance in Fight Outcome",main=paste("Switches, me=",me1,"se=",se1,"b=",signal_feedback[k]))
image.plot(fight_diff,fight_var,error_mat[,,k],xlab="Difference in Fighting Ability",ylab="Variance in Fight Outcome",main=paste("Incorrect Signaling, me=",me1,"se=",se1,"b=",signal_feedback[k]))
image.plot(fight_diff,fight_var,time_mat[,,k],xlab="Difference in Fighting Ability",ylab="Variance in Fight Outcome",main=paste("Time to Signaling, me=",me1,"se=",se1,"b=",signal_feedback[k]))

k=1
image.plot(fight_var,signal_feedback,switch_mat[k,,],xlab="Variance in Fight Outcome",ylab="Signaling Feedback",main=paste("Switches, me=",me1,"se=",se1,"mf=",fight_diff[k]))
image.plot(fight_var,signal_feedback,error_mat[k,,],xlab="Variance in Fight Outcome",ylab="Signaling Feedback",main=paste("Incorrect Signaling, me=",me1,"se=",se1,"mf=",fight_diff[k]))
image.plot(fight_var,signal_feedback,time_mat[k,,],xlab="Variance in Fight Outcome",ylab="Signaling Feedback",main=paste("Time to Signaling, me=",me1,"se=",se1,"mf=",fight_diff[k]))

#embedding it in a network:

f<-function(x){
	y=ifelse(x<=-T,1,0)
	return(y)
	}

T<-.5  #threshold

me1<--.5 #mean error in animal 1
me2<-me1 #mean error in animal 2
se1<-.5 #variance of error in animal 1
se2<-se1 #variance of error in animal 2

sf<-.5 #variance of fight outcome

bs1=0 #boost to animal 2 when animal 1 signals
bs2=bs1 #boost to animal 1 when animal 2 signals

n=10
abilities=sort(runif(20,0,10),decreasing=TRUE)
runs=1
timesteps=1000 #how many reactions to have
deltat=0.05
time=seq(0,timesteps*deltat,by=deltat)

hitting=array(0,c(runs,n,n))
decision=array(0,c(runs,n,n))
switches=array(0,c(runs,n,n))

X=array(0,c(runs,n,n,timesteps+1))
Snow=array(0,c(runs,n,n,timesteps+1))
Scum=array(0,c(runs,n,n,timesteps+1))

for(k in 1:runs){

X[k,,,1]=0
S[k,,,1]=0

for(i in 1:timesteps+1){
	for(m in 1:(n-1)){
		for(p in (m+1):n){
	z=rnorm(5,0,1)
	X[k,m,p,i]=X[k,m,p,i-1]+(me1*X[k,m,p,i-1]+abilities[m]-abilities[p]+bs2*f(X[k,m,p,i-1]))*deltat+sf*sqrt(deltat)*z[1]+se1*sqrt(deltat)*z[2]#+bs2*sqrt(f(X[2,i-1]))*z[3]
	X[k,p,m,i]=X[k,p,m,i-1]+(me2*X[k,p,m,i-1]+abilities[p]-abilities[m]+bs1*f(X[k,p,m,i-1]))*deltat+sf*sqrt(deltat)*z[1]+se2*sqrt(deltat)*z[4]#+bs1*sqrt(f(X[1,i-1]))*z[5]
	Snow[k,m,p,i]=f(X[k,m,p,i])
	Snow[k,p,m,i]=f(X[k,p,m,i])
	Scum[k,m,p,i]=Scum[k,m,p,i-1]+Snow[k,m,p,i]
	Scum[k,p,m,i]=Scum[k,p,m,i-1]+Snow[k,p,m,i]
		}
		}
	}

for(m in 1:(n-1)){
	for(p in (m+1):n){
	pair<-rbind(X[k,m,p,],X[k,p,m,])
	v=apply(pair,2,min)
	dt=min(which(v<=(-T)))
	hitting[k,m,p]=time[dt]
	decision[k,m,p]=which.min(pair[,dt])
	u=apply(pair,2,max)
	w=intersect(which(v<(-T)),which(u>T))
	count=0
if(length(w)>1){
for(i in 2:length(w)){
	if(which.max(X[k,,,w[i]])!=which.max(X[k,,,w[i-1]])){
		count=count+1
		}
	}
	}
switches[k,m,p]=count
			}
	}
 

} #end of k loop

layout(matrix(1:9,nrow=3,byrow=T))

toplot<-seq(100,timesteps,by=20)
cols=heat.colors(length(toplot))

for(i in 1:9){
	plot(metrics(i,Scum[1,,,1000])/max(metrics(i,Scum[1,,,1000])),type="l",xlab="",ylab="",main=metric_names(i))
	for(j in 1:length(toplot)){
		lines(metrics(i,Scum[1,,,toplot[j]])/max(metrics(i,Scum[1,,,toplot[j]])),col=cols[j])
		}
	}
	
for(i in 1:9){
	hist(metrics(i,Scum[1,,,1000]),xlab="",main=metric_names(i))
	}
