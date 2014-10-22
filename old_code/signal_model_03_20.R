lambda=-.01
alpha=0
a1=1
a2=1
sigma2=2
T=-5

steps=100000
runs=1
x1=array(0,c(runs,steps))
x2=array(0,c(runs,steps))
s1=array(0,c(runs,steps))
s2=array(0,c(runs,steps))

for(j in 1:runs){
for(i in 3:steps){
 	if(x2[j,i-2]>T&&x2[j,i-1]<=T){s2[j,i]=1} else s2[j,i]=0
 	if(x1[j,i-2]>T&&x1[j,i-1]<=T){s1[j,i]=1} else s1[j,i]=0
 	fight=rnorm(1,a1-a2,sqrt(sigma2))
	x1[j,i]=(1+lambda)*x1[j,i-1]+fight+alpha*s2[j,i]
	x2[j,i]=(1+lambda)*x2[j,i-1]-fight+alpha*s1[j,i]
	}
}
received1=apply(s2,1,sum)
received2=apply(s1,1,sum)

layout(matrix(1:1,nrow=1))

#plot(x1[1,],type="l",ylim=range(cbind(x1,x2)))
#for(j in 2:runs){lines(x1[j,])}
#lines(x2[1,],col="red")
#for(j in 2:runs){lines(x2[j,],col="red")}

#plot(s1[1,])
#points(s2[1,],col="red")
#hist(received1)
#hist(received2,col="red")
#hist(received1-received2)

both_exchanged=intersect(which(received1!=0),which(received2!=0))

lags=NULL

for(j in both_exchanged){
	
	exchanged=array(0,steps)
	exchanged[which(s1[j,]==1)]=1
	exchanged[which(s2[j,]==1)]=-1
	
	k=min(which(exchanged!=0))
	direction=exchanged[k]
	walked=0
	
	while(k<steps){
		k=k+1
		if(exchanged[k]==0){walked=walked+1}
		if(exchanged[k]==direction){walked=walked+1}
		if(exchanged[k]==-1*direction){
			direction=-1*direction
			lags=c(lags,walked)
			walked=0			
			}
		}
	
	}
	
hist(lags)
print(mean(lags))


##adrian's model
#
#delta=.01
#psi=1
#alpha=50
#lambda=-0.01
#
#runs=100
#max=1000
#
#dv=list()
#dv_leak=list()
#lengths=array(0,runs)
#lengths_leak=array(0,runs)
#endpt=array(0,runs)
#endpt_leak=array(0,runs)
#
#for(i in 1:runs){
#	dv[[i]]=array()
#	dv_leak[[i]]=array()
#	dv[[i]][1]=0
#	dv_leak[[i]][1]=0
#	y=rnorm(max,mean=delta,sd=2*psi^2)
#	j=1
#	while(abs(dv[[i]][j])<alpha&&j<=max){
#		j=j+1
#		dv[[i]][j]=dv[[i]][j-1]+y[j-1]
#		}
#	j=1
#	while(abs(dv_leak[[i]][j])<alpha&&j<=max){
#		j=j+1
#		dv_leak[[i]][j]=(1+lambda)*dv_leak[[i]][j-1]+y[j-1]
#		}
#	lengths[i]=length(dv[[i]])
#	endpt[i]=dv[[i]][length(dv[[i]])]
#	lengths_leak[i]=length(dv_leak[[i]])
#	endpt_leak[i]=dv_leak[[i]][length(dv_leak[[i]])]
#	}
# 
#layout(matrix(1:6,nrow=2,byrow=T))
#
#hist(lengths,xlab="Time until Signal",main="Without Leak")
#abline(v=mean(lengths),col="red")
#hist(endpt,xlab="Which Threshold is Reached",main="Without Leak")
#abline(v=mean(endpt),col="red")
#plot(dv[[1]],type="l",ylim=c(-alpha,alpha),xlim=c(0,max),xlab="Time",ylab="Estimate")
#for(i in 2:runs){lines(dv[[i]])}
#abline(h=-alpha)
#abline(h=alpha)
#
#hist(lengths_leak,xlab="Time Until Signal",main="With Leak")
#abline(v=mean(lengths_leak),col="red")
#hist(endpt_leak,xlab="Which Threshold is Reached",main="With Leak")
#abline(v=mean(endpt_leak),col="red")
#plot(dv_leak[[1]],type="l",ylim=c(-alpha,alpha),xlim=c(0,max),xlab="Time",ylab="Estimate")
#for(i in 2:runs){lines(dv_leak[[i]])}
#abline(h=-alpha)
#abline(h=alpha)


 