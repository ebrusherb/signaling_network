eps=.5
e=.5
T=1
b=5


makesignet<-function(N,b,T,sd=NULL,dist=NULL){
	if(!is.null(sd)){sd=sd
		} else(sd=1/(2*sqrt(3)))
# abilities=rnorm(N,mean=0,sd=sd)
if(is.null(dist)){
abilities=runif(N,-sd*sqrt(3),sd*sqrt(3))
} else{ abilities=rnorm(N,mean=0,sd=sd)}
abilities=sort(abilities,decreasing=TRUE)
fmat<-array(,c(N,N))
for(i in 1:N){
	fmat[i,]=abilities[i]-abilities
}

if(T>=eps){
	thresh=e*(T-eps)
	} else{thresh=e*(eps-T)*abs((2*eps*e-b)/(2*eps*e+b))}

sigmat<-array(0,c(N,N))
for(i in 1:(N-1)){
	for(j in (i+1):N){
		f=fmat[i,j]
		if(f>e*(T+eps)){
			sigmat[i,j]=0
			sigmat[j,i]=b
		}
		if(f<e*(T+eps)&&f>thresh){
			sigmat[i,j]=0
			sigmat[j,i]=b*(f/(2*eps*e)-T/(2*eps)+1/2)
		}
		if(f<thresh&&T>=eps){
			sigmat[i,j]=sigmat[j,i]=0
		}
		if(f<thresh&&T<eps){
		if(b<2*eps*e){
			sigmat[i,j]=b*(b*(T-eps)/(2*eps*(2*eps*e+b)+2*eps*f/(2*eps*(b-2*eps*e))-T/2*eps+1/2))
			sigmat[j,i]=b*(b*(T-eps)/(2*eps*(2*eps*e+b)-2*eps*f/(2*eps*(b-2*eps*e))-T/2*eps+1/2))
		} else {
			x=runif(1,0,1)
			prop=(.5*b*(eps-T)/(2*eps*e+b)+3*eps*f/(b-2*eps*e))/(b*(eps-T)/(2*eps*e+b)+2*eps*f/(b-2*eps*e))
			if(x<prop){ 
				sigmat[i,j]=0
				sigmat[j,i]=b*(f/(2*eps*e)-T/(2*eps)+1/2)
			} else {
				sigmat[i,j]=b*(-f/(2*eps*e)-T/(2*eps)+1/2)
				sigmat[j,i]=0
				}
		}	
		}
	}
}
return(sigmat)
}

binsig<-function(name){ #puts 1 where there nonzero entries of name and 0s elsewhere
	n<-dim(name)[1]
	x<-array(0,c(n,n))
for(i in 1:n){for(j in 1:n){if(name[i,j]!=0){x[i,j]=1}}}
return(x)
}

myhist<-function(v){
	u=unique(sort(v,decreasing=FALSE))
	p=array(,length(u))
	for(i in 1:length(u)){
		p[i]=length(which(v==u[i]))/length(v)
	}
	return(p)
}


reps=500
numsgs<-array(,N*reps)
for(i in 1:reps){
	# M=makesignet(N,b,T,dist='norm')
	M=makesignet(N,b,T)
	v=apply(binsig(M),2,sum)
	numsgs[(1:N)+(i-1)*N]=v
}

plot(sort(unique(numsgs)),myhist(numsgs))



N=50

Tvals=seq(eps,eps+1,by=0.5)
layout(matrix(1:6,nrow=2,byrow=TRUE))
e=0.1
for(i in 1:length(Tvals)){
numsgs<-array(,N*reps)
for(j in 1:reps){
	M=makesignet(N,b,T,dist='norm')
	# M=makesignet(N,b,T)
	v=apply(binsig(M),2,sum)
	numsgs[(1:N)+(j-1)*N]=v
}

plot(sort(unique(numsgs)),myhist(numsgs),xlab="Number of Signalers",ylab="Probability",main=paste("e =",e,", eps =",eps,", T =",Tvals[i]),ylim=c(0,.1))
lines(0:(N-1),probs(N,T))
}
e=0.9
for(i in 1:length(Tvals)){
numsgs<-array(,N*reps)
for(j in 1:reps){
	M=makesignet(N,b,T,dist='norm')
	# M=makesignet(N,b,T)
	v=apply(binsig(M),2,sum)
	numsgs[(1:N)+(j-1)*N]=v
}

plot(sort(unique(numsgs)),myhist(numsgs),xlab="Number of Signalers",ylab="Probability",main=paste("e =",e,", eps =",eps,", T =",Tvals[i]),ylim=c(0,.1))
lines(0:(N-1),probs(N,T))
}