rate<-function(x){
	if(x<=(-T-eps)){
		a=b
	} 
	if(x>=(-T+eps)){
		a=0
	} 
	if(x>(-T-eps)&&x<(-T+eps)){
	a=b*(-1/(2*eps)*(x+T)+1/2)
	}
	return(a)
}

# p<-function(n,N,T,sd=NULL){
	# if(!is.null(sd)){sd=sd
		# } else(sd=1/(2*sqrt(3)))
	# L=sd*sqrt(3)
	# K=0:n
	# J=0:(N-1-n)
	# mat=array(,c(n+1,N-n))
	# for(k in 1:(n+1)){
		# for(j in 1:(N-n)){
			# mat[k,j]=choose(n,K[k])*choose(N-1-n,J[j])*L^(K[k]+J[j])*(-1)^(N-1-n-J[j])*1/(N-K[k]-J[j])*((L-e*(T-eps))^(N-K[k]-J[j])-(-L)^(N-K[k]-J[j]))
		# }
	# }
	# s=sum(mat)
	# s=1/(2*L)^N*choose(N-1,n)*s
	# if(n==0){
		# s=s+e*(T-eps)/(2*L)
	# }
	# return(s)
# }

p<-function(n,N,T,sd=NULL){ #given N animals, threshold T, sd of uniform distribution, return probability of receiving signals from n animals
	if(!is.null(sd)){sd=sd
		} else(sd=1/(2*sqrt(3)))
	L=sd*sqrt(3)
	c=1-e*(T-eps)/(2*L)
	s=pbeta(c,n+1,N-n)
	s=s/N
	if(n==0){
		s=s+e*(T-eps)/(2*L)
	}
	return(s)
}

step_size=0.01
int<-function(v){
	l=length(v)
	step_size/2*(sum(2*v)-v[1]-v[l])
}

N=50
L=1
eps=.5
e=.1
T=eps+.1
b=5
f=.1


probs<-function(N,T,sd=NULL){ #return vector of probabilities of receiving signals from 0 to N-1 animals
	if(!is.null(sd)){sd=sd
		} else(sd=1/(2*sqrt(3)))
	L=sd*sqrt(3)
	poss=matrix(0:(N-1),nrow=1)
	v=apply(poss,2,p,N=N,T=T,sd=sd)
	return(v)
	}
	
	
well<-function(x,y){
	w=f*(y-x)+e/2*(x^2+y^2)-b*(x*rate(y)+y*rate(x))
	return(w)
}

eq1<-c(f/e+b/e*rate(-f/e),-f/e)
eq2<-c(f/e,-f/e+b/e*rate(f/e))
eq3<-c(b*(eps-T)/(2*eps*e+b)-2*eps*f/(b-2*eps*e),b*(eps-T)/(2*eps*e+b)+2*eps*f/(b-2*eps*e))


N=50

Tvals=seq(eps,eps+1,by=0.5)
layout(matrix(1:6,nrow=2,byrow=TRUE))
e=0.1
for(i in 1:length(Tvals)){

plot(0:(N-1),probs(N,Tvals[i],1/sqrt(3)),type="l",xlab="Number of Signalers",ylab="Probability",main=paste("e =",e,", eps =",eps,", T =",Tvals[i]),ylim=c(0,.1))
}
e=0.9
for(i in 1:length(Tvals)){

plot(0:(N-1),probs(N,Tvals[i],1/sqrt(3)),type="l",xlab="Number of Signalers",ylab="Probability",main=paste("e =",e,", eps =",eps,", T =",Tvals[i]),ylim=c(0,.5))
}

