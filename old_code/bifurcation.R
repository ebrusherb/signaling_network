e=0.5
eps=0.5

bifurc<-function(b,f){ #if T<bifurc(b) and b>2*eps*e a bifurcation occurs
	if(b>2*eps*error){
	a=(2*eps*error+b)*f/error+2*eps^2*error-eps*b
	a=a/(2*eps*error-b)
	} else{ a=-1}
	return(max(a,-1))
}


both<-function(b,f){ #if T<both(b) and b<2*eps*e both signal at non zero rates
	if(b<2*eps*error){
		a=(2*eps*error)*(eps-f/error)+b*(-f/error-eps)
		a=a/(2*eps*error-b)
	} else{a=-1}
	return(max(a,-1))
}

lowerline<-function(b,f){
	a=f/error-eps
	return(a)
}


upperline<-function(b,f){
	a=eps+f/error
	return(a)
}

error=0.5
eps=.1
lb=200
lf=200
bvals=seq(0.1,10,length.out=lb)
fvals=seq(0,2,length.out=lf)

bifurcmat<-array(,c(lb,lf))
bothmat<-array(,c(lb,lf))
lowermat<-array(,c(lb,lf))
uppermat<-array(,c(lb,lf))

for(i in 1:lb){
	for(j in 1:lf){
		bifurcmat[i,j]=bifurc(bvals[i],fvals[j])
		bothmat[i,j]=both(bvals[i],fvals[j])
		lowermat[i,j]=lowerline(bvals[i],fvals[j])
		uppermat[i,j]=upperline(bvals[i],fvals[j])
	}
}



error=0.5
eps=.1
lb=500
fmax=.3
bvals=seq(0.01,1,length.out=lb)

quartz()
pdf(file='bifurc_diagram.pdf')

layout(matrix(1:2,nrow=1))
f=0.01
plot(bvals,apply(matrix(bvals,nrow=1),2,bifurc,f=f),type='l',ylim=c(0,eps+eps*error+.25),xlab="Feedback level, b",ylab="Threshold, T",cex.lab=1.5,mgp=c(2.2,1,0))
# ,main=paste("f = ",round(f,2),", e = ",e,", eps = ",eps)
title(main=paste("f =",f),cex.main=2,line=1)
lines(bvals,apply(matrix(bvals,nrow=1),2,both,f=f))
lines(bvals,apply(matrix(bvals,nrow=1),2,lowerline,f=f))
lines(bvals,apply(matrix(bvals,nrow=1),2,upperline,f=f))
legend(x=.1,y=eps+f/error+.1,bty="n",legend="neither animal signals")
legend(x=.5,y=.03,bty="n",legend="2 stable equilibria")
legend(x=-.095,y=0.02,bty="n",legend="both signal")
legend(x=0.07,y=0.12,bty="n",legend="weaker animal signals, stronger does not")

for(f in c(eps*error+.15)){
plot(bvals,apply(matrix(bvals,nrow=1),2,lowerline,f=f),type='l',ylim=c(0,eps+eps*error+.25),xlab="Feedback level, b",ylab="Threshold, T",cex.lab=1.5,mgp=c(2.2,1,0))
title(main=paste("f =",f),cex.main=2,line=1)
lines(bvals,apply(matrix(bvals,nrow=1),2,upperline,f=f))
legend(x=.08,y=eps+f/error+.1,bty="n",legend="neither animal signals")
legend(x=.08,y=eps+f/error-.1,bty="n",legend="weaker animal signals, stronger does not")
legend(x=.08,y=-eps+f/error-.1,bty="n",legend="weaker animal signals maximally, stronger does not")
}

graphics.off()

library(fields)



bvals<-seq(0,5,length.out=500)
tvals<-seq(0,5,length.out=500)
fvals<-array(,c(length(bvals),length(tvals)))
for(i in 1:length(bvals)){
	for(j in 1:length(tvals)){
		if(tvals[j]<eps){
			fvals[i,j]=abs((2*eps*e-bvals[i])*e*(tvals[j]-eps)/(2*eps*e+bvals[i]))
		} else{
			fvals[i,j]=e*(tvals[j]-eps)}
	}
}

quartz()
setwd("/Users/eleanorbrush/Documents/research/signaling_network")
postscript("heatmap_bifurcation_v1.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 6.83,bg="white",height=6.83)
image.plot(bvals,tvals,fvals,xlab="b",ylab="T")
abline(h=eps,col="white")
lines(rep(2*eps*e,500),seq(0,eps,length.out=500),col='white')
legend(eps*e-.75,eps/2+.25,'both',bty='n',text.col='white',cex=1.5)
legend(2*eps*e+.75,eps/2+.25,'bistable',bty='n',text.col='white',cex=1.5)
legend(eps*e,eps+1,'neither',bty='n',text.col='white',cex=1.5)
title(paste("e = ",e,", eps = ",eps))
graphics.off()