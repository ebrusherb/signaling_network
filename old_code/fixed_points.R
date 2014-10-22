library(fields)
library(gplots)

rate<-function(x){
	y=(1-tanh(x))/2
	return(y)
}

rate<-function(x){
	y=0
	if(x<=(-eps)){y=1}
	if(x>(-eps)&&x<eps){y=-x/(2*eps)+.5}
	return(y)
}

fixedy<-function(x,b,t){
	r=-f/error+b/error*rate(x+t)
	return(r)
}

findy<-function(y,b,t){
	r=-f/error+b/error*rate(f/error+b/error*rate(y+t)+t)-y
	return(r)
}

fixedx<-function(y,b,t){
	r=f/error+b/error*rate(y+t)
	return(r)
}

yvals<-function(b,t){
	yvec=matrix(seq(-f/error+b/error,-f/error,length.out=1000),nrow=1)
	seek=apply(yvec,2,findy,b=b,t=t)
	seek=round(seek,5)
	zeros=array(0,3)
	if(sum(abs(diff(seek<=0)))<=1){
	# if(length(setdiff(unique(sign(seek)),0))==1){
		zeros[1]=uniroot(findy,lower=-f/error,upper=-f/error+b/error,b=b,t=t)$root
		zeros[2]=zeros[1]
		zeros[3]=zeros[1]
	} else{
	Y1=min(which(seek>0))
	Y2=min(which(seek[-(1:Y1)]<0))+Y1
	zeros[1]=uniroot(findy,upper=-f/error+b/error,lower=yvec[Y1],b=b,t=t)$root
	zeros[2]=uniroot(findy,upper=yvec[Y1],lower=yvec[Y2],b=b,t=t)$root
	zeros[3]=uniroot(findy,upper=yvec[Y2],lower=-f/error,b=b,t=t)$root
	}
	return(zeros)
}

# type<-function(x,y,b,t){
	# k=0
	# if(x[1]!=x[3]){k=7
	# } else{
	# x=x[1]
	# y=y[1]
	# if(fixedy(x,b,t)>(-f/error)&& fixedy(x,b,t)<(-f/error)+b/error&& fixedx(y,b,t)==f/error){k=1}	#1 kind of, 2 not
	# if(fixedy(x,b,t)==(-f/error)&&fixedx(y,b,t)==f/error){k=2} #1 not, 2 not
	# if(fixedy(x,b,t)>(-f/error)&& fixedy(x,b,t)<(-f/error)+b/error&&fixedx(y,b,t)>(f/error)&& fixedx(y,b,t)<(f/error)+b/error){k=3} #1 kind of, 2 kind of
	# if(fixedy(x,b,t)==(-f/error)&&fixedx(y,b,t)>(f/error)&& fixedx(y,b,t)<(f/error)+b/error){k=4} #1 not, 2 kind of
	# if(fixedy(x,b,t)>(-f/error)&& fixedy(x,b,t)<(-f/error)+b/error&&fixedx(y,b,t)==f/error+b/error){k=5} #1 kind of, 2 signaling
	# if(fixedy(x,b,t)==(-f/error)&&fixedx(y,b,t)==f/error+b/error){k=6} #1 not, 2 signaling
	# }
	# return(k)
# }

type<-function(x,y,b,t){
	k=0
	if(x[1]!=x[3]){k=7
	} else{
	x=x[1]
	y=y[1]
	if(y>(-f/error)&& y<(-f/error)+b/error&& x==f/error){k=1}	#1 kind of, 2 not
	if(y==(-f/error)&&x==f/error){k=2} #1 not, 2 not
	if(y>(-f/error)&& y<(-f/error)+b/error&&x>(f/error)&& fixedx(y,b,t)<(f/error)+b/error){k=3} #1 kind of, 2 kind of
	if(y==(-f/error)&&x>(f/error)&& x<(f/error)+b/error){k=4} #1 not, 2 kind of
	if(y>(-f/error)&& y<(-f/error)+b/error&&x==f/error+b/error){k=5} #1 kind of, 2 signaling
	if(y==(-f/error)&&x==f/error+b/error){k=6} #1 not, 2 signaling
	}
	return(k)
}

typenames<-function(v){
	l=length(v)
	names=array(,l)
	for(i in 1:l){
		if(v[i]==1){names[i]="1 kind of, 2 not"}
		if(v[i]==2){names[i]="1 not, 2 not"}
		if(v[i]==3){names[i]="1 kind of, 2 kind of"}
		if(v[i]==4){names[i]="1 not, 2 kind of"}
		if(v[i]==5){names[i]="1 kind of, 2 signaling"}
		if(v[i]==6){names[i]="1 not, 2 signaling"}
		if(v[i]==7){names[i]="2 stable nodes"}
	}
	return(names)
}

f=.5
error=.5
eps=.5
bvals=seq(0.01,5,by=0.1)
tvals=seq(0,10,by=0.5)
types<-array(0,c(length(bvals),length(tvals)))
for(i in 1:length(bvals)){
	for(j in 1:length(tvals)){
		y=yvals(bvals[i],tvals[j])
		x=apply(matrix(y,nrow=1),2,fixedx,b=bvals[i],t=tvals[j])
		types[i,j]=type(x,y,bvals[i],tvals[j])
	}
}

bifurc<-function(b){ #if T<bifurc(b) and b>2*eps*e a bifurcation occurs
	# if(b>2*eps*error){
	a=(2*eps*error+b)*f/error+2*eps^2*error-eps*b
	a=a/(2*eps*error-b)
	# } else{ a=0}
	return(a)
}
bifurcline<-apply(matrix(bvals,nrow=1),2,bifurc)

both<-function(b){ #if T<both(b) and b<2*eps*e both signal at non zero rates
	if(b<2*eps*error){
		a=(2*eps*error)*(eps-f/error)+b*(-f/error-eps)
		a=a/(2*eps*error-b)
	} else{a=0}
	return(a)
}
bothline<-apply(matrix(bvals,nrow=1),2,both)

h<-function(b){
	a=eps+f/error
	return(a)
}
hline<-apply(matrix(bvals,nrow=1),2,h)

weak<-function(b){
	a=eps-f/error-b/error
	return(a)
}
wline=apply(matrix(bvals,nrow=1),2,weak)

v=sort(unique(array(types)))
l=length(v)
mat=types
for(i in 1:l){
	mat[which(types==v[i])]=i
}
m=max(which(tvals<=min(tvals)))
M=min(which(tvals>=max(tvals)))
par(oma=c(0,0,0,4))
filled.contour(bvals, tvals, mat, levels=(0:l)+.5, col=colorpanel(l+1, "white", "blue")[-1], nlevels=l, key.axes=axis(4,at=1:l,labels=typenames(v)), plot.axes={ axis(1); axis(2); lines(bvals,bifurcline,col='red'); lines(bvals,bothline,col='red'); lines(bvals,hline,col='red'); abline(h=f/error-eps,col='red');abline(h=f/error+eps,col="red")}, main=paste("f = ",as.numeric(f),", e = ",as.numeric(error),", epsilon = ",as.numeric(eps)), xlab="Level of Feedback, b", ylab="Threshold, T")

# ; abline(v=2*eps*error,col="red")
# ; lines(bvals,eps-f/error-bvals/error,col="red")

