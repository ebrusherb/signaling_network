a=.1
b=10
c=1



delta=0
x=seq(a-b+delta,a+b-delta,length.out=100)
y=seq(-a-b+delta,-a+b-delta,length.out=100)

fixedy<-function(y,a,b,c){
	r=-a+b-b*tanh(c*(a+b-b*tanh(c*(y+1))+1))-y
	return(r)
}

fixedx<-function(y,a,b,c){
	r=a+b-b*tanh(c*(y-1))
	return(r)
}

# yvals<-function(a,b,c){
	# d=2*b/10
	# zeros=array(0,3)
	# zeros[1]=uniroot(fixedy,lower=-a-b-1,upper=-a-b+d,a=a,b=b,c=c)$root
	# if(sign(fixedy(-a-b+d,a,b,c))==sign(fixedy(-a+b-d,a,b,c))){
		# zeros[2]=zeros[1]
		# zeros[3]=zeros[1]
	# } else{
	# zeros[2]=uniroot(fixedy,lower=-a-b+d,upper=-a+b-d,a=a,b=b,c=c)$root
	# zeros[3]=uniroot(fixedy,lower=-a+b-d,upper=-a+b+1,a=a,b=b,c=c)$root
	# }
	# return(zeros)
# }

yvals<-function(a,b,c){
	yvec=matrix(seq(-a-1,2*b,length.out=100),nrow=1)
	seek=apply(yvec,2,fixedy,a=a,b=b,c=c)
	d=.1
	zeros=array(0,3)
	if(sum(abs(diff(seek>0)))==1){
		zeros[1]=uniroot(fixedy,lower=-a-1,upper=2*b,a=a,b=b,c=c)$root
		zeros[2]=zeros[1]
		zeros[3]=zeros[1]
	} else{
	Y1=min(which(seek<0))
	Y2=min(which(seek[-(1:Y1)]>0))+Y1
	zeros[1]=uniroot(fixedy,lower=-a-1,upper=yvec[Y1],a=a,b=b,c=c)$root
	zeros[2]=uniroot(fixedy,lower=yvec[Y1],upper=yvec[Y2],a=a,b=b,c=c)$root
	zeros[3]=uniroot(fixedy,lower=yvec[Y2],upper=2*b,a=a,b=b,c=c)$root
	}
	return(zeros)
}

# xvals=apply(matrix(yvals(a,b,c),nrow=1),2,fixedx,a=a,b=b,c=c)

# bvals=seq(12.6,12.8,length.out=1000)
# l=length(bvals)
# bifurcation<-array(0,c(l,2,3))
# for(i in 1:l){
	# bifurcation[i,2,]=yvals(a,bvals[i],c)
	# bifurcation[i,1,]=apply(matrix(bifurcation[i,2,],nrow=1),2,fixedx,a=a,b=bvals[i],c=c)
# }

# # layout(matrix(1:2,nrow=2))
# plot(bvals,bifurcation[,2,1],type="l",ylim=range(bifurcation[,2,]))
# lines(bvals,bifurcation[,2,2])
# lines(bvals,bifurcation[,2,3])
# # plot(bvals,bifurcation[,1,1],type="l",ylim=range(bifurcation[,1,]))
# # lines(bvals,bifurcation[,1,2])
# # lines(bvals,bifurcation[,1,3])

bvals=seq(1,20,length.out=100)
lb=length(bvals)
cvals=seq(.5,20,length.out=100)
lc=length(cvals)
bifurcation<-array(0,c(lb,lc,2,3))
for(i in 1:lc){
	for(j in 1:lb){
	bifurcation[j,i,2,]=yvals(a,bvals[j],cvals[i])
	bifurcation[j,i,1,]=apply(matrix(bifurcation[j,i,2,],nrow=1),2,fixedx,a=a,b=bvals[j],c=cvals[i])
	print(c(i,j))
	}
}

# layout(matrix(1:2,nrow=2))
k=100
plot(bifurcation[k,,2,1],type="l",ylim=range(bifurcation[k,,2,]))
# lines(bifurcation[k,,2,2])
lines(bifurcation[k,,2,3])

X=rep(bvals,lc)
Y=array(,lb*lc)
for(i in 1:lc){
	Y[(1:lb)+lb*(i-1)]=cvals[i]
}
Z=array(bifurcation[,,2,3])
scatterplot3d(X,Y,Z,color="red")

s3d<-scatterplot3d(bvals,cvals,seq(min(bifurcation[,,2,]),max(bifurcation[,,2,]),length.out=lb),type="n",grid=FALSE, angle = 70,zlab = "Fixed Point",xlab = "Feedback", ylab = "Threshold",main = "")

for(i in lb:1)
s3d$points3d(rep(bvals[i], lc), cvals, bifurcation[i,,2,1], type = "l",pch=8)
for(i in lb:1)
s3d$points3d(rep(bvals[i], lc), cvals, bifurcation[i,,2,3], type = "l",pch=8)


# for(i in lc:1)
# s3d$points3d(bvals, rep(cvals[i], lb), bifurcation[,i,2,1], type = "l")
# for(i in lc:1)
# s3d$points3d(bvals, rep(cvals[i], lb), bifurcation[,i,2,3], type = "l")