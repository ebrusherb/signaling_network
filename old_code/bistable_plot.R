f=0.05
eps=2
e=.5
t=1
b=5

rate<-function(x){
	if(x<=(-t-eps)){
		a=b
	} 
	if(x>=(-t+eps)){
		a=0
	} 
	if(x>(-t-eps)&&x<(-t+eps)){
	a=b*(-1/(2*eps)*(x+t)+1/2)
	}
	return(a)
}


timefun<-function(y){
	c=1/e*(f-b/(2*eps)*(b-f)/e-b*t/(2*eps)+b/2)
	ty=1/e*log((y-(b-f)/e)/(-t-eps+4*eps*f/(b-2*eps*e)-(b-f)/e))
	x=(-t-eps-c)*(y-(b-f)/e)/(-t-eps+4*eps*f/(b-2*eps*e)-(b-f)/e)+c+b/(2*eps)*(y-(b-f)/e)*ty
	return(x)
}

timefun2<-function(x){
	c=1/e*(-f-b/(2*eps)*f/e-b*t/(2*eps)+b/2)
	tx=1/e*log((x-f/e)/(-t+eps-4*eps*f/(b-2*eps*e)-f/e))
	y=(-t+eps-c)*(x-f/e)/(-t+eps-4*eps*f/(b-2*eps*e)-f/e)+c+b/(2*eps)*(x-f/e)*tx
	return(y)
}

xrange=c(-t-eps-1,f/e+rate(-f/e)/e+1)
yrange=c(-t-eps-1,-f/e+rate(f/e)/e+1)

plot(xrange,yrange,col='white',xlab="",ylab="")

xvals=matrix(seq(min(xrange),max(xrange),length.out=500),nrow=1)
yvals=matrix(seq(min(yrange),max(yrange),length.out=500),nrow=1)

xisocline=f/e+apply(yvals,2,rate)/e
yisocline=-f/e+apply(xvals,2,rate)/e
# lines(xisocline,yvals)
# lines(xvals,yisocline)

xvals=matrix(seq(-t-eps,-t+eps-4*eps*f/(b-2*eps*e),length.out=500),nrow=1)
midline=xvals+4*eps*f/(b-2*eps*e)
lines(xvals,midline)

yvals=matrix(seq(-t-eps,-t-eps+4*eps*f/(b-2*eps*e),length.out=500),nrow=1)
line=apply(yvals,2,timefun)
lines(line,yvals)

xvals=matrix(seq(-t+eps-4*eps*f/(b-2*eps*e),-t+eps,length.out=500),nrow=1)
line=apply(xvals,2,timefun2)
lines(xvals,line)

xvals=matrix(seq(min(xrange)-1,timefun(-t-eps),length.out=500),nrow=1)
xhat=timefun(-t-eps)
yhat=-t-eps
line=(yhat-(b-f)/e)/(xhat-(b+f)/e)*(xvals-(b+f)/e)+(b-f)/e
lines(xvals,line)

xvals=matrix(seq(-t+eps,max(xrange)+1,length.out=500),nrow=1)
xhat=-t+eps
yhat=timefun2(-t+eps)
line=(yhat+f/e)/(xhat-f/e)*(xvals-f/e)-f/e
lines(xvals,line)

points(f/e,-f/e+rate(f/e)/e)
points(f/e+rate(-f/e)/e,-f/e)
points(b*(eps-t)/(2*eps*e+b)-2*eps*f/(b-2*eps*e),b*(eps-t)/(2*eps*e+b)+2*eps*f/(b-2*eps*e),col='red')
lines(seq(-t-eps,-t+eps,length.out=500),rep(-t+eps,500),lty=2)
lines(seq(min(xrange)-1,-t+eps,length.out=500),rep(-t-eps,500),lty=2)
lines(rep(-t+eps,500),seq(-t-eps,max(yrange)+1,length.out=500),lty=2)
lines(rep(-t-eps,500),seq(-t-eps,-t+eps,length.out=500),lty=2)
# abline(h=-t-eps+4*eps*f/(b-2*eps*e),lty=2)
lines(seq(min(xrange)-1,-t-eps,length.out=500),rep(-t-eps+4*eps*f/(b-2*eps*e),500),lty=2)
lines(rep(-t+eps-4*eps*f/(b-2*eps*e),500),seq(-t+eps,max(yrange)+1,length.out=500),lty=2)
abline(h=0,v=0)