f=.5
e=.5
T=1
eps=.1

c=1

sd=c^2/(2*e)
mu=f/e
h=2

Sx=(f-e*(-T+eps))*dnorm(-T+eps,mean=mu,sd=sd)*(pnorm(h,mean=-mu,sd=sd)-pnorm(-T+eps,mean=-mu,sd=sd))+(-f-e*h)*dnorm(h,mean=-mu,sd=sd)*(pnorm(h,mean=mu,sd=sd)-pnorm(h,mean=mu,sd=sd))

Sy=(f-e*h)*dnorm(h,mean=mu,sd=sd)*(pnorm(h,mean=-mu,sd=sd)-pnorm(-T+eps,mean=-mu,sd=sd))+(-f-e*(-T-eps))*dnorm(-T-eps,mean=-mu,sd=sd)*(dnorm(h,mean=mu,sd=sd)-dnorm(-T+eps,mean=mu,sd=sd))