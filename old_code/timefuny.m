function y=timefuny(x)
global f b eps t e 
c=1/e*(-f-b/(2*eps)*f/e-b*t/(2*eps)+b/2);
tx=1/e*log((x-f/e)/(-t+eps-4*eps*f/(b-2*eps*e)-f/e));
y = (-t+eps-c)*(x-f/e)/(-t+eps-4*eps*f/(b-2*eps*e)-f/e)+c+b/(2*eps)*(x-f/e).*tx;
end