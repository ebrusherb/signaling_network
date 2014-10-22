function x=timefunx(y)
global f b eps t e 
c = 1/e*(f-b/(2*eps)*(b-f)/e-b*t/(2*eps)+b/2);
ty = 1/e*log((y-(b-f)/e)/(-t-eps+4*eps*f/(b-2*eps*e)-(b-f)/e));
x = (-t-eps-c)*(y-(b-f)/e)/(-t-eps+4*eps*f/(b-2*eps*e)-(b-f)/e)+c+b/(2*eps)*(y-(b-f)/e).*ty;
end