function a=rate(x)
a=zeros(length(x),1);
v=1:length(x);
global eps t b
a(x<=-t-eps)=b;
middle=intersect(v(x>-t-eps),v(x<-t+eps));
a(middle)=b*(-1/(2*eps)*(x(middle)+t)+1/2);
a(x>=-t+eps)=0;
end