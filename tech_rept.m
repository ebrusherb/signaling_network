function [theirs]=tech_rept(z,d,l,b)

ztilde=z/(2*b*(2*d-1));
atilde=(2*d-1)^2;
lambda=-l;

if atilde/lambda>0
    erfplus=erf(sqrt(atilde/lambda)*(1+lambda*ztilde));
    erfminus=erf(sqrt(atilde/lambda)*(1-lambda*ztilde));
    erfdot=erf(sqrt(atilde/lambda));
else 
    erfplus=erfi(-sqrt(atilde/lambda)*(1+lambda*ztilde)/1j)/1j;
    erfminus=erfi(-sqrt(atilde/lambda)*(1-lambda*ztilde)/1j)/1j;
    erfdot=erfi(-sqrt(atilde/lambda)/1j)/1j;
end

theirs=(erfplus-erfdot)/(erfplus-erfminus);    

end