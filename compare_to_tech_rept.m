function [theirs mine]=compare_to_tech_rept(z,d,l,b,deltax)

ztilde=z/(2*b*(2*d-1));
atilde=(2*d-1)^2;
lambda=-l;

theirs=(erf(sqrt(atilde/lambda)*(1+lambda*ztilde))-erf(sqrt(atilde/lambda)))/(erf(sqrt(atilde/lambda)*(1+lambda*ztilde))-erf(sqrt(atilde/lambda)*(1-lambda*ztilde)));    

[x,~,~]=solve_pde_1d(-z,d,z,deltax,l,b);
mine=x;

end