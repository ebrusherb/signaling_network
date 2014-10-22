function [u]=pdesolver(a,b,c,z,deltax,bl,br)
%PDE written as 0=(a(y+b)df/dy+c*f/dy^2


%% setting up vectors I need

xvals=-z:deltax:z;

Nj=length(xvals);

jvals = reshape(1:(Nj-2),1,Nj-2);

%% construct adjacency matrix
A=sparse(Nj-2,Nj-2);
bvec=zeros(1,Nj-2);

left=1;
right=Nj-2;
interior=2:(Nj-3);

for k=interior
    x=xvals(k+1);
    fx=a*x+b;
    A(k,k)=-fx/deltax-2*c/(deltax^2);
    A(k,k+1)=fx/deltax+c/(deltax^2);
    A(k,k-1)=c/(deltax^2);
end

for k=left
    x=xvals(k+1);
    fx=a*x+b;
    A(k,k)=-fx/deltax-2*c/(deltax^2);
    A(k,k+1)=fx/deltax+c/(deltax^2);
    bvec(k)=c/(deltax^2)*bl;
end

for k=right
    x=xvals(k+1);
    fx=a*x+b;
    A(k,k)=-fx/deltax-2*c/(deltax^2);
    A(k,k-1)=c/(deltax^2);
    bvec(k)=(fx/deltax+c/(deltax^2))*br;
end

%% solve for u
%solve A*u+b=f

fval=0;
f=fval*ones(1,(Nj-2));

rightvec=reshape(f-bvec,(Nj-2),1);
u=A\rightvec;
u=[bl; u; br];

check=zeros(1,Nj-2);
for j=2:(Nj-1)
    check(j-1)=abs((a*xvals(j)+b)*(u(j+1)-u(j))/deltax+c*(u(j+1)+u(j-1)-2*u(j))/deltax^2-fval);
end


end


