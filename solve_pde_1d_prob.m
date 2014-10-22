function prob =solve_pde_1d_prob(t,d,l,b,deltax)
%dY=(-lY+2*b*(2d-1))dt+2*b*sqrt(d)dW1t-2*b*sqrt(1-d)dW2t
%simplified PDE written as 0=(-ly+2*a)df/dy+4D1d^2f/dy^2
Lx=-t;
Ux=t;

if nargin==1 
    d=.5;
    l=1;
    b=1;
end

if nargin==2
    l=1;
    b=1;
end
%%
D1=1/2*b^2;
a=2*b*(2*d-1);

%% setting up vectors I need

xvals=Lx:deltax:Ux;

Nj=length(xvals);

jvals = reshape(1:(Nj-2),1,Nj-2);

%% construct adjacency matrix
A=sparse(Nj-2,Nj-2);

bl_prob=0;
br_prob=1;
b_prob=zeros(1,Nj-2);

bl_time=0;
br_time=0;
b_time=zeros(1,(Nj-2));

left=1;
right=Nj-2;
interior=2:(Nj-3);

for k=interior
    x=xvals(k+1);
    fx=-l*x+a;
    A(k,k)=-fx/deltax-8*D1/(deltax^2);
    A(k,k+1)=fx/deltax+4*D1/(deltax^2);
    A(k,k-1)=4*D1/(deltax^2);
end

for k=left
    x=xvals(k+1);
    fx=-l*x+a;
    A(k,k)=-fx/deltax-8*D1/(deltax^2);
    A(k,k+1)=fx/deltax+4*D1/(deltax^2);
    b_prob(k)=4*D1/(deltax^2)*bl_prob;
    b_time(k)=4*D1/(deltax^2)*bl_time;
end

for k=right
    x=xvals(k+1);
    fx=-l*x+a;
    A(k,k)=-fx/deltax-8*D1/(deltax^2);
    A(k,k-1)=4*D1/(deltax^2);
    b_prob(k)=(fx/deltax+4*D1/(deltax^2))*br_prob;
    b_time(k)=(fx/deltax+4*D1/(deltax^2))*br_time;
end

%% solve for u
%solve A*u+b=f

fval_prob=0;
f_prob=fval_prob*ones(1,(Nj-2));

rightvec_prob=reshape(f_prob-b_prob,(Nj-2),1);
u_prob=A\rightvec_prob;
u_prob=[bl_prob; u_prob; br_prob];

fval_time=-1;
f_time=fval_time*ones(1,(Nj-2));

rightvec_time=reshape(f_time-b_time,(Nj-2),1);
u_time=A\rightvec_time;
u_time=[bl_time; u_time; br_time];

check_prob=zeros(1,Nj-2);
check_time=zeros(1,Nj-2);
for j=2:(Nj-1)
    check_prob(j-1)=abs((-l*xvals(j)+a)*(u_prob(j+1)-u_prob(j))/deltax+4*D1*(u_prob(j+1)+u_prob(j-1)-2*u_prob(j))/deltax^2-fval_prob);
    check_time(j-1)=abs((-l*xvals(j)+a)*(u_time(j+1)-u_time(j))/deltax+4*D1*(u_time(j+1)+u_time(j-1)-2*u_time(j))/deltax^2-fval_time);
end


%% find prob and time
prob=interp1(xvals,u_prob,0);

end


