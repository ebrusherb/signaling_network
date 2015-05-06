function [prob time perf]=solve_pde_2d(Lx,Ly,d,Ux,Uy,deltax,deltay,l,b)
% dX1=(-lX1+b*(2d-1))dt+b*sqrt(d)dW1t-b*sqrt(1-d)dW2t, similar for dX2
% simplified PDE written as
% 0=(-lX+a)df/dx+(-ly-a)df/dy+D1d^2f/dx^2-2Ddd^2f/dxdy+D2d^2f/dy^2
if nargin==2 
    d=.5;
    Ux=5;
    Uy=5;
    deltax=.01;
    deltay=.05;
    l=1;
    b=1;
end

if nargin==3
    Ux=5;
    Uy=5;
    deltax=.05;
    deltay=.05;
    l=1;
    b=1;
end

D1=1/2*b^2;
Dd=-1/2*b^2;
D2=1/2*b^2;
a=b*(2*d-1);

%% setting up vectors I need

xvals=Lx:deltax:Ux;

yvals=Ly:deltay:Uy;

Ni=length(yvals);
Nj=length(xvals);

c=xvals(2:(end))';
cc=c(:,ones(Ni-1,1));
cc=transpose(cc);
xvec = cc(:)';

c=[1:(Nj-1)]';
cc=c(:,ones(Ni-1,1));
cc=transpose(cc);
jvals = cc(:)';

c=yvals(2:(end))';
cc=c(:,ones(Nj-1,1));
yvec = cc(:)';

c=[1:(Ni-1)]';
cc=c(:,ones(Nj-1,1));
ivals = cc(:)';

%% construct adjacency matrix
A=sparse((Ni-1)*(Nj-1),(Ni-1)*(Nj-1));

Aleft=D1/(deltax^2);
Adown=D2/(deltay^2);
Aupright=2*Dd/deltax/deltay;

bl_prob=0;
br_prob=1;
bd_prob=1;
bu_prob=0;
b_prob=zeros(1,(Ni-1)*(Nj-1));

bl_time=0;
br_time=0;
bd_time=0;
bu_time=0;
b_time=zeros(1,(Ni-1)*(Nj-1));

down=1+(1:(Nj-3))*(Ni-1);
up=Ni-1+(1:(Nj-3))*(Ni-1);
left=2:(Ni-2);
right=(Nj-2)*(Ni-1)+(2:(Ni-2));
SW=1;
SE=(Ni-1)*(Nj-2)+1;
NW=Ni-1;
NE=(Ni-1)*(Nj-1);
interior=1:((Ni-1)*(Nj-1));
interior([down up left right SW SE NW NE])=[];

for k=interior
    x=xvec(k);
    y=yvec(k);
    fx=-l*x+a;
    gy=-l*y-a;
    A(k,k)=-fx/deltax-gy/deltay-2*D1/(deltax^2)+2*Dd/deltax/deltay-2*D2/(deltay^2);
    A(k,k+Ni-1)=fx/deltax+D1/(deltax^2)-2*Dd/deltax/deltay;
    A(k,k+1)=gy/deltay-2*Dd/deltax/deltay+D2/(deltay^2);
    A(k,k-(Ni-1))=Aleft;
    A(k,k-1)=Adown;
    A(k,k+Ni)=Aupright;
end

for k=down
    x=xvec(k);
    y=yvec(k);
    fx=-l*x+a;
    gy=-l*y-a;
    A(k,k)=-fx/deltax-gy/deltay-2*D1/(deltax^2)+2*Dd/deltax/deltay-2*D2/(deltay^2);
    A(k,k+Ni-1)=fx/deltax+D1/(deltax^2)-2*Dd/deltax/deltay;
    A(k,k+1)=gy/deltay-2*Dd/deltax/deltay+D2/(deltay^2);
    A(k,k-(Ni-1))=Aleft;
    b_prob(k)=Adown*bd_prob;
    b_time(k)=Adown*bd_time;
    A(k,k+Ni)=Aupright;
end

for k=up
    x=xvec(k);
    y=yvec(k);
    fx=-l*x+a;
    gy=-l*y-a;
    A(k,k)=-fx/deltax-gy/deltay-2*D1/(deltax^2)+2*Dd/deltax/deltay-2*D2/(deltay^2);
    A(k,k+Ni-1)=fx/deltax+D1/(deltax^2)-2*Dd/deltax/deltay;
    b_prob(k)=(gy/deltay-2*Dd/deltax/deltay+D2/(deltay^2))*bu_prob;
    b_time(k)=(gy/deltay-2*Dd/deltax/deltay+D2/(deltay^2))*bu_time;
    A(k,k-(Ni-1))=Aleft;
    A(k,k-1)=Adown;
    b_prob(k)=b_prob(k)+Aupright*bu_prob;
    b_time(k)=b_time(k)+Aupright*bu_time;
end

for k=left
    x=xvec(k);
    y=yvec(k);
    fx=-l*x+a;
    gy=-l*y-a;
    A(k,k)=-fx/deltax-gy/deltay-2*D1/(deltax^2)+2*Dd/deltax/deltay-2*D2/(deltay^2);
    A(k,k+Ni-1)=fx/deltax+D1/(deltax^2)-2*Dd/deltax/deltay;
    A(k,k+1)=gy/deltay-2*Dd/deltax/deltay+D2/(deltay^2);
    b_prob(k)=Aleft*bl_prob;
    b_time(k)=Aleft*bl_time;
    A(k,k-1)=Adown;
    A(k,k+Ni)=Aupright;
end

for k=right
    x=xvec(k);
    y=yvec(k);
    fx=-l*x+a;
    gy=-l*y-a;
    A(k,k)=-fx/deltax-gy/deltay-2*D1/(deltax^2)+2*Dd/deltax/deltay-2*D2/(deltay^2);
    b_prob(k)=(fx/deltax+D1/(deltax^2)-2*Dd/deltax/deltay)*br_prob;
    b_time(k)=(fx/deltax+D1/(deltax^2)-2*Dd/deltax/deltay)*br_time;
    A(k,k+1)=gy/deltay-2*Dd/deltax/deltay+D2/(deltay^2);
    A(k,k-(Ni-1))=Aleft;
    A(k,k-1)=Adown;
    b_prob(k)=b_prob(k)+Aupright*br_prob;
    b_time(k)=b_time(k)+Aupright*br_time;
end

for k=SW
    x=xvec(k);
    y=yvec(k);
    fx=-l*x+a;
    gy=-l*y-a;
    A(k,k)=-fx/deltax-gy/deltay-2*D1/(deltax^2)+2*Dd/deltax/deltay-2*D2/(deltay^2);
    A(k,k+Ni-1)=fx/deltax+D1/(deltax^2)-2*Dd/deltax/deltay;
    A(k,k+1)=gy/deltay-2*Dd/deltax/deltay+D2/(deltay^2);
    b_prob(k)=Aleft*bl_prob;
    b_time(k)=Aleft*bl_time;
    b_prob(k)=b_prob(k)+Adown*bd_prob;
    b_time(k)=b_time(k)+Adown*bd_time;
    A(k,k+Ni)=Aupright;
end

for k=SE
    x=xvec(k);
    y=yvec(k);
    fx=-l*x+a;
    gy=-l*y-a;
    A(k,k)=-fx/deltax-gy/deltay-2*D1/(deltax^2)+2*Dd/deltax/deltay-2*D2/(deltay^2);
    b_prob(k)=(fx/deltax+D1/(deltax^2)-2*Dd/deltax/deltay)*br_prob;
    b_time(k)=(fx/deltax+D1/(deltax^2)-2*Dd/deltax/deltay)*br_time;
    A(k,k+1)=gy/deltay-2*Dd/deltax/deltay+D2/(deltay^2);
    A(k,k-(Ni-1))=Aleft;
    b_prob(k)=b_prob(k)+Adown*bd_prob;
    b_time(k)=b_time(k)+Adown*bd_time;
    b_prob(k)=b_prob(k)+Aupright*br_prob;
    b_time(k)=b_time(k)+Aupright*br_time;
end

for k=NW
    x=xvec(k);
    y=yvec(k);
    fx=-l*x+a;
    gy=-l*y-a;
    A(k,k)=-fx/deltax-gy/deltay-2*D1/(deltax^2)+2*Dd/deltax/deltay-2*D2/(deltay^2);
    A(k,k+Ni-1)=fx/deltax+D1/(deltax^2)-2*Dd/deltax/deltay;
    b_prob(k)=(gy/deltay-2*Dd/deltax/deltay+D2/(deltay^2))*bu_prob;
    b_time(k)=(gy/deltay-2*Dd/deltax/deltay+D2/(deltay^2))*bu_time;
    b_prob(k)=b_prob(k)+Aleft*bl_prob;
    b_time(k)=b_time(k)+Aleft*bl_time;
    A(k,k-1)=Adown;
    b_prob(k)=b_prob(k)+Aupright*bu_prob;
    b_time(k)=b_time(k)+Aupright*bu_time;
end

for k=NE
    x=xvec(k);
    y=yvec(k);
    fx=-l*x+a;
    gy=-l*y-a;
    A(k,k)=-fx/deltax-gy/deltay-2*D1/(deltax^2)+2*Dd/deltax/deltay-2*D2/(deltay^2);
    b_prob(k)=(fx/deltax+D1/(deltax^2)-2*Dd/deltax/deltay)*br_prob;
    b_time(k)=(fx/deltax+D1/(deltax^2)-2*Dd/deltax/deltay)*br_time;
    b_prob(k)=b_prob(k)+(gy/deltay-2*Dd/deltax/deltay+D2/(deltay^2))*bu_prob;
    b_time(k)=b_time(k)+(gy/deltay-2*Dd/deltax/deltay+D2/(deltay^2))*bu_time;
    A(k,k-(Ni-1))=Aleft;
    A(k,k-1)=Adown;
    b_prob(k)=b_prob(k)+Aupright*br_prob;
    b_time(k)=b_time(k)+Aupright*br_time;
end

%% solve for u
%solve A*u+b=f
fval_prob=0;
f_prob=fval_prob*ones(1,(Ni-1)*(Nj-1));

rightvec_prob=reshape(f_prob-b_prob,(Ni-1)*(Nj-1),1);
u_prob=A\rightvec_prob;

fval_time=-1;
f_time=fval_time*ones(1,(Ni-1)*(Nj-1));

rightvec_time=reshape(f_time-b_time,(Ni-1)*(Nj-1),1);
u_time=A\rightvec_time;

umat_prob=zeros(Ni-1,Nj-1);
umat_time=zeros(Ni-1,Nj-1);
xmat=zeros(Ni-1,Nj-1);
ymat=zeros(Ni-1,Nj-1);
for k=1:((Ni-1)*(Nj-1))
    i=ivals(k);
    j=jvals(k);
    umat_prob(i,j)=u_prob(k);
    umat_time(i,j)=u_time(k);
    xmat(i,j)=xvec(k);
    ymat(i,j)=yvec(k);
end
umat_prob=[umat_prob; bu_prob*ones(1,Nj-1)];
umat_prob=[bd_prob*ones(1,Nj-1); umat_prob];
umat_prob=[bl_prob*ones(Ni+1,1) umat_prob];
umat_prob=[umat_prob br_prob*ones(Ni+1,1)];

umat_time=[umat_time; bu_time*ones(1,Nj-1)];
umat_time=[bd_time*ones(1,Nj-1); umat_time];
umat_time=[bl_time*ones(Ni+1,1) umat_time];
umat_time=[umat_time br_time*ones(Ni+1,1)];

xmat=[xmat; xvals(2:Nj)];
xmat=[xvals(2:Nj); xmat];
xmat=[(xvals(1))*ones(Ni+1,1) xmat];
xmat=[xmat (xvals(end)+deltax)*ones(Ni+1,1)];

ymat=[ymat; (yvals(end)+deltay)*ones(1,Nj-1)];
ymat=[yvals(1)*ones(1,Nj-1); ymat];
ymat=[[yvals yvals(end)+deltay]' ymat];
ymat=[ymat [yvals yvals(end)+deltay]'];

check_prob=zeros(Ni+1,Nj+1);
check_time=zeros(Ni+1,Nj+1);
for i=2:(Ni)
    for j=2:(Nj)
    check_prob(i,j)=abs((-l*xmat(i,j)+a)*(umat_prob(i,j+1)-umat_prob(i,j))/deltax+(-l*ymat(i,j)-a)*(umat_prob(i+1,j)-umat_prob(i,j))/deltay+D1*(umat_prob(i,j+1)+umat_prob(i,j-1)-2*umat_prob(i,j))/deltax^2+D2*(umat_prob(i+1,j)+umat_prob(i-1,j)-2*umat_prob(i,j))/deltay^2+2*Dd*(umat_prob(i+1,j+1)-umat_prob(i+1,j)-umat_prob(i,j+1)+umat_prob(i,j))/deltax/deltay-fval_prob);
    check_time(i,j)=abs((-l*xmat(i,j)+a)*(umat_time(i,j+1)-umat_time(i,j))/deltax+(-l*ymat(i,j)-a)*(umat_time(i+1,j)-umat_time(i,j))/deltay+D1*(umat_time(i,j+1)+umat_time(i,j-1)-2*umat_time(i,j))/deltax^2+D2*(umat_time(i+1,j)+umat_time(i-1,j)-2*umat_time(i,j))/deltay^2+2*Dd*(umat_time(i+1,j+1)-umat_time(i+1,j)-umat_time(i,j+1)+umat_time(i,j))/deltax/deltay-fval_time);
    end
end

%% find prob and time
% iindex=1:(Ni+1);
% zeroi=(abs(ymat(:,1))<=0.00001);
% zeroi=iindex(zeroi);
% 
% jindex=1:(Nj+1);
% zeroj=(abs(xmat(1,:))<=0.00001);
% zeroj=jindex(zeroj);
% 
% prob=umat_prob(zeroi,zeroj);
% time=umat_time(zeroi,zeroj);

Y=reshape(ymat(:,1),1,[]);
X=reshape(xmat(1,:),1,[]);

prob=interp2(X,Y,umat_prob,0,0);
time=interp2(X,Y,umat_time,0,0);

perf=prob/time;
end


