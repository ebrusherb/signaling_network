%% define function from dominance to probability of signaling
Ly=-2;
Lx=Ly;

deps=.05;
dvals=0:deps:1;
Nd=length(dvals);
fun_d=zeros(3,Nd);

for j=1:Nd
    d=dvals(j);
    [x,y,z]=solve_pde(Lx,Ly,d,5,5,.05,.05,1,1);
    fun_d(:,j)=[x; y; z];
end

%% generate 'true' dominance matrix 
N=50;

% varval=.1;
% fighting_abilities=randn(1,N)*sqrt(varval);

fighting_abilities=rand(1,N);

fighting_abilities=sort(fighting_abilities,'descend');

dvec=zeros(1,N*(N-1)/2);
probs=dvec;
times=dvec;

k=0;
for i=1:N
    for j=(i+1):N
        k=k+1;
        diff=fighting_abilities(i)-fighting_abilities(j);
        d=exp(diff)/(exp(diff)+1);
        d=floor(round(d/deps))*deps;
        dvec(k)=d;
        probs(k)=fun_d(1,abs(d-dvals)<0.001);
        times(k)=fun_d(2,abs(d-dvals)<0.001);
    end
end

%% simulate signaling matrices

its=1000;
metrics=zeros(2,its);

for r=1:its
    flips=rand(1,N*(N-1)/2);
    right=zeros(1,N*(N-1)/2);
    right(flips<=probs)=1;
    metrics(1,r)=sum(right)/length(right);
end

mean_frac_right=mean(metrics(1,:));
mean_time=mean(times);


