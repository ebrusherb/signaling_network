N=10;

% varval=.1;
% fighting_abilities=randn(1,N)*sqrt(varval);

fighting_abilities=rand(1,N);
fighting_abilities=sort(fighting_abilities,'descend');

Lvals=[-.5:.1:0];
NL=length(Lvals);
group_learning=zeros(3,NL);

    Lx=-3;
    Ly=Lx;
%% define function from dominance to probability of signaling

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

obs_ranks=int16.empty(N+1,0);
metrics=zeros(1,its);

for r=1:its
    flips=rand(1,N*(N-1)/2);
    sigmat=zeros(N,N);
    k=1;
    for j=1:N
        for i=1:(j-1)
            if flips(k)<=probs(k)
                sigmat(i,j)=1;
            else sigmat(j,i)=1;
            end
            k=k+1;
        end
    end
    v=tiedrank(sum(sigmat,2));
    l=size(obs_ranks,2);
    replicate=0;
    while i<=l
        if isequal(v,obs_ranks(:,i))
            obs_ranks(N+1,i)=obs_ranks(N+1,i);
            replicate=1;
            i=l+1;
        else i=i+1;
        end
    end
    if replicate==0
        obs_ranks=[obs_ranks, [v;1]];
    end
    metrics(r)=corr(fighting_abilities',sum(sigmat,2),'type','Kendall');
end


