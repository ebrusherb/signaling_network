
time=0.05;
c11=0;c21=time;c31=1-time;
c12=1-time;c22=time;c32=0;


% time=0.0;
% c11=0;c21=time;c31=1-time;
% c12=1-time;c22=time;c32=0;

cvals=[c11 c21 c31; c12 c22 c32];

its=100;
[optthresh, dmatssaved, pairwiseaccuracy, pairwisetime] = example_nash_optimization(cvals,domvals,twomat);

Nc=size(cvals,1);
ER=zeros(its,Nc+5,N,N);
DT=zeros(its,Nc+5,N,N);
DP=zeros(its,Nc+5,N,N);

for i=1:its
    for j=1:(Nc+5)
        for k=1:N
            for l=setdiff(1:N,k)
                m1=min(k,l);
                m2=max(k,l);
            ER(i,j,k,l)=pairwiseaccuracy(i,j,m2,m1);
            DT(i,j,k,l)=pairwisetime(i,j,k,l);
            DP(i,j,k,l)=pairwiseaccuracy(i,j,k,l);
            end
        end
    end
end

ERmean=sum(ER,4)/N;
DTmean=sum(DT,4)/N;
DPmean=sum(DP,4)/N;

i=1;
c1=cvals(i,1);
c2=cvals(i,2);
c3=cvals(i,3);
totalperf1=c1*ERmean+c2*DTmean+c3*(1-DPmean);

i=2;
c1=cvals(i,1);
c2=cvals(i,2);
c3=cvals(i,3);
totalperf2=c1*ERmean+c2*DTmean+c3*(1-DPmean);

%%
i=i+1;
k=16;
M=max(max(max(ER(i,1,k:end,k:end))),max(max(ER(i,2,k:end,k:end))));
subplot(1,2,1)
imagesc(reshape(ER(i,1,k:end,k:end),N-k+1,[]));colorbar
caxis manual;caxis([0 M])
subplot(1,2,2)
imagesc(reshape(ER(i,2,k:end,k:end),N-k+1,[]));colorbar
caxis manual;caxis([0 M])

