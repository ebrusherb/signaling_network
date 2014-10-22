
%%
opt_its=10;
N=10;

fighting_abilities=20*rand(1,N);
% fighting_abilities=normrnd(0,10,1,N);
% fighting_abilities=lognrnd(0,1,1,N);
fighting_abilities=sort(fighting_abilities,'descend');

dmat=zeros(N); %turn difference in fighting abilities into probabilitise of winning
extendeddomvals=[1-domvals(end:-1:2) domvals];
v=1:(Nd+Nd-1);
deps=.1;

for i=1:N
    for j=[1:(i-1),(i+1):N]
        diff=fighting_abilities(i)-fighting_abilities(j);
        d=exp(diff)/(exp(diff)+1);
        d=floor(round(d/deps))*deps;
        dmat(i,j)=v(abs(d-extendeddomvals)<=deps/2);
    end
end

c1=0;
c2=.1;
c3=1;

%twomat=zeros(Nl,Nd,Nt,Nt,2);
perf=zeros(2,Nl,Nd,Nt,Nt); %individual, leak, dominance, thresholds

perf(1,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2))+c3*(1-twomat(:,2:Nd,:,:,1));
perf(2,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2))+c3*(twomat(:,2:Nd,:,:,1));

perf(1,:,1,:,:)=c2*(twomat(:,1,:,:,2))+c3*(1-twomat(:,1,:,:,1));
perf(2,:,1,:,:)=c2*(twomat(:,1,:,:,2))+c3*(twomat(:,1,:,:,1));


%%
l=2;

Tvals=zeros(N,opt_its);
Tvals(:,1)=2*ones(N,1);
% Tvals(:,1)=randi([1 Nt],N,1);
perfvals=zeros(N,opt_its);

for i=1:N
    opp_thresh=Tvals([1:(i-1),(i+1):N],1);
    ds=dmat(i,[1:(i-1),(i+1):N]);
    perfsum=0;
            for q=1:(N-1)
                if ds(q)<=Nd-1
                    perfsum=perfsum+perf(2,l,Nd-ds(q)+1,Tvals(q,1),Tvals(i,1)); 
                else 
                    perfsum=perfsum+perf(1,l,ds(q)-Nd+1,Tvals(i,1),Tvals(q,1));
                end
            end  
    perfvals(i,1)=perfsum;
end

count=1;
while count<=opt_its
    for i=1:N
        opp_thresh=Tvals([1:(i-1),(i+1):N],count);
        ds=dmat(i,[1:(i-1),(i+1):N]);
        perftest=zeros(1,Nt);
        for j=1:Nt
            perfsum=0;
            for q=1:(N-1)
                if ds(q)<=Nd-1
                    perfsum=perfsum+perf(2,l,Nd-ds(q)+1,Tvals(q,count),j);
                else 
                    perfsum=perfsum+perf(1,l,ds(q)-Nd+1,j,Tvals(q,count));
                end
            end 
            perftest(j)=perfsum;
        end
        [p,n]=min(perftest);
        Tvals(i,count+1)=n;
        perfvals(i,count+1)=p;
    end
    if sum(Tvals(:,count+1)==Tvals(:,count))==N
        maxit=count+2;
        count=opt_its+1;
    end
    count=count+1;
end

Tvals(:,maxit:end)=[];
perfvals(:,maxit:end)=[];

probmat=zeros(N,N);
timemat=zeros(N,N);

%twomat=zeros(Nl,Nd,Nt,Nt,2);

for i=1:N
    for j=[1:(i-1),(i+1):N]
        if dmat(i,j)<=Nd-1
            probmat(i,j)=1-twomat(l,Nd-dmat(i,j)+1,Tvals(j,end),Tvals(i,end),1);
            timemat(i,j)=twomat(l,Nd-dmat(i,j)+1,Tvals(j,end),Tvals(i,end),2);
        else
            probmat(i,j)=twomat(l,dmat(i,j)-Nd+1,Tvals(i,end),Tvals(j,end),1);
            timemat(i,j)=twomat(l,dmat(i,j)-Nd+1,Tvals(i,end),Tvals(j,end),2);
        end
    end
end

%%
r=ceil(max(fighting_abilities))-floor(min(fighting_abilities));
bins=floor(min(fighting_abilities)):floor(r/5):ceil(max(fighting_abilities));
figure
subplot(2,2,1)
imagesc(probmat)
colorbar

subplot(2,2,2)
imagesc(timemat)
colorbar

subplot(2,2,3)
% hist(Tvals(:,end))
hist(sum(probmat,2))

subplot(2,2,4)
plot(fighting_abilities,sum(probmat,2))