Lxvals=-1:.05:-.1;
Lyvals=-1:.05:-.1;
deps=.1;
dvals=0:deps:1;
NLx=length(Lxvals);
NLy=length(Lyvals);
Nd=length(dvals);

l=1;
b=1;

%%
f=zeros(NLx,NLy,Nd,3);

for i=1:NLx
    for j=1:NLy
        for k=1:Nd
            Lx=Lxvals(i);
            Ly=Lyvals(j);
            d=dvals(k);
            [x,y,~]=solve_pde(Lx,Ly,d,5,5,.05,.05,l,b);
            f(i,j,k,:)=[x,1-x,y];
        end
    end
end


%%
Lx_best=zeros(NLy,Nd,2);

for j=1:NLy
    for k=1:Nd
        d=dvals(k);
        perfvals1=-f(:,j,k,3)+f(:,j,k,1).*(1-f(:,j,k,3)); %receiving a signal 
        if d<.5
            perfvals2=-f(:,j,k,3)+f(:,j,k,2).*(1-f(:,j,k,3)); %getting it 'right'
        else perfvals2=perfvals1;
        end
        [~,w1]=max(perfvals1);
        [~,w2]=max(perfvals2);
        Lx_best(j,k,:)=Lxvals([w1 w2]);
    end
end

%%
NEx={};
NEy={};

for i=0:(Nd-1)
    pp1=interp1(Lyvals,Lx_best(:,1+i,1),'linear','pp');
    pp2=interp1(Lyvals,Lx_best(:,Nd-i,1),'linear','pp');
    Tystar=Lyvals(abs(ppval(pp2,ppval(pp1,Lyvals))-Lyvals)<0.001);
    Txstar=ppval(pp1,Tystar);
    NEx{i+1}=Txstar;
    NEy{i+1}=Tystar;
end


%%
its=100;
N=50;

fighting_abilities=20*rand(1,N);
fighting_abilities=sort(fighting_abilities,'descend');

dmat=zeros(N);
v=1:Nd;

for i=1:N
    for j=[1:(i-1),(i+1):N]
        diff=fighting_abilities(i)-fighting_abilities(j);
        d=exp(diff)/(exp(diff)+1);
        d=floor(round(d/deps))*deps;
        dmat(i,j)=min(v(abs(d-dvals)<=deps+0.001));
    end
end

Lvals=zeros(N,its);
Lvals(:,1)=11*ones(N,1);
% Lvals(:,1)=randi([1 NLx],N,1);
perfvals=zeros(N,its);

for i=1:N
    Lys=Lvals([1:(i-1),(i+1):N]);
    ds=dmat(i,[1:(i-1),(i+1):N]);
            times=zeros(1,N-1);
            probs=zeros(1,N-1);
            for q=1:(N-1)
                times(q)=f(Lvals(i),Lys(q),ds(q),3);
                probs(q)=f(Lvals(i),Lys(q),ds(q),1);
            end
        perf=sum(-times+probs.*(1-times));
    perfvals(i,1)=perf;
end

count=1;
while count<=its
    for i=1:N
        Lys=Lvals([1:(i-1),(i+1):N]);
        ds=dmat(i,[1:(i-1),(i+1):N]);
        perf=zeros(1,NLx);
        for j=1:NLx
            times=zeros(1,N-1);
            probs=zeros(1,N-1);
            for q=1:(N-1)
                times(q)=f(j,Lys(q),ds(q),3);
                probs(q)=f(j,Lys(q),ds(q),1);
            end
            perf(j)=sum(-times+probs.*(1-times));
        end
        [p,n]=max(perf);
        Lvals(i,count+1)=n;
        perfvals(i,count+1)=p;
    end
    if sum(Lvals(:,count+1)==Lvals(:,count))==N
        maxit=count+2;
        count=its+1;
    end
    count=count+1;
end

Lvals(:,maxit:end)=[];
perfvals(:,maxit:end)=[];

probmat=zeros(N,N);
timemat=zeros(N,N);

for i=1:N
    for j=[1:(i-1),(i+1):N]
        probmat(i,j)=f(Lvals(i,end),Lvals(j,end),dmat(i,j),1);
        timemat(i,j)=f(Lvals(i,end),Lvals(j,end),dmat(i,j),3);
    end
end

%%
figure
hist(fighting_abilities)
xlabel('Fighting Ability')
print('/Users/eleanorbrush/Dropbox/signaling_network/ability_histogram.eps','-depsc','-r300')

figure
hold on
for i=1:Nd
    plot(dvals(i),-NEx{i},'or',dvals(i),-NEy{i},'ob')
end
xlabel('Dominance');
ylabel('Nash Equilibrium Threshold');
leg=legend('focal','other');
set(leg,'Position',[.75 .5 .1 .2])
box on
print('/Users/eleanorbrush/Dropbox/signaling_network/pair_nasheq_threshold.eps','-depsc','-r300')


figure
plot(fighting_abilities,-Lxvals(Lvals(:,end)))
box on
xlabel('Ability')
ylabel('Optimal Threshold in Group Setting')
print('/Users/eleanorbrush/Dropbox/signaling_network/group_nasheq_threshold.eps','-depsc','-r300')

figure
imagesc(probmat)
colorbar
print('/Users/eleanorbrush/Dropbox/signaling_network/prob_mat.eps','-depsc','-r300')

figure
imagesc(timemat)
colorbar
print('/Users/eleanorbrush/Dropbox/signaling_network/time_mat.eps','-depsc','-r300')
