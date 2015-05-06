%% coordinates
load variables.mat
load twomat.mat

N=20;
its=1000;
cstep=.1;
vec=0:cstep:1;
Nc=length(vec);
endvals=zeros(1,Nc);

% X=repmat(0:cstep:1,1,11);
% Y=reshape(repmat(0:cstep:1,11,1),1,[]);

X=[];
Y=[];
layers=cell(Nc,1);
for i=0:(Nc-1)
    layers{i+1}=length(X)+(1:(Nc-i));
    endvals(i+1)=layers{i+1}(end);
    X=[X, vec(1:(Nc-i))];
    Y=[Y,vec(i+1)*ones(1,Nc-i)];   
end
Z=1-X-Y;

transformed=transform([X',Z',Y']);
Nc=length(X);

indices=1:Nt;
leak=2;
powerfuns={'D','R','Delta','H','Pi','C','probs'};
numpowerfuns=length(powerfuns);

dist='unif';
ivals=[1 11 12 39 52 57];
for i=1:length(ivals)
    cvecs(i,:)=[X(ivals(i)), Y(ivals(i)), Z(ivals(i))];
end

Nb1=10;
Nb2=10;
l=1;
weight=.05;

Nc=size(cvecs,1);
%% finding info mats
% function toreturn=group_props_parallelized(cvecs,its,N,dist,Nd,Nt,twomat,domvals,threshvals)
% 
% if nargin==3
%     N=30;
%     its=25;
% end

keeptrackofinfoD=zeros(its,Nc,Nb1*Nb2);
keeptrackofinfoR=zeros(its,Nc,Nb1*Nb2);
keeptrackofinfoDelta=zeros(its,Nc,Nb1*Nb2);
keeptrackofinfoH=zeros(its,Nc,Nb1*Nb2);
keeptrackofinfoPi=zeros(its,Nc,Nb1*Nb2);
keeptrackofinfoC=zeros(its,Nc,Nb1*Nb2);
keeptrackofinfoprobs=zeros(its,Nc,Nb1*Nb2);
infomatD=zeros(Nc,Nb1,Nb2);
infomatR=zeros(Nc,Nb1,Nb2);
infomatDelta=zeros(Nc,Nb1,Nb2);
infomatH=zeros(Nc,Nb1,Nb2);
infomatPi=zeros(Nc,Nb1,Nb2);
infomatC=zeros(Nc,Nb1,Nb2);
infomatprobs=zeros(Nc,Nb1,Nb2);
powervecD=zeros(its,Nc,N);
powervecR=zeros(its,Nc,N);
powervecDelta=zeros(its,Nc,N);
powervecH=zeros(its,Nc,N);
powervecPi=zeros(its,Nc,N);
powervecC=zeros(its,Nc,N);
powervecprobs=zeros(its,Nc,N);

skewvalsD=zeros(its,Nc);
skewvalsR=zeros(its,Nc);
skewvalsDelta=zeros(its,Nc);
skewvalsH=zeros(its,Nc);
skewvalsPi=zeros(its,Nc);
skewvalsprobs=zeros(its,Nc);
skewvalsC=zeros(its,Nc);

finalthresh=zeros(its,Nc,N);

keepprobmat=zeros(its,Nc,N*N);
keeptimemat=zeros(its,Nc,N*N);
keepsigmat=zeros(its,Nc,N*N);

pairwiseaccuracy=zeros(its,Nc);
pairwisetime=zeros(its,Nc);

for c=1:its

    opt_its=10;

    fighting_abilities=zeros(1,N);
    
    switch dist
        case 'unif'
            fighting_abilities=20*rand(1,N);
        case 'norm'
            fighting_abilities=normrnd(0,10,1,N);
        case 'lognorm'
            fighting_abilities=lognrnd(0,1,1,N);
    end
    
    fighting_abilities=sort(fighting_abilities,'descend');

    dmat=zeros(N); %turn difference in fighting abilities into probabilitise of winning
    extendeddomvals=[1-domvals(end:-1:2) domvals];
    v=1:(Nd+Nd-1);
    deps=.05;

    for i=1:N
        for j=[1:(i-1),(i+1):N]
            abilitydiff=fighting_abilities(i)-fighting_abilities(j);
            d=exp(abilitydiff)/(exp(abilitydiff)+1);
            d=floor(round(d/deps))*deps;
            dmat(i,j)=v(abs(d-extendeddomvals)<=deps/2);
        end
    end
    
    for ic=1:Nc
        Tvals=abilities2thresh(cvecs(ic,:),twomat,dmat,opt_its,threshvals);
        finalthresh(c,ic,:)=threshvals(Tvals);
        probmat=zeros(N,N);
        timemat=zeros(N,N);

        %size(twomat)=[Nl,Nd,Nt,Nt,2];

        for i=1:N
            for j=[1:(i-1),(i+1):N]
                if dmat(i,j)<=Nd-1
                    probmat(i,j)=1-twomat(l,Nd-dmat(i,j)+1,Tvals(j),Tvals(i),1);
                    timemat(i,j)=twomat(l,Nd-dmat(i,j)+1,Tvals(j),Tvals(i),2);
                else
                    probmat(i,j)=twomat(l,dmat(i,j)-Nd+1,Tvals(i),Tvals(j),1);
                    timemat(i,j)=twomat(l,dmat(i,j)-Nd+1,Tvals(i),Tvals(j),2);
                end
            end
        end

        binsigmat=zeros(N,N);

        for i=1:N
            for j=(i+1):N
                draw=rand;
                if draw <=probmat(i,j)
                    binsigmat(i,j)=1;
                else binsigmat(j,i)=1;
                end
            end
        end

        pairwiseaccuracy(c,ic)=sum(sum(triu(binsigmat,1)));
        pairwisetime(c,ic)=sum(sum(triu(timemat,1)));

        sigmat=binsigmat.*(max(max(timemat))+1-timemat);

        todivide=sum(sigmat,2);
        todivide(todivide==0)=1;
        entmat=sigmat./repmat(todivide,1,N);
        logmat=entmat;
        logmat(entmat~=0)=-log(entmat(entmat~=0));
        entmat=entmat.*logmat;

    %         powervec=sum(probmat,2);
        powervecD(c,ic,:)=sum(binsigmat,2);
        powervecR(c,ic,:)=sum(sigmat,2);
        powervecDelta(c,ic,:)=(powervecD(c,ic,:)).*(powervecR(c,ic,:));
        powervecH(c,ic,:)=sum(entmat,2);
        powervecH(c,ic,:)=real(powervecH(c,ic,:));
        powervecPi(c,ic,:)=(powervecR(c,ic,:)).*(powervecH(c,ic,:));
        powervecC(c,ic,:)=eigcentrality(sigmat,weight);
        powervecprobs(c,ic,:)=sum(probmat,2);

%         abilitiesdiff=ceil(max(fighting_abilities))-floor(min(fighting_abilities));
%         abilitiesbins=floor(min(fighting_abilities)):abilitiesdiff/Nb1:ceil(max(fighting_abilities));
        q=quantile(fighting_abilities,0:1/Nb1:1);
        diffvec=diff(q);
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        abilitiesbins=q;
        [~,abilitiesindices]=histc(fighting_abilities,abilitiesbins);

        powervec=powervecD(c,ic,:);
        powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
        q=quantile(powervec,0:1/Nb2:1);
        diffvec=diff(q);
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);

        v=zeros(Nb1*Nb2,1);
        for i=1:N
            ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
            v(ind)=v(ind)+1;
        end
        keeptrackofinfoD(c,ic,:)=v;

        skewvalsD(c,ic)=skewness(powervec);

        powervec=powervecR(c,ic,:);
    %     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
        q=quantile(powervec,0:1/Nb2:1);
        diffvec=diff(q);
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);

        v=zeros(Nb1*Nb2,1);
        for i=1:N
            ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
            v(ind)=v(ind)+1;
        end
        keeptrackofinfoR(c,ic,:)=v;

        skewvalsR(c,ic)=skewness(powervec);

        powervec=powervecDelta(c,ic,:);
    %     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
        q=quantile(powervec,0:1/Nb2:1);
        diffvec=diff(q);
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);

        v=zeros(Nb1*Nb2,1);
        for i=1:N
            ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
            v(ind)=v(ind)+1;
        end
        keeptrackofinfoDelta(c,ic,:)=v;

        skewvalsDelta(c,ic)=skewness(powervec);

        powervec=powervecH(c,ic,:);
    %     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
        q=quantile(powervec,0:1/Nb2:1);
        diffvec=diff(q);
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);

        v=zeros(Nb1*Nb2,1);
        for i=1:N
            ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
            v(ind)=v(ind)+1;
        end
        keeptrackofinfoH(c,ic,:)=v;

        skewvalsH(c,ic)=skewness(powervec);

        powervec=powervecPi(c,ic,:);
    %     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
        q=quantile(powervec,0:1/Nb2:1);
        diffvec=diff(q);
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);

        v=zeros(Nb1*Nb2,1);
        for i=1:N
            ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
            v(ind)=v(ind)+1;
        end
        keeptrackofinfoPi(c,ic,:)=v;

        skewvalsPi(c,ic)=skewness(powervec);

        powervec=powervecC(c,ic,:);
    %     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
        q=quantile(powervec,0:1/Nb2:1);
        diffvec=diff(q);
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);

        v=zeros(Nb1*Nb2,1);
        for i=1:N
            ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
            v(ind)=v(ind)+1;
        end
        keeptrackofinfoC(c,ic,:)=v;

        skewvalsC(c,ic)=skewness(powervec);

        powervec=powervecprobs(c,ic,:);
    %     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
        q=quantile(powervec,0:1/Nb2:1);
        diffvec=diff(q);
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);

        v=zeros(Nb1*Nb2,1);
        for i=1:N
            ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
            v(ind)=v(ind)+1;
        end
        keeptrackofinfoprobs(c,ic,:)=v;

        skewvalsprobs(c,ic)=skewness(powervec);

        keepprobmat(c,ic,:)=reshape(probmat,1,[]);
        keeptimemat(c,ic,:)=reshape(timemat,1,[]);
        keepsigmat(c,ic,:)=reshape(sigmat,1,[]);
    end
end

summedinfoD=reshape(sum(keeptrackofinfoD,1),Nc,[]);
summedinfoR=reshape(sum(keeptrackofinfoR,1),Nc,[]);
summedinfoDelta=reshape(sum(keeptrackofinfoDelta,1),Nc,[]);
summedinfoH=reshape(sum(keeptrackofinfoH,1),Nc,[]);
summedinfoPi=reshape(sum(keeptrackofinfoPi,1),Nc,[]);
summedinfoC=reshape(sum(keeptrackofinfoC,1),Nc,[]);
summedinfoprobs=reshape(sum(keeptrackofinfoprobs,1),Nc,[]);

for ic=1:Nc
    for i=1:Nb1
        for j=1:Nb2
            ind=sub2ind([Nb1,Nb2],i,j);
            infomatD(ic,i,j)=summedinfoD(ic,ind);
        end
    end
    for i=1:Nb1
        for j=1:Nb2
            ind=sub2ind([Nb1,Nb2],i,j);
            infomatR(ic,i,j)=summedinfoR(ic,ind);
        end
    end
    for i=1:Nb1
        for j=1:Nb2
            ind=sub2ind([Nb1,Nb2],i,j);
            infomatDelta(ic,i,j)=summedinfoDelta(ic,ind);
        end
    end
    for i=1:Nb1
        for j=1:Nb2
            ind=sub2ind([Nb1,Nb2],i,j);
            infomatH(ic,i,j)=summedinfoH(ic,ind);
        end
    end
    for i=1:Nb1
        for j=1:Nb2
            ind=sub2ind([Nb1,Nb2],i,j);
            infomatPi(ic,i,j)=summedinfoPi(ic,ind);
        end
    end
    for i=1:Nb1
        for j=1:Nb2
            ind=sub2ind([Nb1,Nb2],i,j);
            infomatC(ic,i,j)=summedinfoC(ic,ind);
        end
    end
    for i=1:Nb1
        for j=1:Nb2
            ind=sub2ind([Nb1,Nb2],i,j);
            infomatprobs(ic,i,j)=summedinfoprobs(ic,ind);
        end
    end
end

pairwiseaccuracy=sum(pairwiseaccuracy,1)/its/(N*(N-1)/2);
pairwisetime=sum(pairwisetime,1)/its/(N*(N-1)/2);

% probmat=zeros(ic,N,N);
% timemat=zeros(ic,N,N);
% binsigmat=zeros(ic,N,N);
% sigmat=zeros(ic,N,N);
% powervecD=zeros(ic,N);
% powervecR=zeros(ic,N);
% powervecDelta=zeros(ic,N);
% powervecH=zeros(ic,N);
% powervecPi=zeros(ic,N);
% powervecC=zeros(ic,N);
% powervecprobs=zeros(ic,N);
% 
% for ic=1:Nc
%     probmat(ic,:,:)=reshape(keepprobmat(its,ic,:),N,N);
%     timemat(ic,:,:)=reshape(keeptimemat(its,ic,:),N,N);
%     binsigmat(ic,:,:)=reshape(keepsigmat(its,ic,:),N,N);
% 
%     sigmat(ic,:,:)=binsigmat(ic,:,:).*(max(max(timemat(ic,:,:)))+1-timemat(ic,:,:));
% 
%     todivide=sum(reshape(sigmat(ic,:,:),N,N),2);
%     todivide(todivide==0)=1;
%     entmat=reshape(sigmat(ic,:,:),N,N)./repmat(todivide,1,N);
%     logmat=entmat;
%     logmat(entmat~=0)=-log(entmat(entmat~=0));
%     entmat=entmat.*logmat;
% 
%     %         powervec=sum(probmat,2);
%     powervecD(ic,:)=sum(reshape(binsigmat(ic,:,:),N,N),2);
%     powervecR(ic,:)=sum(reshape(sigmat(ic,:,:),N,N),2);
%     powervecDelta(ic,:)=(powervecD(ic,:)).*(powervecR(ic,:));
%     powervecH(ic,:)=sum(entmat,2);
%     powervecH(ic,:)=real(powervecH(ic,:));
%     powervecPi(ic,:)=(powervecR(ic,:)).*(powervecH(ic,:));
%     powervecC(ic,:)=eigcentrality(reshape(sigmat(ic,:,:),N,N),weight);
%     powervecprobs(ic,:)=sum(reshape(probmat(ic,:,:),N,N),2);
% end

groupmutinfoD=zeros(Nc,1);
groupmutinfoR=zeros(Nc,1);
groupmutinfoDelta=zeros(Nc,1);
groupmutinfoH=zeros(Nc,1);
groupmutinfoPi=zeros(Nc,1);
groupmutinfoC=zeros(Nc,1);
groupmutinfoprobs=zeros(Nc,1);

for ic=1:Nc
    infomat=reshape(infomatD(ic,:,:),Nb1,Nb2);
    for i=1:Nb1
        px=sum(infomat(i,:))/(N*its);
        for j=1:Nb2
            py=sum(infomat(:,j))/(N*its);
            pxy=infomat(i,j)/(N*its);
            if pxy~=0
                groupmutinfoD(ic)=groupmutinfoD(ic)+pxy*log(pxy/px/py);
            end
        end
    end

    infomat=reshape(infomatR(ic,:,:),Nb1,Nb2);
    for i=1:Nb1
        px=sum(infomat(i,:))/(N*its);
        for j=1:Nb2
            py=sum(infomat(:,j))/(N*its);
            pxy=infomat(i,j)/(N*its);
            if pxy~=0
                groupmutinfoR(ic)=groupmutinfoR(ic)+pxy*log(pxy/px/py);
            end
        end
    end

    infomat=reshape(infomatDelta(ic,:,:),Nb1,Nb2);
    for i=1:Nb1
        px=sum(infomat(i,:))/(N*its);
        for j=1:Nb2
            py=sum(infomat(:,j))/(N*its);
            pxy=infomat(i,j)/(N*its);
            if pxy~=0
                groupmutinfoDelta(ic)=groupmutinfoDelta(ic)+pxy*log(pxy/px/py);
            end
        end
    end


    infomat=reshape(infomatH(ic,:,:),Nb1,Nb2);
    for i=1:Nb1
        px=sum(infomat(i,:))/(N*its);
        for j=1:Nb2
            py=sum(infomat(:,j))/(N*its);
            pxy=infomat(i,j)/(N*its);
            if pxy~=0
                groupmutinfoH(ic)=groupmutinfoH(ic)+pxy*log(pxy/px/py);
            end
        end
    end


    infomat=reshape(infomatPi(ic,:,:),Nb1,Nb2);
    for i=1:Nb1
        px=sum(infomat(i,:))/(N*its);
        for j=1:Nb2
            py=sum(infomat(:,j))/(N*its);
            pxy=infomat(i,j)/(N*its);
            if pxy~=0
                groupmutinfoPi(ic)=groupmutinfoPi(ic)+pxy*log(pxy/px/py);
            end
        end
    end


    infomat=reshape(infomatC(ic,:,:),Nb1,Nb2);
    for i=1:Nb1
        px=sum(infomat(i,:))/(N*its);
        for j=1:Nb2
            py=sum(infomat(:,j))/(N*its);
            pxy=infomat(i,j)/(N*its);
            if pxy~=0
                groupmutinfoC(ic)=groupmutinfoC(ic)+pxy*log(pxy/px/py);
            end
        end
    end

    infomat=reshape(infomatprobs(ic,:,:),Nb1,Nb2);
    for i=1:Nb1
        px=sum(infomat(i,:))/(N*its);
        for j=1:Nb2
            py=sum(infomat(:,j))/(N*its);
            pxy=infomat(i,j)/(N*its);
            if pxy~=0
                groupmutinfoprobs(ic)=groupmutinfoprobs(ic)+pxy*log(pxy/px/py);
            end
        end
    end
end
meanskewnessD=mean(skewvalsD,1);
meanskewnessR=mean(skewvalsR,1);
meanskewnessDelta=mean(skewvalsDelta,1);
meanskewnessH=mean(skewvalsH,1);
meanskewnessPi=mean(skewvalsPi,1);
meanskewnessC=mean(skewvalsC,1);
meanskewnessprobs=mean(skewvalsprobs,1);

toreturn={[powervecD,powervecR,powervecDelta,powervecH,powervecPi,powervecC,powervecprobs],[groupmutinfoD,groupmutinfoR,groupmutinfoDelta,groupmutinfoH,groupmutinfoPi,groupmutinfoC,groupmutinfoprobs],[meanskewnessD,meanskewnessR,meanskewnessDelta,meanskewnessH,meanskewnessPi,meanskewnessC,meanskewnessprobs],pairwiseaccuracy,pairwisetime,{infomatD,infomatR,infomatDelta,infomatH,infomatPi,infomatC}};
% end
save('/Users/eleanorbrush/Desktop/group_props_comparison.mat','summedinfoD','summedinfoR','summedinfoC','summedinfoH','summedinfoPi','summedinfoDelta','groupmutinfoD','groupmutinfoR','groupmutinfoDelta','groupmutinfoH','groupmutinfoPi','groupmutinfoC')
%% plot info mats
load variables.mat
load /Users/eleanorbrush/Desktop/group_props_comparison.mat
divcolormap=cbrewer('div','RdBu',11);
interestingfuns=[1 2 4 6];
groupmutinfomat=zeros(Nc,length(interestingfuns));
for i=1:length(interestingfuns)
    groupmutinfomat(:,i)=eval(sprintf('groupmutinfo%s',powerfuns{interestingfuns(i)}));
end
bestfunvec=zeros(Nc,1);
for i=1:Nc
    [~,w]=max(groupmutinfomat(i,:));
    bestfunvec(i)=w;
end

for ic=2
    bestnummat=eval(sprintf('summedinfo%s',powerfuns{interestingfuns(bestfunvec(ic))}));
    bestnummat=bestnummat(ic,:);
    bestnummat=reshape(bestnummat,Nb1,Nb2);
    toplotmat=zeros(4,Nb1,Nb2);
    for i=1:4
        nummat=eval(sprintf('summedinfo%s',powerfuns{interestingfuns(i)}));
        nummat=nummat(ic,:);
        nummat=reshape(nummat,Nb1,Nb2);
%         toplotmat(i,:,:)=nummat2logprobmat(bestnummat)-nummat2logprobmat(nummat);
        toplotmat(i,:,:)=(bestnummat)-(nummat);
    end
    bottom=min(min(min(toplotmat)));
    top=max(max(max(toplotmat)));
    M=max(abs(bottom),abs(top));
%     w=max(abs(bottom),abs(top));
%     bottom=-w;
%     top=w; 
    figure
    for i=1:4
        subplot(2,2,i)   
        imagesc(reshape(toplotmat(i,:,:),Nb1,Nb2))
        caxis manual
        caxis([-M M])
        colormap(divcolormap)
        colorbar
        title(powerfuns{interestingfuns(i)},'FontSize',20)
    end
    set(gcf,'Position',[20+100*(ic-1) 378 560 420])
end

%%

seqcolormap=cbrewer('seq','Reds',10);
for ic=1:5;
    M=0;
figure
for i=1:4
nummat=eval(sprintf('summedinfo%s',powerfuns{interestingfuns(i)}));
% toplot=nummat2logprobmat(reshape(nummat(ic,:,:),Nb1,Nb2));
toplot=(reshape(nummat(ic,:,:),Nb1,Nb2));
M=max(M,max(max(toplot)));
subplot(2,2,i)
imagesc(toplot)
colormap(seqcolormap)
colorbar
title(powerfuns{interestingfuns(i)},'FontSize',20)
end
for i=1:4
    subplot(2,2,i)
    caxis manual
    caxis([0 M]);
end
set(gcf,'Position',[20+100*(ic-1) 20 560 420])
end

%also see test_networks.m
%%
plot(Y(ivals),mutvals(interestingfuns,ivals)','-o')
legend('D','R','H','C')

%%
figure
subplot(1,3,1)

imagesc(nummat2logprobmat(reshape(summedinfoR(1,:,:),Nb1,Nb2)))
colorbar
% colormap(seqcolormap)
caxis manual
caxis([0 .17])
subplot(1,3,2)
imagesc(nummat2logprobmat(reshape(summedinfoR(2,:,:),Nb1,Nb2)))
colorbar
% colormap(seqcolormap)
caxis manual
caxis([0 .17])
subplot(1,3,3)
imagesc(nummat2logprobmat(reshape(summedinfoR(2,:,:),Nb1,Nb2))-nummat2logprobmat(reshape(summedinfoR(1,:,:),Nb1,Nb2)))
colorbar
colormap(divcolormap)
caxis manual
m=min(min(nummat2logprobmat(reshape(summedinfoR(2,:,:),Nb1,Nb2))-nummat2logprobmat(reshape(summedinfoR(1,:,:),Nb1,Nb2))));
M=max(max(nummat2logprobmat(reshape(summedinfoR(2,:,:),Nb1,Nb2))-nummat2logprobmat(reshape(summedinfoR(1,:,:),Nb1,Nb2))));
M=max(abs(m),abs(M));
caxis([-M M])

figure
subplot(1,3,1)
imagesc(nummat2logprobmat(reshape(summedinfoC(1,:,:),Nb1,Nb2)))
colorbar
% colormap(seqcolormap)
caxis manual
caxis([0 .17])
subplot(1,3,2)
imagesc(nummat2logprobmat(reshape(summedinfoC(2,:,:),Nb1,Nb2)))
colorbar
% colormap(seqcolormap)
caxis manual
caxis([0 .17])
subplot(1,3,3)
imagesc(nummat2logprobmat(reshape(summedinfoC(2,:,:),Nb1,Nb2))-nummat2logprobmat(reshape(summedinfoC(1,:,:),Nb1,Nb2)))
colorbar
colormap(divcolormap)
caxis manual
m=min(min(nummat2logprobmat(reshape(summedinfoC(2,:,:),Nb1,Nb2))-nummat2logprobmat(reshape(summedinfoC(1,:,:),Nb1,Nb2))));
M=max(max(nummat2logprobmat(reshape(summedinfoC(2,:,:),Nb1,Nb2))-nummat2logprobmat(reshape(summedinfoC(1,:,:),Nb1,Nb2))));
M=max(abs(m),abs(M));
caxis([-M M])
