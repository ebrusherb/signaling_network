function toreturn=group_props_parallelized(c1,c2,c3,its,N,distribution,Nd,Nt,twomat,domvals,threshvals,l)

if nargin==3
    N=30;
    its=25;
end

Nb1=10;
Nb2=10;
% l=1;
weight=.05;

keeptrackofinfoD=zeros(its,Nb1*Nb2);
keeptrackofinfoR=zeros(its,Nb1*Nb2);
keeptrackofinfoDelta=zeros(its,Nb1*Nb2);
keeptrackofinfoH=zeros(its,Nb1*Nb2);
keeptrackofinfoPi=zeros(its,Nb1*Nb2);
keeptrackofinfoC=zeros(its,Nb1*Nb2);
keeptrackofinfoprobs=zeros(its,Nb1*Nb2);
infomatD=zeros(Nb1,Nb2);
infomatR=zeros(Nb1,Nb2);
infomatDelta=zeros(Nb1,Nb2);
infomatH=zeros(Nb1,Nb2);
infomatPi=zeros(Nb1,Nb2);
infomatC=zeros(Nb1,Nb2);
infomatprobs=zeros(Nb1,Nb2);

%powerfuns={'D','R','Delta','H','Pi''C','probs'};

skewvalsD=zeros(1,its);
skewvalsR=zeros(1,its);
skewvalsDelta=zeros(1,its);
skewvalsH=zeros(1,its);
skewvalsPi=zeros(1,its);
skewvalsprobs=zeros(1,its);
skewvalsC=zeros(1,its);
finalthresh=zeros(N,its);

keepprobmat=zeros(its,N*N);
keeptimemat=zeros(its,N*N);
keepsigmat=zeros(its,N*N);

pairwiseaccuracy=zeros(1,its);
pairwisetime=zeros(1,its);

for c=1:its

    opt_its=10;

    fighting_abilities=zeros(1,N);
    
    switch distribution
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

    % find differences in abilities --> dominance
    for i=1:N
        for j=[1:(i-1),(i+1):N]
            abilitydiff=fighting_abilities(i)-fighting_abilities(j);
            d=exp(abilitydiff)/(exp(abilitydiff)+1);
            d=floor(round(d/deps))*deps;
            dmat(i,j)=v(abs(d-extendeddomvals)<=deps/2);
        end
    end

    %size(twomat)=[Nl,Nd,Nt,Nt,2];
    
    perf=zeros(2,Nd,Nt,Nt); %individual, leak, dominance, thresholds

    perf(1,2:Nd,:,:)=c1*(1-twomat(l,2:Nd,:,:,1))+c2*(twomat(l,2:Nd,:,:,2))+c3*(1-twomat(l,2:Nd,:,:,1));
    perf(2,2:Nd,:,:)=c1*(1-twomat(l,2:Nd,:,:,1))+c2*(twomat(l,2:Nd,:,:,2))+c3*(twomat(l,2:Nd,:,:,1));

    perf(1,1,:,:)=c2*(twomat(l,1,:,:,2))+c3*(1-twomat(l,1,:,:,1));
    perf(2,1,:,:)=c2*(twomat(l,1,:,:,2))+c3*(twomat(l,1,:,:,1));

    Tvals=zeros(N,opt_its);
    Tvals(:,1)=2*ones(N,1);
    % Tvals(:,1)=randi([1 Nt],N,1);
     
    %find optimal thresholds
    for i=1:N
        opp_thresh=Tvals([1:(i-1),(i+1):N],1);
        ds=dmat(i,[1:(i-1),(i+1):N]);
        perfsum=0;
                for q=1:(N-1)
                    if ds(q)<=Nd-1
                        perfsum=perfsum+perf(2,Nd-ds(q)+1,opp_thresh(q),Tvals(i,1)); 
                    else 
                        perfsum=perfsum+perf(1,ds(q)-Nd+1,Tvals(i,1),opp_thresh(q));
                    end
                end  
        
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
                        perfsum=perfsum+perf(2,Nd-ds(q)+1,opp_thresh(q),j);
                    else 
                        perfsum=perfsum+perf(1,ds(q)-Nd+1,j,opp_thresh(q));
                    end
                end 
                perftest(j)=perfsum;
            end
            [~,n]=min(perftest);
            Tvals(i,count+1)=n;
            
        end
        if sum(Tvals(:,count+1)==Tvals(:,count))==N
            maxit=count+2;
            count=opt_its+1;
            Tvals(:,maxit:end)=[];
            
        end
        count=count+1;
    end
    
    finalthresh(:,c)=threshvals(Tvals(:,end));

    probmat=zeros(N,N);
    timemat=zeros(N,N);

    %size(twomat)=[Nl,Nd,Nt,Nt,2];
    
    %find prob and time between each pair
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
    
    binsigmat=zeros(N,N);
    
    %construct signaling network
    for i=1:N
        for j=(i+1):N
            draw=rand;
            if draw <=probmat(i,j)
                binsigmat(i,j)=1;
            else binsigmat(j,i)=1;
            end
        end
    end
    
    pairwiseaccuracy(c)=sum(sum(triu(binsigmat,1)));
    pairwisetime(c)=sum(sum(triu(timemat,1)));
    
    sigmat=binsigmat.*(max(max(timemat))+1-timemat);
    
    todivide=sum(sigmat,2);
    todivide(todivide==0)=1;
    entmat=sigmat./repmat(todivide,1,N);
    logmat=entmat;
    logmat(entmat~=0)=-log(entmat(entmat~=0));
    entmat=entmat.*logmat;

%         powervec=sum(probmat,2);
    powervecD=sum(binsigmat,2);
    powervecR=sum(sigmat,2);
    powervecDelta=(powervecD).*(powervecR);
    powervecH=sum(entmat,2);
    powervecH=real(powervecH);
    powervecPi=(powervecR).*(powervecH);
    powervecC=eigcentrality(sigmat,weight);
    powervecprobs=sum(probmat,2);
    
%     abilitiesdiff=ceil(max(fighting_abilities))-floor(min(fighting_abilities));
%     abilitiesbins=floor(min(fighting_abilities)):abilitiesdiff/Nb1:ceil(max(fighting_abilities));
    q=quantile(fighting_abilities,0:1/Nb1:1);
    diffvec=diff(q);
    eps=min(diffvec(diffvec>=1e-5))/10;
    q(1)=q(1)-eps;q(end)=q(end)+eps;
    abilitiesbins=q;
    [~,abilitiesindices]=histc(fighting_abilities,abilitiesbins);
    
    powervec=powervecD;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    if sum(diffvec)~=0
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);
    else powerindices=ceil(Nb2/2)*ones(N,1);
    end
    
    v=zeros(Nb1*Nb2,1);
    for i=1:N
        ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
        v(ind)=v(ind)+1;
    end
    keeptrackofinfoD(c,:)=v;
    
    skewvalsD(c)=skewness(powervec);
    
    powervec=powervecR;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    if sum(diffvec)~=0
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);
    else powerindices=ceil(Nb2/2)*ones(N,1);
    end
    
    v=zeros(Nb1*Nb2,1);
    for i=1:N
        ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
        v(ind)=v(ind)+1;
    end
    keeptrackofinfoR(c,:)=v;
    
    skewvalsR(c)=skewness(powervec);
    
    powervec=powervecDelta;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    if sum(diffvec)~=0
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);
    else powerindices=ceil(Nb2/2)*ones(N,1);
    end
    
    v=zeros(Nb1*Nb2,1);
    for i=1:N
        ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
        v(ind)=v(ind)+1;
    end
    keeptrackofinfoDelta(c,:)=v;
    
    skewvalsDelta(c)=skewness(powervec);
    
    powervec=powervecH;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    if sum(diffvec)~=0
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);
    else powerindices=ceil(Nb2/2)*ones(N,1);
    end
    
    v=zeros(Nb1*Nb2,1);
    for i=1:N
        ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
        v(ind)=v(ind)+1;
    end
    keeptrackofinfoH(c,:)=v;
    
    skewvalsH(c)=skewness(powervec);
    
    powervec=powervecPi;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    if sum(diffvec)~=0
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);
    else powerindices=ceil(Nb2/2)*ones(N,1);
    end
    
    v=zeros(Nb1*Nb2,1);
    for i=1:N
        ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
        v(ind)=v(ind)+1;
    end
    keeptrackofinfoPi(c,:)=v;
    
    skewvalsPi(c)=skewness(powervec);
    
    powervec=powervecC;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    if sum(diffvec)~=0
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);
    else powerindices=ceil(Nb2/2)*ones(N,1);
    end
    
    v=zeros(Nb1*Nb2,1);
    for i=1:N
        ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
        v(ind)=v(ind)+1;
    end
    keeptrackofinfoC(c,:)=v;
    
    skewvalsC(c)=skewness(powervec);
    
    powervec=powervecprobs;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    if sum(diffvec)~=0
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        powerbins=q;
        [~,powerindices]=histc(powervec,powerbins);
    else powerindices=ceil(Nb2/2)*ones(N,1);
    end
    
    v=zeros(Nb1*Nb2,1);
    for i=1:N
        ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
        v(ind)=v(ind)+1;
    end
    keeptrackofinfoprobs(c,:)=v;
    
    skewvalsprobs(c)=skewness(powervec);
 
    keepprobmat(c,:)=reshape(probmat,1,[]);
    keeptimemat(c,:)=reshape(timemat,1,[]);
    keepsigmat(c,:)=reshape(binsigmat,1,[]);
    
end

keeptrackofinfoD=sum(keeptrackofinfoD,1);
keeptrackofinfoR=sum(keeptrackofinfoR,1);
keeptrackofinfoDelta=sum(keeptrackofinfoDelta,1);
keeptrackofinfoH=sum(keeptrackofinfoH,1);
keeptrackofinfoPi=sum(keeptrackofinfoPi,1);
keeptrackofinfoC=sum(keeptrackofinfoC,1);
keeptrackofinfoprobs=sum(keeptrackofinfoprobs,1);

infomatD=reshape(keeptrackofinfoD,Nb1,Nb2);
infomatR=reshape(keeptrackofinfoR,Nb1,Nb2);
infomatDelta=reshape(keeptrackofinfoDelta,Nb1,Nb2);
infomatH=reshape(keeptrackofinfoH,Nb1,Nb2);
infomatPi=reshape(keeptrackofinfoPi,Nb1,Nb2);
infomatC=reshape(keeptrackofinfoC,Nb1,Nb2);
infomatprobs=reshape(keeptrackofinfoprobs,Nb1,Nb2);

avgthresh=mean(finalthresh,2);
stdthresh=std(finalthresh,0,2);

probmat=reshape(keepprobmat(its,:),N,N);
timemat=reshape(keeptimemat(its,:),N,N);
binsigmat=reshape(keepsigmat(its,:),N,N);
totalsigmat=reshape(sum(keepsigmat,1),N,N);

sigmat=binsigmat.*(max(max(timemat))+1-timemat);
    
todivide=sum(sigmat,2);
todivide(todivide==0)=1;
entmat=sigmat./repmat(todivide,1,N);
logmat=entmat;
logmat(entmat~=0)=-log(entmat(entmat~=0));
entmat=entmat.*logmat;

%         powervec=sum(probmat,2);
powervecD=sum(binsigmat,2);
powervecR=sum(sigmat,2);
powervecDelta=(powervecD).*(powervecR);
powervecH=sum(entmat,2);
powervecH=real(powervecH);
powervecPi=(powervecR).*(powervecH);
powervecC=eigcentrality(sigmat,weight);
powervecprobs=sum(probmat,2);

pairwiseaccuracy=sum(pairwiseaccuracy)/its/(N*(N-1)/2);
pairwisetime=sum(pairwisetime)/its/(N*(N-1)/2);

groupmutinfoD=0;

for i=1:Nb1
    px=sum(infomatD(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomatD(:,j))/(N*its);
        pxy=infomatD(i,j)/(N*its);
        if pxy~=0
            groupmutinfoD=groupmutinfoD+pxy*log(pxy/px/py);
        end
    end
end

groupmutinfoR=0;

for i=1:Nb1
    px=sum(infomatR(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomatR(:,j))/(N*its);
        pxy=infomatR(i,j)/(N*its);
        if pxy~=0
            groupmutinfoR=groupmutinfoR+pxy*log(pxy/px/py);
        end
    end
end

groupmutinfoDelta=0;

for i=1:Nb1
    px=sum(infomatDelta(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomatDelta(:,j))/(N*its);
        pxy=infomatDelta(i,j)/(N*its);
        if pxy~=0
            groupmutinfoDelta=groupmutinfoDelta+pxy*log(pxy/px/py);
        end
    end
end

groupmutinfoH=0;

for i=1:Nb1
    px=sum(infomatH(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomatH(:,j))/(N*its);
        pxy=infomatH(i,j)/(N*its);
        if pxy~=0
            groupmutinfoH=groupmutinfoH+pxy*log(pxy/px/py);
        end
    end
end

groupmutinfoPi=0;

for i=1:Nb1
    px=sum(infomatPi(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomatPi(:,j))/(N*its);
        pxy=infomatPi(i,j)/(N*its);
        if pxy~=0
            groupmutinfoPi=groupmutinfoPi+pxy*log(pxy/px/py);
        end
    end
end

groupmutinfoC=0;

for i=1:Nb1
    px=sum(infomatC(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomatC(:,j))/(N*its);
        pxy=infomatC(i,j)/(N*its);
        if pxy~=0
            groupmutinfoC=groupmutinfoC+pxy*log(pxy/px/py);
        end
    end
end

groupmutinfoprobs=0;
for i=1:Nb1
    px=sum(infomatprobs(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomatprobs(:,j))/(N*its);
        pxy=infomatprobs(i,j)/(N*its);
        if pxy~=0
            groupmutinfoprobs=groupmutinfoprobs+pxy*log(pxy/px/py);
        end
    end
end

meanskewnessD=mean(skewvalsD);
meanskewnessR=mean(skewvalsR);
meanskewnessDelta=mean(skewvalsDelta);
meanskewnessH=mean(skewvalsH);
meanskewnessPi=mean(skewvalsPi);
meanskewnessC=mean(skewvalsC);
meanskewnessprobs=mean(skewvalsprobs);

powervecarray=[powervecD,powervecR,powervecDelta,powervecH,powervecPi,powervecC,powervecprobs];
groupmutarray=[groupmutinfoD,groupmutinfoR,groupmutinfoDelta,groupmutinfoH,groupmutinfoPi,groupmutinfoC,groupmutinfoprobs];
meanskewnessarray=[meanskewnessD,meanskewnessR,meanskewnessDelta,meanskewnessH,meanskewnessPi,meanskewnessC,meanskewnessprobs];
infomatlist={infomatD,infomatR,infomatDelta,infomatH,infomatPi,infomatC};

toreturn={powervecarray,groupmutarray,meanskewnessarray,pairwiseaccuracy,pairwisetime,infomatlist,avgthresh,stdthresh,timemat,probmat,totalsigmat};
end