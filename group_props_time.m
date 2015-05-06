function toreturn=group_props_time(c1,c2,c3,its,N,distribution,Nd,Nt,twomat,domvals,threshvals,l,timevec)

if nargin==3
    N=30;
    its=25;
end

Nb1=10;
Nb2=10;
% l=1;

Ntime=length(timevec);

keeptrackofinfoD=zeros(its,Nb1*Nb2,Ntime);
keeptrackofinfoR=zeros(its,Nb1*Nb2,Ntime);
keeptrackofinfoDelta=zeros(its,Nb1*Nb2,Ntime);
keeptrackofinfoH=zeros(its,Nb1*Nb2,Ntime);
keeptrackofinfoPi=zeros(its,Nb1*Nb2,Ntime);
keeptrackofinfoC=zeros(its,Nb1*Nb2,Ntime);

% infomatD=zeros(Nb1,Nb2,Ntime);
% infomatR=zeros(Nb1,Nb2,Ntime);
% infomatDelta=zeros(Nb1,Nb2,Ntime);
% infomatH=zeros(Nb1,Nb2,Ntime);
% infomatPi=zeros(Nb1,Nb2,Ntime);
% infomatC=zeros(Nb1,Nb2,Ntime);

%powerfunlabs={'D','R','Delta','H','Pi','C','probs'};

skewvalsD=zeros(1,its);
skewvalsR=zeros(1,its);
skewvalsDelta=zeros(1,its);
skewvalsH=zeros(1,its);
skewvalsPi=zeros(1,its);
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
    
    %     abilitiesdiff=ceil(max(fighting_abilities))-floor(min(fighting_abilities));
%     abilitiesbins=floor(min(fighting_abilities)):abilitiesdiff/Nb1:ceil(max(fighting_abilities));
    q=quantile(fighting_abilities,0:1/Nb1:1);
    diffvec=diff(q);
    eps=min(diffvec(diffvec>=1e-5))/10;
    q(1)=q(1)-eps;q(end)=q(end)+eps;
    abilitiesbins=q;
    [~,abilitiesindices]=histc(fighting_abilities,abilitiesbins);
    
    for t=1:(Ntime-1)
        timematnow=zeros(N,N);
        
        timematnow(timemat<=timevec(t))=1;
        numsigsnow=zeros(N,N);
        for i=1:N
            for j=[1:(i-1),(i+1):N]
                if timemat(i,j)<=timevec(t)
                    numsigsnow(i,j)=(timevec(t)-timemat(i,j))/(max(max(timemat))+1-timemat(i,j))*sigmat(i,j); %linear growth of signals
                end
            end
        end
%         numsigsnow=timematnow;
        
        binsigmatnow=binsigmat.*timematnow;
        
        sigmatnow=binsigmatnow.*numsigsnow;
        
        [powervecDnow,powervecRnow,powervecDeltanow,powervecHnow,powervecPinow,powervecCnow] = powerfuns(binsigmatnow,sigmatnow);
        powervec=powervecDnow;
    %     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
        q=quantile(powervec,0:1/Nb2:1);
        diffvec=diff(q);
        if sigfig(sum(diffvec),4)~=0
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
        keeptrackofinfoD(c,:,t)=v;

        powervec=powervecRnow;
    %     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
        q=quantile(powervec,0:1/Nb2:1);
        diffvec=diff(q);
        if sigfig(sum(diffvec),4)~=0
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
        keeptrackofinfoR(c,:,t)=v;

        powervec=powervecDeltanow;
    %     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
        q=quantile(powervec,0:1/Nb2:1);
        diffvec=diff(q);
        if sigfig(sum(diffvec),4)~=0
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
        keeptrackofinfoDelta(c,:,t)=v;

        powervec=powervecHnow;
    %     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
        q=quantile(powervec,0:1/Nb2:1);
        diffvec=diff(q);
        if sigfig(sum(diffvec),4)~=0
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
        keeptrackofinfoH(c,:,t)=v;

        powervec=powervecPinow;
    %     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
        q=quantile(powervec,0:1/Nb2:1);
        diffvec=diff(q);
        if sigfig(sum(diffvec),4)~=0
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
        keeptrackofinfoPi(c,:,t)=v;

        powervec=powervecCnow;
    %     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
        q=quantile(powervec,0:1/Nb2:1);
        diffvec=diff(q);
        if sigfig(sum(diffvec),4)~=0
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
        keeptrackofinfoC(c,:,t)=v;
    end
    
    [powervecD,powervecR,powervecDelta,powervecH,powervecPi,powervecC] = powerfuns(binsigmat,sigmat);
   
    powervec=powervecD;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    if sigfig(sum(diffvec),4)~=0
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
    keeptrackofinfoD(c,:,end)=v;
    
    skewvalsD(c)=skewness(powervec);
    
    powervec=powervecR;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    if sigfig(sum(diffvec),4)~=0
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
    keeptrackofinfoR(c,:,end)=v;
    
    skewvalsR(c)=skewness(powervec);
    
    powervec=powervecDelta;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    if sigfig(sum(diffvec),4)~=0
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
    keeptrackofinfoDelta(c,:,end)=v;
    
    skewvalsDelta(c)=skewness(powervec);
    
    powervec=powervecH;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    if sigfig(sum(diffvec),4)~=0
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
    keeptrackofinfoH(c,:,end)=v;
    
    skewvalsH(c)=skewness(powervec);
    
    powervec=powervecPi;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    if sigfig(sum(diffvec),4)~=0
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
    keeptrackofinfoPi(c,:,end)=v;
    
    skewvalsPi(c)=skewness(powervec);
    
    powervec=powervecC;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    if sigfig(sum(diffvec),4)~=0
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
    keeptrackofinfoC(c,:,end)=v;
    
    skewvalsC(c)=skewness(powervec);
 
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

infomatD=reshape(keeptrackofinfoD,Nb1,Nb2,Ntime);
infomatR=reshape(keeptrackofinfoR,Nb1,Nb2,Ntime);
infomatDelta=reshape(keeptrackofinfoDelta,Nb1,Nb2,Ntime);
infomatH=reshape(keeptrackofinfoH,Nb1,Nb2,Ntime);
infomatPi=reshape(keeptrackofinfoPi,Nb1,Nb2,Ntime);
infomatC=reshape(keeptrackofinfoC,Nb1,Nb2,Ntime);

% avgthresh=mean(finalthresh,2);
% stdthresh=std(finalthresh,0,2);

% probmat=reshape(keepprobmat(its,:),N,N);
% timemat=reshape(keeptimemat(its,:),N,N);
% binsigmat=reshape(keepsigmat(its,:),N,N);
% totalsigmat=reshape(sum(keepsigmat,1),N,N);
% 
% sigmat=binsigmat.*(max(max(timemat))+1-timemat);
%    
% [powervecD,powervecR,powervecDelta,powervecH,powervecPi,powervecC] = powerfuns(binsigmat,sigmat);
% 
% pairwiseaccuracy=sum(pairwiseaccuracy)/its/(N*(N-1)/2);
% pairwisetime=sum(pairwisetime)/its/(N*(N-1)/2);

groupmutinfoD=zeros(1,Ntime);

for t=1:Ntime
    for i=1:Nb1
        px=sum(infomatD(i,:,t))/(N*its);
        for j=1:Nb2
            py=sum(infomatD(:,j,t))/(N*its);
            pxy=infomatD(i,j,t)/(N*its);
            if pxy~=0
                groupmutinfoD(t)=groupmutinfoD(t)+pxy*log(pxy/px/py);
            end
        end
    end
end

groupmutinfoR=zeros(1,Ntime);

for t=1:Ntime
    for i=1:Nb1
        px=sum(infomatR(i,:,t))/(N*its);
        for j=1:Nb2
            py=sum(infomatR(:,j,t))/(N*its);
            pxy=infomatR(i,j,t)/(N*its);
            if pxy~=0
                groupmutinfoR(t)=groupmutinfoR(t)+pxy*log(pxy/px/py);
            end
        end
    end
end

groupmutinfoDelta=zeros(1,Ntime);

for t=1:Ntime
    for i=1:Nb1
        px=sum(infomatDelta(i,:,t))/(N*its);
        for j=1:Nb2
            py=sum(infomatDelta(:,j,t))/(N*its);
            pxy=infomatDelta(i,j,t)/(N*its);
            if pxy~=0
                groupmutinfoDelta(t)=groupmutinfoDelta(t)+pxy*log(pxy/px/py);
            end
        end
    end
end

groupmutinfoH=zeros(1,Ntime);

for t=1:Ntime
    for i=1:Nb1
        px=sum(infomatH(i,:,t))/(N*its);
        for j=1:Nb2
            py=sum(infomatH(:,j,t))/(N*its);
            pxy=infomatH(i,j,t)/(N*its);
            if pxy~=0
                groupmutinfoH(t)=groupmutinfoH(t)+pxy*log(pxy/px/py);
            end
        end
    end
end

groupmutinfoPi=zeros(1,Ntime);

for t=1:Ntime
    for i=1:Nb1
        px=sum(infomatPi(i,:,t))/(N*its);
        for j=1:Nb2
            py=sum(infomatPi(:,j,t))/(N*its);
            pxy=infomatPi(i,j,t)/(N*its);
            if pxy~=0
                groupmutinfoPi(t)=groupmutinfoPi(t)+pxy*log(pxy/px/py);
            end
        end
    end
end

groupmutinfoC=zeros(1,Ntime);

for t=1:Ntime
    for i=1:Nb1
        px=sum(infomatC(i,:,t))/(N*its);
        for j=1:Nb2
            py=sum(infomatC(:,j,t))/(N*its);
            pxy=infomatC(i,j,t)/(N*its);
            if pxy~=0
                groupmutinfoC(t)=groupmutinfoC(t)+pxy*log(pxy/px/py);
            end
        end
    end
end

% meanskewnessD=mean(skewvalsD);
% meanskewnessR=mean(skewvalsR);
% meanskewnessDelta=mean(skewvalsDelta);
% meanskewnessH=mean(skewvalsH);
% meanskewnessPi=mean(skewvalsPi);
% meanskewnessC=mean(skewvalsC);

% powervecarray=[powervecD,powervecR,powervecDelta,powervecH,powervecPi,powervecC];
groupmutarray=[groupmutinfoD;groupmutinfoR;groupmutinfoDelta;groupmutinfoH;groupmutinfoPi;groupmutinfoC];
% meanskewnessarray=[meanskewnessD,meanskewnessR,meanskewnessDelta,meanskewnessH,meanskewnessPi,meanskewnessC];
% infomatlist={infomatD,infomatR,infomatDelta,infomatH,infomatPi,infomatC};

toreturn=groupmutarray;
end