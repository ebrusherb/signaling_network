function toreturn=group_props_parallelized(c1,c2,c3,its,N,dist,Nl,Nd,Nt,twomat,domvals,threshvals)

if nargin==3
    N=30;
    its=25;
end

Nb1=10;
Nb2=10;
infomat=zeros(Nb1,Nb2);
skewvals=zeros(1,its);
finalthresh=zeros(N,its);

indinfomats=zeros(N,2,2);
indmutinfos=zeros(N,its);

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
            diff=fighting_abilities(i)-fighting_abilities(j);
            d=exp(diff)/(exp(diff)+1);
            d=floor(round(d/deps))*deps;
            dmat(i,j)=v(abs(d-extendeddomvals)<=deps/2);
        end
    end

    %size(twomat)=[Nl,Nd,Nt,Nt,2];
    
    perf=zeros(2,Nl,Nd,Nt,Nt); %individual, leak, dominance, thresholds

    perf(1,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2))+c3*(1-twomat(:,2:Nd,:,:,1));
    perf(2,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2))+c3*(twomat(:,2:Nd,:,:,1));

    perf(1,:,1,:,:)=c2*(twomat(:,1,:,:,2))+c3*(1-twomat(:,1,:,:,1));
    perf(2,:,1,:,:)=c2*(twomat(:,1,:,:,2))+c3*(twomat(:,1,:,:,1));

    l=1;

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
                        perfsum=perfsum+perf(2,l,Nd-ds(q)+1,opp_thresh(q),Tvals(i,1)); 
                    else 
                        perfsum=perfsum+perf(1,l,ds(q)-Nd+1,Tvals(i,1),opp_thresh(q));
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
                        perfsum=perfsum+perf(2,l,Nd-ds(q)+1,opp_thresh(q),j);
                    else 
                        perfsum=perfsum+perf(1,l,ds(q)-Nd+1,j,opp_thresh(q));
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
            Tvals(:,maxit:end)=[];
            perfvals(:,maxit:end)=[];
        end
        count=count+1;
    end
    
    finalthresh(:,c)=threshvals(Tvals(:,end));
    if c==its
        threshtohist=threshvals(Tvals(:,end));
    end

    probmat=zeros(N,N);
    timemat=zeros(N,N);

    %size(twomat)=[Nl,Nd,Nt,Nt,2];

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
    
    sigmat=zeros(N,N);
    
    for i=1:N
        for j=(i+1):N
            draw=rand;
            if draw <+probmat(i,j)
                sigmat(i,j)=1;
            else sigmat(j,i)=1;
            end
        end
    end
    
%     power=sum(probmat,2);
    power=sum(sigmat,2);
    
    abilitiesdiff=ceil(max(fighting_abilities))-floor(min(fighting_abilities));
    abilitiesbins=floor(min(fighting_abilities)):abilitiesdiff/Nb1:ceil(max(fighting_abilities));
    [~,abilitiesindices]=histc(fighting_abilities,abilitiesbins);
    
    sigdiff=ceil(max(power))-floor(min(power))+.1;
    sigbins=(floor(min(power))-.05):(sigdiff/Nb2):(ceil(max(power))+.05);
    [~,sigindices]=histc(power,sigbins);
    
    for i=1:N
        infomat(abilitiesindices(i),sigindices(i))=infomat(abilitiesindices(i),sigindices(i))+1;
    end
    
    skewvals(c)=skewness(power);
    
    for i=1:N
        for j=1:(i-1) 
            indinfomats(i,1,1)=indinfomats(i,1,1)+probmat(j,i);
            indinfomats(i,1,2)=indinfomats(i,1,2)+probmat(i,j);
        end
        for j=(i+1):N
            indinfomats(i,2,1)=indinfomats(i,2,1)+probmat(j,i);
            indinfomats(i,2,2)=indinfomats(i,2,2)+probmat(i,j);
        end
    end
    indinfomats=indinfomats/(N-1);
%     indinfomats(25,:)

    for i=1:N
        for j=1:2
            px=sum(indinfomats(i,j,:));
            for k=1:2
               py=sum(indinfomats(i,:,k));
               pxy=indinfomats(i,j,k);
               if pxy~=0
                   indmutinfos(i,c)=indmutinfos(i,c)+pxy*log(pxy/px/py);
               end
            end
        end
    end
end

indmutinfos=mean(indmutinfos,2);

groupmutinfo=0;

for i=1:Nb1
    px=sum(infomat(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomat(:,j))/(N*its);
        pxy=infomat(i,j)/(N*its);
        if pxy~=0
            groupmutinfo=groupmutinfo+pxy*log(pxy/px/py);
        end
    end
end

meanskewness=mean(skewvals);

avgthresh=mean(finalthresh,2);

toreturn={avgthresh,threshtohist,probmat,timemat,indmutinfos,groupmutinfo,meanskewness};
end