its=10;
N=20;
dist='unif';

c11=0;c21=.55;c31=.45;
c12=.85;c22=.2;c32=.05;

Nb1=5;
Nb2=5;
infomat1=zeros(Nb1,Nb2);
infomat2=zeros(Nb1,Nb2);
skewvals1=zeros(1,its);
skewvals2=zeros(1,its);
finalthresh1=zeros(N,its);
finalthresh2=zeros(N,its);

indinfomats1=zeros(N,2,2);
indinfomats2=zeros(N,2,2);
indmutinfos1=zeros(N,its);
indmutinfos2=zeros(N,its);

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
    
    perf1=zeros(2,Nl,Nd,Nt,Nt); %individual, leak, dominance, thresholds

    perf1(1,:,2:Nd,:,:)=c11*(1-twomat(:,2:Nd,:,:,1))+c21*(twomat(:,2:Nd,:,:,2))+c31*(1-twomat(:,2:Nd,:,:,1));
    perf1(2,:,2:Nd,:,:)=c11*(1-twomat(:,2:Nd,:,:,1))+c21*(twomat(:,2:Nd,:,:,2))+c31*(twomat(:,2:Nd,:,:,1));

    perf1(1,:,1,:,:)=c21*(twomat(:,1,:,:,2))+c31*(1-twomat(:,1,:,:,1));
    perf1(2,:,1,:,:)=c21*(twomat(:,1,:,:,2))+c31*(twomat(:,1,:,:,1));
    
    perf2=zeros(2,Nl,Nd,Nt,Nt); %individual, leak, dominance, thresholds

    perf2(1,:,2:Nd,:,:)=c12*(1-twomat(:,2:Nd,:,:,1))+c22*(twomat(:,2:Nd,:,:,2))+c32*(1-twomat(:,2:Nd,:,:,1));
    perf2(2,:,2:Nd,:,:)=c12*(1-twomat(:,2:Nd,:,:,1))+c22*(twomat(:,2:Nd,:,:,2))+c32*(twomat(:,2:Nd,:,:,1));

    perf2(1,:,1,:,:)=c22*(twomat(:,1,:,:,2))+c32*(1-twomat(:,1,:,:,1));
    perf2(2,:,1,:,:)=c22*(twomat(:,1,:,:,2))+c32*(twomat(:,1,:,:,1));

    l=1;

    Tvals1=zeros(N,opt_its);
    Tvals1(:,1)=2*ones(N,1);
    % Tvals(:,1)=randi([1 Nt],N,1);
    perfvals1=zeros(N,opt_its);
    
    Tvals2=zeros(N,opt_its);
    Tvals2(:,1)=2*ones(N,1);
    perfvals2=zeros(N,opt_its);

    for i=1:N
        opp_thresh=Tvals1([1:(i-1),(i+1):N],1);
        ds=dmat(i,[1:(i-1),(i+1):N]);
        perfsum=0;
                for q=1:(N-1)
                    if ds(q)<=Nd-1
                        perfsum=perfsum+perf1(2,l,Nd-ds(q)+1,opp_thresh(q),Tvals1(i,1)); 
                    else 
                        perfsum=perfsum+perf1(1,l,ds(q)-Nd+1,Tvals1(i,1),opp_thresh(q));
                    end
                end  
        perfvals1(i,1)=perfsum;
        
        opp_thresh=Tvals2([1:(i-1),(i+1):N],1);
        ds=dmat(i,[1:(i-1),(i+1):N]);
        perfsum=0;
                for q=1:(N-1)
                    if ds(q)<=Nd-1
                        perfsum=perfsum+perf2(2,l,Nd-ds(q)+1,opp_thresh(q),Tvals2(i,1)); 
                    else 
                        perfsum=perfsum+perf2(1,l,ds(q)-Nd+1,Tvals2(i,1),opp_thresh(q));
                    end
                end  
        perfvals2(i,1)=perfsum;
    end

    count=1;
    while count<=opt_its
        for i=1:N
            opp_thresh=Tvals1([1:(i-1),(i+1):N],count);
            ds=dmat(i,[1:(i-1),(i+1):N]);
            perftest=zeros(1,Nt);
            for j=1:Nt
                perfsum=0;
                for q=1:(N-1)
                    if ds(q)<=Nd-1
                        perfsum=perfsum+perf1(2,l,Nd-ds(q)+1,opp_thresh(q),j);
                    else 
                        perfsum=perfsum+perf1(1,l,ds(q)-Nd+1,j,opp_thresh(q));
                    end
                end 
                perftest(j)=perfsum;
            end
            [p,n]=min(perftest);
            Tvals1(i,count+1)=n;
            perfvals1(i,count+1)=p;
            
            opp_thresh=Tvals2([1:(i-1),(i+1):N],count);
            ds=dmat(i,[1:(i-1),(i+1):N]);
            perftest=zeros(1,Nt);
            for j=1:Nt
                perfsum=0;
                for q=1:(N-1)
                    if ds(q)<=Nd-1
                        perfsum=perfsum+perf2(2,l,Nd-ds(q)+1,opp_thresh(q),j);
                    else 
                        perfsum=perfsum+perf2(1,l,ds(q)-Nd+1,j,opp_thresh(q));
                    end
                end 
                perftest(j)=perfsum;
            end
            [p,n]=min(perftest);
            Tvals2(i,count+1)=n;
            perfvals2(i,count+1)=p;
        end
        if sum(Tvals1(:,count+1)==Tvals1(:,count))==N && sum(Tvals2(:,count+1)==Tvals2(:,count))==N
            maxit=count+2;
            count=opt_its+1;
            Tvals1(:,maxit:end)=[];
            perfvals1(:,maxit:end)=[];
            
            Tvals2(:,maxit:end)=[];
            perfvals2(:,maxit:end)=[];
        end
        count=count+1;
    end
    
    finalthresh1(:,c)=threshvals(Tvals1(:,end));
    if c==its
        threshtohist1=threshvals(Tvals1(:,end));
    end
    finalthresh2(:,c)=threshvals(Tvals2(:,end));
    if c==its
        threshtohist2=threshvals(Tvals2(:,end));
    end

    probmat1=zeros(N,N);
    timemat1=zeros(N,N);
    
    probmat2=zeros(N,N);
    timemat2=zeros(N,N);

    %size(twomat)=[Nl,Nd,Nt,Nt,2];

    for i=1:N
        for j=[1:(i-1),(i+1):N]
            if dmat(i,j)<=Nd-1
                probmat1(i,j)=1-twomat(l,Nd-dmat(i,j)+1,Tvals1(j,end),Tvals1(i,end),1);
                timemat1(i,j)=twomat(l,Nd-dmat(i,j)+1,Tvals1(j,end),Tvals1(i,end),2);
            else
                probmat1(i,j)=twomat(l,dmat(i,j)-Nd+1,Tvals1(i,end),Tvals1(j,end),1);
                timemat1(i,j)=twomat(l,dmat(i,j)-Nd+1,Tvals1(i,end),Tvals1(j,end),2);
            end
        end
    end
    
    for i=1:N
        for j=[1:(i-1),(i+1):N]
            if dmat(i,j)<=Nd-1
                probmat2(i,j)=1-twomat(l,Nd-dmat(i,j)+1,Tvals2(j,end),Tvals2(i,end),1);
                timemat2(i,j)=twomat(l,Nd-dmat(i,j)+1,Tvals2(j,end),Tvals2(i,end),2);
            else
                probmat2(i,j)=twomat(l,dmat(i,j)-Nd+1,Tvals2(i,end),Tvals2(j,end),1);
                timemat2(i,j)=twomat(l,dmat(i,j)-Nd+1,Tvals2(i,end),Tvals2(j,end),2);
            end
        end
    end
    
    numsig1=sum(probmat1,2);
    numsig2=sum(probmat2,2);
    
    abilitiesdiff=ceil(max(fighting_abilities))-floor(min(fighting_abilities));
    abilitiesbins=floor(min(fighting_abilities)):abilitiesdiff/Nb1:ceil(max(fighting_abilities));
    [~,abilitiesindices]=histc(fighting_abilities,abilitiesbins);
    
    sigdiff1=ceil(max(numsig1))-floor(min(numsig1));
    sigbins1=floor(min(numsig1)):sigdiff1/Nb2:ceil(max(numsig1));
    [~,sigindices1]=histc(numsig1,sigbins1);
    
    sigdiff2=ceil(max(numsig2))-floor(min(numsig2));
    sigbins2=floor(min(numsig2)):sigdiff2/Nb2:ceil(max(numsig2));
    [~,sigindices2]=histc(numsig2,sigbins2);
    
    for i=1:N
        infomat1(abilitiesindices(i),sigindices1(i))=infomat1(abilitiesindices(i),sigindices1(i))+1;
        infomat2(abilitiesindices(i),sigindices2(i))=infomat2(abilitiesindices(i),sigindices2(i))+1;
    end
    
    skewvals1(c)=skewness(numsig1);
    skewvals2(c)=skewness(numsig2);
    
    for i=1:N
        for j=1:(i-1) 
            indinfomats1(i,1,1)=indinfomats1(i,1,1)+probmat1(j,i);
            indinfomats1(i,1,2)=indinfomats1(i,1,2)+probmat1(i,j);
            
            indinfomats2(i,1,1)=indinfomats2(i,1,1)+probmat2(j,i);
            indinfomats2(i,1,2)=indinfomats2(i,1,2)+probmat2(i,j);
        end
        for j=(i+1):N
            indinfomats1(i,2,1)=indinfomats1(i,2,1)+probmat1(j,i);
            indinfomats1(i,2,2)=indinfomats1(i,2,2)+probmat1(i,j);
            
            indinfomats2(i,2,1)=indinfomats2(i,2,1)+probmat1(j,i);
            indinfomats2(i,2,2)=indinfomats2(i,2,2)+probmat1(i,j);
        end
    end
    indinfomats1=indinfomats1/(N-1);
    indinfomats2=indinfomats2/(N-1);
%     indinfomats(25,:)

    for i=1:N
        for j=1:2
            px=sum(indinfomats1(i,j,:));
            for k=1:2
               py=sum(indinfomats1(i,:,k));
               pxy=indinfomats1(i,j,k);
               if pxy~=0
                   indmutinfos1(i,c)=indmutinfos1(i,c)+pxy*log(pxy/px/py);
               end
            end
        end
    end
    
    for i=1:N
        for j=1:2
            px=sum(indinfomats2(i,j,:));
            for k=1:2
               py=sum(indinfomats2(i,:,k));
               pxy=indinfomats2(i,j,k);
               if pxy~=0
                   indmutinfos2(i,c)=indmutinfos2(i,c)+pxy*log(pxy/px/py);
               end
            end
        end
    end
end

indmutinfos1=mean(indmutinfos1,2);
indmutinfos2=mean(indmutinfos2,2);

groupmutinfo1=0;
groupmutinfo2=0;

mutmat1=zeros(Nb1,Nb2);
mutmat2=zeros(Nb1,Nb2);

for i=1:Nb1
    px=sum(infomat1(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomat1(:,j))/(N*its);
        pxy=infomat1(i,j)/(N*its);
        if pxy~=0
            groupmutinfo1=groupmutinfo1+pxy*log(pxy/px/py);
            mutmat1(i,j)=pxy*log(pxy/px/py);
        end
    end
end

for i=1:Nb1
    px=sum(infomat2(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomat2(:,j))/(N*its);
        pxy=infomat2(i,j)/(N*its);
        if pxy~=0
            groupmutinfo2=groupmutinfo2+pxy*log(pxy/px/py);
            mutmat2(i,j)=pxy*log(pxy/px/py);
        end
    end
end

avgthresh1=mean(finalthresh1,2);
avgthresh2=mean(finalthresh2,2);

figure
subplot(2,4,1)
plot(avgthresh1)
subplot(2,4,2)
imagesc(infomat1)
subplot(2,4,3)
imagesc(mutmat1)
colorbar
subplot(2,4,4)
imagesc(probmat1)
colorbar

subplot(2,4,5)
plot(avgthresh2)
subplot(2,4,6)
imagesc(infomat2)
subplot(2,4,7)
imagesc(mutmat2)
colorbar
subplot(2,4,8)
imagesc(probmat2)
colorbar

groupmutinfo1
groupmutinfo2



