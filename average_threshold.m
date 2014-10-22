function threshvec=average_threshold(c1,c2,c3,its,N,dist)

if nargin==3
    N=30;
    its=25;
end

Nl=evalin('base','Nl');
Nd=evalin('base','Nd');
Nt=evalin('base','Nt');
twomat=evalin('base','twomat');
domvals=evalin('base','domvals');
threshvals=evalin('base','threshvals');

Nb1=10;
Nb2=10;
infomat=zeros(Nb1,Nb2);
skewvals=zeros(1,its);
finalthresh=zeros(N,its);
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
end
    threshvec=mean(finalthresh,2);
end
