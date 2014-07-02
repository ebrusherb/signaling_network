its=1;
N=20;
dist='unif';
extendeddomvals=[1-domvals(end:-1:2) domvals];
v=1:(Nd+Nd-1);
deps=.05;

% c11=.4;c21ç=.6;c31=0;
% c12=.9;c22=.1;c32=0;
% c11=0;c21=.4;c31=.6;
% c12=.9;c22=.1;c32=0;
% c11=0;c21=.55;c31=.45;
% c12=.85;c22=.2;c32=.05;

c11=0;c21=.3;c31=7;
c12=.7;c22=.3;c32=0;
c13=0;c23=.3;c33=.7;

% cvals=[c11 c21 c31; c12 c22 c32];
cvals=[c12,c22,c32];
totalc=size(cvals,1);

Nb1=10;
Nb2=10;

seqoptthresh=zeros(its,totalc,N);
seqprobmat=zeros(its,totalc,N,N);
seqaccmat=zeros(its,totalc,N,N);
seqtimemat=zeros(its,totalc,N,N);
seqsigmat=zeros(its,totalc,N,N);
perfvalssaved=zeros(its,totalc,N);

for c=1:its

    for ic=1:totalc
 
        c1=cvals(ic,1);
        c2=cvals(ic,2);
        c3=cvals(ic,3);
        l=1;
    
        opt_its=10;

        fighting_abilities=zeros(1,N);
        
        fighting_abilities=20*rand(1,N);

        fighting_abilities=sort(fighting_abilities,'descend');

        dmat=zeros(N); %turn difference in fighting abilities into probabilitise of winning
        
        for i=1:N
            for j=[1:(i-1),(i+1):N]
                diff=fighting_abilities(i)-fighting_abilities(j);
                d=exp(diff)/(exp(diff)+1);
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
                [p,n]=min(perftest);
                Tvals(i,count+1)=n;
                
            end
            
            if sum(Tvals(:,count+1)==Tvals(:,count))==N
                Tvals(:,count+2)=Tvals(:,count+1);
                perfvec=zeros(N,1);
                for i=1:N
                    opp_thresh=Tvals([1:(i-1),(i+1):N],count+1);
                    ds=dmat(i,[1:(i-1),(i+1):N]);
                    perfsum=0;
                            for q=1:(N-1)
                                if ds(q)<=Nd-1
                                    perfsum=perfsum+perf(2,Nd-ds(q)+1,opp_thresh(q),Tvals(i,count+1)); 
                                else 
                                    perfsum=perfsum+perf(1,ds(q)-Nd+1,Tvals(i,count+1),opp_thresh(q));
                                end
                            end  
                    perfvec(i)=perfsum;
                end
                
                maxit=count+3;
                count=opt_its+1;
                Tvals(:,maxit:end)=[];
                seqoptthresh(c,ic,:)=Tvals(:,end);
            end
            count=count+1;
        end
        
       perfvalssaved(c,ic,:)=perfvec;
      
    end
  
end

% avgoptthresh=reshape(mean(seqoptthresh,1),2,[])';

%%
[probmat1,timemat1,accmat1]=thresh2mat(seqoptthresh(end,1,:),twomat,dmat,l);
perfrightthresh1=cvals(1,1)*(1-accmat1)+cvals(1,2)*timemat1+cvals(1,3)*(1-probmat1);
[reshape(perfvalssaved(end,1,:),[],1) sum(perfrightthresh1,2)]

%%
[probmat2,timemat2,accmat2]=thresh2mat(seqoptthresh(end,2,:),twomat,dmat,l);
perfrightthresh2=cvals(2,1)*(1-accmat2)+cvals(2,2)*timemat2+cvals(2,3)*(1-probmat2);
[reshape(perfvalssaved(end,2,:),[],1)-sum(perfrightthresh2,2)]
%%
[probmat3,timemat3,accmat3]=thresh2mat(seqoptthresh(end,3,:),twomat,dmat,l);
perfrightthresh3=cvals(3,1)*(1-accmat3)+cvals(3,2)*timemat3+cvals(3,3)*(1-probmat3);
[reshape(perfvalssaved(end,2,:),[],1) sum(perfrightthresh3,2)]
%%

[probmat1,timemat1,accmat1]=thresh2mat(seqoptthresh(end,1,:),twomat,dmat,l);
[probmat2,timemat2,accmat2]=thresh2mat(seqoptthresh(end,2,:),twomat,dmat,l);

perfrightthresh1=cvals(1,1)*(1-accmat1)+cvals(1,2)*timemat1+cvals(1,3)*(1-probmat1);
perfwrongthresh1=cvals(1,1)*(1-accmat2)+cvals(1,2)*timemat2+cvals(1,3)*(1-probmat2);

perfrightthresh2=cvals(2,1)*(1-accmat2)+cvals(2,2)*timemat2+cvals(2,3)*(1-probmat2);
perfwrongthresh2=cvals(2,1)*(1-accmat1)+cvals(2,2)*timemat1+cvals(2,3)*(1-probmat1);


%%
accuracylookup=zeros(2,Nd,Nt,Nt); %individual, leak, dominance, thresholds

accuracylookup(1,2:Nd,:,:)=(1-twomat(l,2:Nd,:,:,1));
accuracylookup(2,2:Nd,:,:)=(1-twomat(l,2:Nd,:,:,1));

accuracylookup(1,1,:,:)=0;
accuracylookup(2,1,:,:)=0;

dominancelookup=zeros(2,Nd,Nt,Nt); %individual, leak, dominance, thresholds

dominancelookup(1,1:Nd,:,:)=(1-twomat(l,1:Nd,:,:,1));
dominancelookup(2,1:Nd,:,:)=(twomat(l,1:Nd,:,:,1));

Tvals=reshape(seqoptthresh(its,:,:),totalc,[]); 
accperfvals=zeros(totalc,N);
domperfvals=zeros(totalc,N);

for ic=1:totalc
    
for i=1:N
    opp_thresh=Tvals(ic,[1:(i-1),(i+1):N]);
    ds=dmat(i,[1:(i-1),(i+1):N]);
    perfsum=0;
        for q=1:(N-1)
            if ds(q)<=Nd-1
                perfsum=perfsum+accuracylookup(2,Nd-ds(q)+1,opp_thresh(q),Tvals(ic,i)); 
            else 
                perfsum=perfsum+accuracylookup(1,ds(q)-Nd+1,Tvals(ic,i),opp_thresh(q));
            end
        end  
    accperfvals(ic,i)=perfsum;
    perfsum=0;
        for q=1:(N-1)
            if ds(q)<=Nd-1
                perfsum=perfsum+dominancelookup(2,Nd-ds(q)+1,opp_thresh(q),Tvals(ic,i)); 
            else 
                perfsum=perfsum+dominancelookup(1,ds(q)-Nd+1,Tvals(ic,i),opp_thresh(q));
            end
        end 
    domperfvals(ic,i)=perfsum;
end
end
