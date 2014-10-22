% c11=.4;c21ç=.6;c31=0;
% c12=.9;c22=.1;c32=0;
% c11=0;c21=.4;c31=.6;
% c12=.9;c22=.1;c32=0;
% c11=0;c21=.55;c31=.45;
% c12=.85;c22=.2;c32=.05;

% c11=0;c21=.3;c31=.7;
% c12=.7;c22=.3;c32=0;
% c13=0;c23=.3;c33=.7;
time=0.05;
c11=0;c21=time;c31=1-time;
c12=1-time;c22=time;c32=0;

cvals=[c11 c21 c31; c12 c22 c32];
% cvals=[c12,c22,c32];
% cvals=[c11 c21 c31];

[seqoptthresh, dmatssaved, pairwiseaccuracy] = paradox_of_mutual_info(cvals,domvals,twomat);

%%
c=1;
ic=1;
otheric=2;
[thresh2perf(cvals(ic,:),seqoptthresh(c,ic,:),twomat,dmatssaved(c,:,:),l), thresh2perf(cvals(ic,:),seqoptthresh(c,otheric,:),twomat,dmatssaved(c,:,:),l)]

%%
c=1;
ic=1;
otheric=2;
dmat=reshape(dmatssaved(c,:,:),N,N);
        c1=cvals(ic,1);
        c2=cvals(ic,2);
        c3=cvals(ic,3);
        

        perf=zeros(2,Nd,Nt,Nt); %individual, leak, dominance, thresholds

        perf(1,2:Nd,:,:)=c1*(1-twomat(leak,2:Nd,:,:,1))+c2*(twomat(leak,2:Nd,:,:,2))+c3*(1-twomat(leak,2:Nd,:,:,1));
        perf(2,2:Nd,:,:)=c1*(1-twomat(leak,2:Nd,:,:,1))+c2*(twomat(leak,2:Nd,:,:,2))+c3*(twomat(leak,2:Nd,:,:,1));

        perf(1,1,:,:)=c2*(twomat(leak,1,:,:,2))+c3*(1-twomat(leak,1,:,:,1));
        perf(2,1,:,:)=c2*(twomat(leak,1,:,:,2))+c3*(twomat(leak,1,:,:,1));

        Tvals=zeros(N,opt_its);
%         Tvals(:,1)=2*ones(N,1);
        % Tvals(:,1)=randi([1 Nt],N,1);
%         Tvals(:,1)=reshape(seqoptthresh(c,ic,:),[],1);
%         Tvals(16,1)=seqoptthresh(c,otheric,3);
        Tvals(:,1)=col(seqoptthresh(c,otheric,:));
        perfvals=zeros(N,opt_its);

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
                            perfsum=perfsum+perf(2,Nd-ds(q)+1,opp_thresh(q),j);
                        else 
                            perfsum=perfsum+perf(1,ds(q)-Nd+1,j,opp_thresh(q));
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
        
%%
c=1;
tvec1=seqoptthresh(c,1,:);
tvec2=seqoptthresh(c,2,:);
cvecs=[cvals(1,:); cvals(2,:); [1 0 0]; [0 1 0]; [0 0 1]];
d=size(cvecs,1);
diffmat=zeros(N,d);
for i=1:d
    cvec=cvecs(i,:);
    diffmat(:,i)=thresh2perf(cvec,tvec1,twomat,dmatssaved(c,:,:),l)-thresh2perf(cvec,tvec2,twomat,dmatssaved(c,:,:),l);
end
 
%%
c=1;
tvec1=seqoptthresh(c,1,:);
tvec2=tvec1;
k=3;
tvec2(k)=seqoptthresh(c,2,k);
cvecs=[cvals(1,:); cvals(2,:); [1 0 0]; [0 1 0]; [0 0 1]];
d=size(cvecs,1);
diffmat=zeros(N,d);
for i=1:d
    cvec=cvecs(i,:);
    diffmat(:,i)=thresh2perf(cvec,tvec1,twomat,dmatssaved(c,:,:),leak)-thresh2perf(cvec,tvec2,twomat,dmatssaved(c,:,:),leak);
end


%%
c=1;
ic=2;
otheric=1;
tvec1=seqoptthresh(c,ic,:);

k=12;

cvecs=[cvals(ic,:)];
d=size(cvecs,1);
diffmat=zeros(N,N);
for i=1:N
    cvec=cvecs(1,:);
    tvec2=tvec1;
    tvec2(i)=seqoptthresh(c,otheric,i);
    diffmat(:,i)=thresh2perf(cvec,tvec1,twomat,dmatssaved(c,:,:),leak)-thresh2perf(cvec,tvec2,twomat,dmatssaved(c,:,:),leak);
end

%%
c=1;
ic=1;
otheric=2;
dmat=reshape(dmatssaved(c,:,:),N,N);
c1=cvals(ic,1);
c2=cvals(ic,2);
c3=cvals(ic,3);
tvec=reshape(seqoptthresh(c,ic,:),[],1);
k=6;        
perfmat=zeros(N,Nt,3);
id=eye(3);

for j=1:Nt
    for c=1:3
    tvec(k)=j;
    perfmat(:,j,c)=reshape(thresh2perf(id(c,:),tvec,twomat,dmatssaved(c,:,:),leak),[],1);
    end
end

%%
t1=group_props_parallelized(0,0,1,its,N,'unif',Nd,Nt,twomat,domvals,threshvals);
t2=group_props_parallelized(1,0,0,its,N,'unif',Nd,Nt,twomat,domvals,threshvals);

    