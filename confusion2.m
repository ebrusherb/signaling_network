its=2;
N=20;
dist='unif';
extendeddomvals=[1-domvals(end:-1:2) domvals];
v=1:(Nd+Nd-1);
deps=.05;
l=1;
opt_its=50;

% c11=.4;c21ç=.6;c31=0;
% c12=.9;c22=.1;c32=0;
% c11=0;c21=.4;c31=.6;
% c12=.9;c22=.1;c32=0;
% c11=0;c21=.55;c31=.45;
% c12=.85;c22=.2;c32=.05;

c11=0;c21=.3;c31=.7;
c12=.7;c22=.3;c32=0;
c13=0;c23=.3;c33=.7;

cvals=[c11 c21 c31; c12 c22 c32];
% cvals=[c12,c22,c32];
% cvals=[c11 c21 c31];
totalc=size(cvals,1);

seqoptthresh=zeros(its,totalc,N);
dmatssaved=zeros(its,N,N);
perfvalssaved=zeros(its,totalc,N);

for c=1:its

    fighting_abilities=20*rand(1,N);
    fighting_abilities=sort(fighting_abilities,'descend');

    dmat=zeros(N);
    
    for i=1:N
        for j=[1:(i-1),(i+1):N]
            diff=fighting_abilities(i)-fighting_abilities(j);
            d=exp(diff)/(exp(diff)+1);
            d=floor(round(d/deps))*deps;
            dmat(i,j)=v(abs(d-extendeddomvals)<=deps/2);
        end
    end

    dmatssaved(c,:,:)=dmat;
    
    for ic=1:totalc
 
        c1=cvals(ic,1);
        c2=cvals(ic,2);
        c3=cvals(ic,3);
        
        Tvals=zeros(N,opt_its);
        Tvals(:,1)=2*ones(N,1);

        count=1;
        while count<=opt_its
            for i=1:N
                testthreshvals=Tvals(:,count);
                perftest=zeros(Nt,1);
                for j=1:Nt
                    testthreshvals(i)=j;
                    testperfvec=thresh2perf([c1 c2 c3],testthreshvals,twomat,dmat,l);
                    perftest(j)=testperfvec(i);
                end
                
                [p,n]=min(perftest);
                Tvals(i,count+1)=n;  
                
            end
            
            if sum(Tvals(:,count+1)==Tvals(:,count))==N
                Tvals(:,count+2)=Tvals(:,count+1);
                maxit=count+3;
                count=opt_its+1;
                Tvals(:,maxit:end)=[];    
            else 
                count=count+1;
            end
        end
        seqoptthresh(c,ic,:)=Tvals(:,end);
        perfvalssaved(c,ic,:)=thresh2perf(cvals(ic,:),Tvals(:,end),twomat,dmat,l);
        
    end
  
end