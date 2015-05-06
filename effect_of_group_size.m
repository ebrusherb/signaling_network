numworkers=str2num(getenv('PROCS'));
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

load variables.mat
load twomat.mat

sizes=[500];
Ns=length(sizes);

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

Nlays=Nc;
Nc=length(X);

leak=2;

accuracy=zeros(Ns,Nlays,2);
bottomaccuracy=zeros(Ns,Nlays,2);
topaccuracy=zeros(Ns,Nlays,2);

for s=1:Ns
    N=sizes(s);
    q=N/4;
    topquart=1:q;
    bottomquart=fliplr(N+1-(1:q));
    for l=1:Nlays
        for k=1:2
        if k==1
            i=layers{l}(1);
        else i=layers{l}(end);
        end
        c1=X(i);
        c2=Y(i);
        c3=Z(i);
        t=group_props_parallelized(c1,c2,c3,its,N,'unif',Nd,Nt,twomat,domvals,threshvals,leak);
        accuracy(s,l,k)=t{4};
        totalsigmat=t{11};
        bottomaccuracy(s,l,k)=sum(sum(triu(totalsigmat(bottomquart,bottomquart),1)))/sum(sum(totalsigmat(bottomquart,bottomquart)));
        topaccuracy(s,l,k)=sum(sum(triu(totalsigmat(topquart,topquart),1)))/sum(sum(totalsigmat(topquart,topquart)));
        end
    end
    save('/home/brush/signaling_network/groupsize.mat','accuracy','topaccuracy','bottomaccuracy','layers','Nlays','sizes','Ns')
end

delete(dellapool);

exit ;