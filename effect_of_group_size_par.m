numworkers=str2num(getenv('PROCS'));
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

load variables.mat
load twomat.mat

sizes=[100 500 1000];
Ns=length(sizes);

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

accuracy=zeros(Ns,Nlays*2);
bottomaccuracy=zeros(Ns,Nlays*2);
topaccuracy=zeros(Ns,Nlays*2);

% N=20;
% c1=0;c2=0;c3=1;
% t=group_props_parallelized(c1,c2,c3,its,N,'unif',Nd,Nt,twomat,domvals,threshvals,leak);

for s=1:Ns
    N=sizes(s);
    q=N/4;
    topquart=1:q;
    bottomquart=fliplr(N+1-(1:q));
    parfor ind=1:Nlays*2
        [l,k]=ind2sub([Nlays 2],ind);
%     for l=1:Nlays
%         for k=1:2
        if k==1
            i=layers{l}(1);
        else i=layers{l}(end);
        end
        c1=X(i);
        c2=Y(i);
        c3=Z(i);
        t=group_props_parallelized(c1,c2,c3,its,N,'unif',Nd,Nt,twomat,domvals,threshvals,leak);
        accuracy(s,ind)=t{4};
        totalsigmat=t{11};
        bottomaccuracy(s,ind)=sum(sum(triu(totalsigmat(bottomquart,bottomquart),1)))/sum(sum(totalsigmat(bottomquart,bottomquart)));
        topaccuracy(s,ind)=sum(sum(triu(totalsigmat(topquart,topquart),1)))/sum(sum(totalsigmat(topquart,topquart)));
%         end
%     end
    end
    save('/home/brush/signaling_network/groupsize.mat','accuracy','topaccuracy','bottomaccuracy','layers','Nlays','sizes','Ns')
end

accuracy=reshape(accuracy,Ns,Nlays,2);
bottomaccuracy=reshape(bottomaccuracy,Ns,Nlays,2);
topaccuracy=reshape(topaccuracy,Ns,Nlays,2);

save('/home/brush/signaling_network/groupsize.mat','accuracy','topaccuracy','bottomaccuracy','layers','Nlays','sizes','Ns')

delete(dellapool);

exit ;