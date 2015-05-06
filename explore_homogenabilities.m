numworkers=str2num(getenv('PROCS'));
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

load variables.mat
load twomat.mat

N=20;
its=1000;
cstep=.1;
vec=0:cstep:1;
Nc=length(vec);
endvals=zeros(1,Nc);
leak=2;

% X=repmat(0:cstep:1,1,11);
% Y=reshape(repmat(0:cstep:1,11,1),1,[]);

X=[];
Y=[];
layers=cell(Nc,1);
for i=0:(Nc-1)
    layers{i+1}=length(X)+(1:(Nc-i));
    endvals(i+1)=layers{i+1}(end);
    X=[X, vec(1:(Nc-i))]; %#ok<*AGROW>
    Y=[Y,vec(i+1)*ones(1,Nc-i)];   
end
Z=1-X-Y;
Nc=length(X);

avgthreshmat=zeros(N,Nc);
finalthreshmat=zeros(N,its,Nc);
threshvalscell=cell(1,Nc);

parfor i=1:Nc
    c1=X(i);
    c2=Y(i);
    c3=Z(i);
% toreturn={avgthresh,finalthresh,threshvalscell,probmat,timemat,meanskewness};
    t=group_props_homogen(c1,c2,c3,its,N,'homogen',Nd,Nt,twomat,domvals,threshvals,leak);
    avgthreshmat(:,i)=t{1};
    finalthreshmat(:,:,i)=t{2};
    threshvalscell{i}=t{3};
end

save('/home/brush/signaling_network/homogenabilities.mat','avgthreshmat','finalthreshmat','threshvalscell')

delete(dellapool);

exit ;
