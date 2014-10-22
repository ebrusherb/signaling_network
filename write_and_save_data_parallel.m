numworkers=str2num(getenv('PROCS'));
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

lvals=0:.1:1;
Nl=length(lvals);

domvals=.5:.05:1;
Nd=length(domvals);

threshvals=.5:.1:2;
Nt=length(threshvals);

b=1;
deltax=0.025;

% twomat=zeros(Nl,Nd,Nt,Nt,2); %leak, dominance, thresholds, prob, time
twomat_par=zeros((Nd-3)*Nt*Nt,2);

k=2;
l=lvals(k);

parfor loopvar=1:(Nd-3)*Nt*Nt
% parfor loopvar=1:3
% for loopvar=1
    [u,i,j]=ind2sub([(Nd-3),Nt,Nt],loopvar);
    u=u+3;
    
    d=domvals(u);
    
    T1=threshvals(i);

    T2=threshvals(j);
    [x,y,~]=solve_pde_2d(-T1,-T2,d,5,5,deltax,deltax,l,b);
%     storevar=sub2ind([Nd,Nt,Nt],u,i,j);
    twomat_par(loopvar,:)=[x y];
       
end

twomat_par=reshape(twomat_par,Nd-3,Nt,Nt,2);

save('/home/brush/signaling_network/twomat_par.mat','twomat_par')

delete(dellapool);

exit ;