numworkers=str2num(getenv('PROCS'));
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

lvals=0:.1:1;
Nl=length(lvals);

domvals=.5:.05:1;
Nd=length(domvals);

threshvals=[.1 .5:.1:2 3:5];
Nt=length(threshvals);

b=1;
deltax=0.025;

x0=.1;
y0=1;

% twomat=zeros(Nl,Nd,Nt,Nt,2); %leak, dominance, thresholds, prob, time
twomat_par=zeros(Nd*Nt*Nt,2);
onemat_par=zeros(Nd*Nt*Nt,2);

k=2;
l=lvals(k);

parfor loopvar=1:Nd*Nt*Nt
% parfor loopvar=1:3
% for loopvar=1
    [u,i,j]=ind2sub([Nd,Nt,Nt],loopvar);
    
    d=domvals(u);
    
    T1=threshvals(i);

    T2=threshvals(j);
    [x,y,~]=solve_pde_2d_init(-T1,-T2,d,5,5,deltax,deltax,l,b,x0,y0);
%     storevar=sub2ind([Nd,Nt,Nt],u,i,j);
    twomat_par(loopvar,:)=[x y];
    
    [x,y,~]=solve_pde_1d_init(-T1,d,T2,deltax,l,b,x0-y0);
    onemat_par(loopvar,:)=[x y];
end

twomat_par=reshape(twomat_par,Nd,Nt,Nt,2);
onemat_par=reshape(onemat_par,Nd,Nt,Nt,2);

save('/home/brush/signaling_network/solve_pde_init.mat','twomat_par','onemat_par','threshvals','Nt','x0','y0')

delete(dellapool);

exit ;