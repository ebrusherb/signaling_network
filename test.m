numworkers=str2num(getenv('PROCS'));
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

load variables.mat
load twomat2.mat

v=[Nd,Nl,Nt];

save('/home/brush/signaling_network/test.mat','v')

delete(dellapool);

exit ;
