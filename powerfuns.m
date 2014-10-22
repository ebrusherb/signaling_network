function [powervecD,powervecR,powervecDelta,powervecH,powervecPi,powervecC] = powerfuns(binsigmat,sigmat)
weight=.05;
N=size(sigmat,2);
todivide=sum(sigmat,2);
todivide(todivide==0)=1;
entmat=sigmat./repmat(todivide,1,N);
logmat=entmat;
logmat(entmat~=0)=-log(entmat(entmat~=0));
entmat=entmat.*logmat;

%         powervec=sum(probmat,2);
powervecD=sum(binsigmat,2);
powervecR=sum(sigmat,2);
powervecDelta=(powervecD).*(powervecR);
powervecH=sum(entmat,2);
powervecH=real(powervecH);
powervecPi=(powervecR).*(powervecH);
powervecC=eigcentrality(sigmat,weight);
% powervecprobs=sum(probmat,2);
% power=[powervecD,powervecR,powervecDelta,powervecH,powervecPi,powervecC];