
ivals=[1 11 12 39 52 57];
cvecs=zeros(length(ivals),3);
for i=1:length(ivals)
    cvecs(i,:)=[X(ivals(i)), Y(ivals(i)), Z(ivals(i))];
end

mycols={'red','blue','green'};
theta=0:pi/3:2*pi;
figure
hold on
runs=1000;
I=6;

powerindicesC=zeros(N,runs,I);
powerindicesR=zeros(N,runs,I);
powermatR=zeros(N,runs,I);
powermatC=zeros(N,runs,I);
powermatH=zeros(N,runs,I);
powermatD=zeros(N,runs,I);
keepsigmat=zeros(N,N,runs,I);

for ic=1:I
for r=1:runs
fighting_abilities=20*rand(1,N);
fighting_abilities=sort(fighting_abilities,'descend');

dmat=zeros(N); %turn difference in fighting abilities into probabilitise of winning
extendeddomvals=[1-domvals(end:-1:2) domvals];
v=1:(Nd+Nd-1);
deps=.05;

for i=1:N
    for j=[1:(i-1),(i+1):N]
        abilitydiff=fighting_abilities(i)-fighting_abilities(j);
        d=exp(abilitydiff)/(exp(abilitydiff)+1);
        d=floor(round(d/deps))*deps;
        dmat(i,j)=v(abs(d-extendeddomvals)<=deps/2);
    end
end

Tvals=abilities2thresh(cvecs(ic,:),twomat,dmat,opt_its,threshvals);

probmat=zeros(N,N);
    timemat=zeros(N,N);

    %size(twomat)=[Nl,Nd,Nt,Nt,2];

    for i=1:N
        for j=[1:(i-1),(i+1):N]
            if dmat(i,j)<=Nd-1
                probmat(i,j)=1-twomat(l,Nd-dmat(i,j)+1,Tvals(j),Tvals(i),1);
                timemat(i,j)=twomat(l,Nd-dmat(i,j)+1,Tvals(j),Tvals(i),2);
            else
                probmat(i,j)=twomat(l,dmat(i,j)-Nd+1,Tvals(i),Tvals(j),1);
                timemat(i,j)=twomat(l,dmat(i,j)-Nd+1,Tvals(i),Tvals(j),2);
            end
        end
    end

    binsigmat=zeros(N,N);

    for i=1:N
        for j=(i+1):N
            draw=rand;
            if draw <=probmat(i,j)
                binsigmat(i,j)=1;
            else binsigmat(j,i)=1;
            end
        end
    end

    sigmat=binsigmat.*(max(max(timemat))+1-timemat);
    keepsigmat(:,:,r,ic)=sigmat;
    
    todivide=sum(sigmat,2);
    todivide(todivide==0)=1;
    entmat=sigmat./repmat(todivide,1,N);
    logmat=entmat;
    logmat(entmat~=0)=-log(entmat(entmat~=0));
    entmat=entmat.*logmat;

%         powervec=sum(probmat,2);
    powervecD=col(sum(binsigmat,2));
    powervecR=col(sum(sigmat,2));
    powervecDelta=col((powervecD).*(powervecR));
    powervecH=(sum(entmat,2));
    powervecH=col(real(powervecH));
    powervecPi=col((powervecR).*(powervecH));
    powervecC=col(eigcentrality(sigmat,weight));
    powervecprobs=col(sum(probmat,2));
    
    powermatR(:,r,ic)=powervecR;
    powermatC(:,r,ic)=powervecC;
    powermatH(:,r,ic)=powervecH;
    powermatD(:,r,ic)=powervecD;
    
    q=quantile(fighting_abilities,0:1/Nb1:1);
        diffvec=diff(q);
        eps=min(diffvec(diffvec>=1e-5))/10;
        q(1)=q(1)-eps;q(end)=q(end)+eps;
        abilitiesbins=q;
        [~,abilitiesindices]=histc(fighting_abilities,abilitiesbins);
    
    powervec=powervecR;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    eps=min(diffvec(diffvec>=1e-5))/10;
    q(1)=q(1)-eps;q(end)=q(end)+eps;
    powerbins=q;
    [~,powerindicesR(:,r,ic)]=histc(powervec,powerbins);

    powervec=powervecC;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    eps=min(diffvec(diffvec>=1e-5))/10;
    q(1)=q(1)-eps;q(end)=q(end)+eps;
    powerbins=q;
    [~,powerindicesC(:,r,ic)]=histc(powervec,powerbins);
    
    powervec=powervecH;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    eps=min(diffvec(diffvec>=1e-5))/10;
    q(1)=q(1)-eps;q(end)=q(end)+eps;
    powerbins=q;
    [~,powerindicesH(:,r,ic)]=histc(powervec,powerbins);
    
    powervec=powervecD;
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
%     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    q=quantile(powervec,0:1/Nb2:1);
    diffvec=diff(q);
    eps=min(diffvec(diffvec>=1e-5))/10;
    q(1)=q(1)-eps;q(end)=q(end)+eps;
    powerbins=q;
    [~,powerindicesD(:,r,ic)]=histc(powervec,powerbins);
    
%     plot(1:20,[powervecR/max(powervecR) powervecC/max(powervecC)],'-o')
% plot(1:20,[powervecR/powervecR(1) powervecC/powervecC(1) powervecD/powervecD(1)])
% plot(powervecR/powervecR(end),powervecC/powervecC(end),'-o','Color',mycols{ic})
% plot(powervecR,powervecC,'-o','Color',mycols{ic})
% plot(powerindicesR,powerindicesC,'o','Color',mycols{ic},'MarkerFaceColor',mycols{ic},'MarkerAlpha',.1)
% for i=1:N
%     patch(powerindicesR(i,r)+.5*sin(theta),powerindicesC(i,r)+.5*cos(theta),mycols{ic},'FaceAlpha',1/runs*5,'EdgeColor','none')
% end
[~,Ro]=sort(powervecR);
[~,Co]=sort(powervecC);
[~,Do]=sort(powervecD);
% plot(1:20,Ro,'-','Color',mycols{ic})
% plot(1:20,Co,'--','Color',mycols{ic})
end
end

%%
i=1;
fc=find(powerindicesC(end-2,:,i)==1);
fr=find(powerindicesR(end-2,:,i)==1);
% fc=find(powerindicesC(end-5,:,i)==2);
% fr=find(powerindicesR(end-5,:,i)==2);
f=setdiff(fc,fr);
%%
i=2;
fc=find(powerindicesC(end,:,i)==1);
f2=find(powerindicesH(end,:,i)==1);
f=setdiff(fc,f2);

%%
i=3;
fc=find(powerindicesC(6,:,i)==9);
fr=find(powerindicesR(6,:,i)==9);
f=setdiff(fr,fc);

%%
i=4;
fd=find(powerindicesD(1,:,i)<10);
f2=find(powerindicesR(1,:,i)<10);
f=setdiff(f2,fd);

%%
i=5;
fh=find(powerindicesH(end,:,i)>2);
f2=find(powerindicesR(end,:,i)>2);
f=setdiff(f2,fh);

%%
i=6;
fr=find(powerindicesR(end,:,i)==1);
f2=find(powerindicesH(end,:,i)==1);
f=setdiff(fr,f2);