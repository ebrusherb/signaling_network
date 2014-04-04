lvals=[0 .5];
Nl=length(lvals);

domvals=.5:.1:1;
Nd=length(domvals);

threshvals=.5:.5:2;
Nt=length(threshvals);

twomat=zeros(Nl,Nd,Nt,Nt,2);

for k=1:Nl
    l=lvals(k);
    for u=1:Nd
        d=domvals(u);
        for i=1:Nt
            T1=threshvals(i);
            for j=1:Nt
                T2=threshvals(j);
                [x,y,~]=solve_pde_2d(-T1,-T2,d,5,5,deltax2,deltax2,l,b);
                twomat(k,u,i,j,:)=[x y];
            end
        end
    end
end

%%
Nc=2;

perf=zeros(Nc,2,Nl,Nd,Nt,Nt); %criteria, individual, leak, dominance, thresholds

perf(1,:,:,:,:,:)=twomat(:,:,:,:,1)+twomat(:,:,:,:,2);
perf(2,1,:,:,:,:)=twomat(:,:,:,:,1)+twomat(:,:,:,:,2)+.5*(1-twomat(:,:,:,:,1));
perf(2,2,:,:,:,:)=twomat(:,:,:,:,1)+twomat(:,:,:,:,2)+.5*(twomat(:,:,:,:,1));

nasheq=zeros(Nc,Nl,Nd);

T1_best=zeros(Nc,NL,Nd,Nt);
T2_best=zeros(Nc,Nl,Nd,Nt);

for i=1:Nc
    for j=1:NL
        for k=1:Nd
            for u=1:Nt
            [~,m]=max(perf(i,1,j,k,:,u));
            T1_best(i,j,k,u)=m;
            [~,m]=max(perf(i,2,j,k,u,:));
            T2_best(i,j,k,u)=m;
            end
        end
    end
end
%%

NE1={};
NE2={};

for i=1:Nc
    for j=1:NL
        for k=1:Nd
    pp1=interp1(Lyvals,T1_best(:,1+i,1),'linear','pp');
    pp2=interp1(Lyvals,T1_best(:,Nd-i,1),'linear','pp');
    Tystar=Lyvals(abs(ppval(pp2,ppval(pp1,Lyvals))-Lyvals)<0.001);
    Txstar=ppval(pp1,Tystar);
    NE1{i+1}=Txstar;
    NE2{i+1}=Tystar;
        end
    end
end


