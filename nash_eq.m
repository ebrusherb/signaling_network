green=importdata('green.txt');
blue=importdata('blue.txt');
red=importdata('red.txt');
black=importdata('black.txt');
%%
axislw=1.25;
lw=3;
labfontsz=15 ;
textfontsz=15;
alpha=.7;
ticklength=[.02,.0250];
markersz=5;
markercol='white';
%%
lvals=0:.1:1;
Nl=length(lvals);

domvals=.5:.05:1;
Nd=length(domvals);

threshvals=.5:.1:2;
Nt=length(threshvals);

indices=1:Nt;

b=1;
deltax2=0.025;

threshvals1=.1:.1:8;
Nt1=length(threshvals1);

%%
twomat=zeros(Nl,Nd,Nt,Nt,2); %leak, dominance, thresholds, prob, time

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
onemat=zeros(Nl,Nd,Nt1,2);

for k=1:Nl
    l=lvals(k);
    for i=1:Nd
        d=domvals(i);
        for j=1:Nt1
            t=threshvals1(j);
            [x,y,~]=solve_pde_1d(-t,d,t,deltax1,l,b);
            onemat(k,i,j,:)=[x y];
        end
    end
end
%%
filename=strcat('twomat_long_run','.mat');
save(filename,'twomat','lvals','domvals','threshvals','threshvals');
%%

c1=0;
c3=1;

c2vals=0:.1:1;
Nc2=length(c2vals);

nasheq_ind=cell(Nc2,Nl,Nd);
giveup=zeros(Nc2,Nl,Nd);

for i=1:Nc2
    c2=c2vals(i);
    perf=zeros(2,Nl,Nd,Nt,Nt); %criteria, individual, leak, dominance, thresholds

    perf(1,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2))+c3*(1-twomat(:,2:Nd,:,:,1));
    perf(2,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2))+c3*(twomat(:,2:Nd,:,:,1));

    perf(1,:,1,:,:)=c2*(twomat(:,1,:,:,2))+c3*(1-twomat(:,1,:,:,1));
    perf(2,:,1,:,:)=c2*(twomat(:,1,:,:,2))+c3*(twomat(:,1,:,:,1));
    for j=1:Nl
        for k=1:Nd
            T1_best=zeros(1,Nt); 
            T2_best=zeros(1,Nt);
            for u=1:Nt
            [~,m]=min(perf(1,j,k,:,u));
            T1_best(u)=m;
            [~,m]=min(perf(2,j,k,u,:));
            T2_best(u)=m;
            end
            pp1=interp1(indices,reshape(T1_best,1,[]),'linear','pp');
            pp2=interp1(indices,reshape(T2_best,1,[]),'linear','pp');
            v1=ppval(pp2,indices);
            v2=ppval(pp1,v1);
            T1=indices(abs((indices)-v2)<0.01);
            T2=ppval(pp2,T1);
            nasheq_ind{i,j,k}=[T1;T2];
            if size(nasheq_ind{i,j,k},2)==1
                giveup(i,j,k)=(sum(nasheq_ind{i,j,k}==[Nt;1])==2);
            end
        end
    end
end

%%

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83*2;
h=w/2;
set(gcf,'Units','inches');
set(gcf,'Position',[v(1) v(2) w h]);
hold on
set(gca,'FontSize',labfontsz)



filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','when_to_give_up','.pdf');
print(filename,'-dpdf','-r300');
