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
lvals=[0 .5];
Nl=length(lvals);

domvals=.5:.1:1;
Nd=length(domvals);

threshvals=.5:.5:2;
Nt=length(threshvals);

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
Nc=2;
c1=0;
c2=.1;
c3=1;

perf=zeros(Nc,2,Nl,Nd,Nt,Nt); %criteria, individual, leak, dominance, thresholds

perf(1,1,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2));
perf(1,2,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2));
perf(2,1,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2))+c3*(1-twomat(:,2:Nd,:,:,1));
perf(2,2,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2))+c3*(twomat(:,2:Nd,:,:,1));

perf(1,1,:,1,:,:)=c2*(twomat(:,1,:,:,2));
perf(1,2,:,1,:,:)=c2*(twomat(:,1,:,:,2));
perf(2,1,:,1,:,:)=c2*(twomat(:,1,:,:,2))+c3*(1-twomat(:,1,:,:,1));
perf(2,2,:,1,:,:)=c2*(twomat(:,1,:,:,2))+c3*(twomat(:,1,:,:,1));

nasheq=cell(Nc,Nl,Nd);

T1_best=zeros(Nc,Nl,Nd,Nt); %criteria, leak, dominance, threshold
T2_best=zeros(Nc,Nl,Nd,Nt);

for i=1:Nc
    for j=1:Nl
        for k=1:Nd
            for u=1:Nt
            [~,m]=min(perf(i,1,j,k,:,u));
            T1_best(i,j,k,u)=m;
            [~,m]=min(perf(i,2,j,k,u,:));
            T2_best(i,j,k,u)=m;
            end
            pp1=interp1(threshvals,reshape(threshvals(T1_best(i,j,k,:)),1,[]),'linear','pp');
            pp2=interp1(threshvals,reshape(threshvals(T2_best(i,j,k,:)),1,[]),'linear','pp');
            v1=ppval(pp1,threshvals);
            v2=ppval(pp2,v1);
            x=threshvals(abs(threshvals-v2)<0.01);
            y=ppval(pp1,x);
            nasheq{i,j,k}=[x;y];
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

p=zeros(1,Nl*2);

for i=Nl:Nl
    for j=1:Nd
        L=size(nasheq{1,i,j},2);
        d=domvals(j);
        num=plot(d*ones(1,L),nasheq{2,i,j}(1,:),'o','Color',blue(6-i,:),'LineWidth',lw,'MarkerSize',markersz+5);
        p(2*i)=num(1);
        num=plot(d*ones(1,L),nasheq{2,i,j}(2,:),'o','Color',red(6-i,:),'LineWidth',lw);
        p(2*i-1)=num(1);
    end
end

p(end-1)=p(1);
p(end)=p(1);

xlabel('Probability of stronger winning a fight','FontSize',textfontsz)
ylabel('Nash equilibrium threshold','FontSize',textfontsz)

leg=legend([p(1) p(2)],'Weaker','Stronger');
set(leg,'Position',[.85 .5 .1 .1])
legend boxoff
% leg_line=findobj(leg,'type','Line');
% for i = 4+1:Nl
%      set(leg_line(2*i), 'Color', black(i+1,:));
%      set(leg_line(2*i-1),'Color',black(i+1,:));
% end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','2d_nash_eq_thresholds','.pdf');
print(filename,'-dpdf','-r300');
