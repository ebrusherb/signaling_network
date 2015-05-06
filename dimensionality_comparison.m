greens=importdata('green.txt');
blues=importdata('blue.txt');
reds=importdata('red.txt');
blacks=importdata('black.txt');
load twomat.mat
load onemat.mat
load variables.mat

axislw=1.25;
lw=3;
labfontsz=15 ;
textfontsz=15;
fontname='Times';
alpha=.7;
ticklength=[.02,.0250];
markersz=5;
markercol='white';

%%

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83*2;
h=w/2;
set(gcf,'Units','inches');
set(gcf,'Position',[v(1) v(2) w h]);

textlabels={'(a)' '(b)'};

p=zeros(1,2*Nd);

currblue=mycolors(Nd,'blue');
currred=mycolors(Nd,'red');
currblack=mycolors(Nd,'black');
% currblue=cbrewer('seq','Blues',Nd);
% currred=cbrewer('seq','Reds',Nd);
% currblack=cbrewer('seq','Grays',Nd);

for k=1:2
    subplot_tight(1,4,(1:2)+2*(k-1),[.1 .07])
    set(gca,'FontSize',labfontsz)
    hold on
    switch k
        case 2
            for i=2:Nd
            p(i+Nd)=plot(reshape(diag(reshape(onemat(k,i,:,:,2),Nt1,[])),1,[]),reshape(diag(reshape(onemat(k,i,:,:,1),Nt1,[])),1,[]),'Color',currred(i,:),'LineWidth',lw);
            p(i)=plot(reshape(diag(reshape(twomat(k,i,:,:,2),Nt,[])),1,[]),reshape(diag(reshape(twomat(k,i,:,:,1),Nt,[])),1,[]),'--','Color',currblue(i,:),'LineWidth',lw);            
            end
        otherwise 
            for i=2:Nd
            plot(reshape(diag(reshape(onemat(k,i,:,:,2),Nt1,[])),1,[]),reshape(diag(reshape(onemat(k,i,:,:,1),Nt1,[])),1,[]),'Color',currred(i,:),'LineWidth',lw);
            plot(reshape(diag(reshape(twomat(k,i,:,:,2),Nt,[])),1,[]),reshape(diag(reshape(twomat(k,i,:,:,1),Nt,[])),1,[]),'--','Color',currblue(i,:),'LineWidth',lw);            
            end
    end
    set(gca,'xlim',[0 5],'ylim',[.5 1]);
    ylabel('Probability of correct output','FontSize',textfontsz,'FontName',fontname)
    xlabel('Expected time to output','FontSize',textfontsz,'FontName',fontname)
text(-.4, 1.03,textlabels(k),'HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz,'FontName',fontname)
set(gca,'FontName',fontname,'FontSize',labfontsz)
set(gca,'Position',[.1+.43*(k-1) .1 .33 .8])
end

domlabs=cell(1,Nd-1);
for i=1:(Nd-1)
    domlabs{i}=num2str(domvals(i+1));
end

leg=legend([p(Nd) p(Nd+Nd) p(Nd+(2:Nd))],horzcat('Two','One',domlabs));
legend boxoff
leg_line=findobj(leg,'type','Line');
for i = 1:(Nd-1)
    set(leg_line(2*i),'Color',currblack(Nd+1-i,:))
%      set(leg_line(2*1), 'Color', currblack(Nd,:));
%      set(leg_line(2*1-1),'Color',currblack(Nd,:));
%      set(leg_line(2*2), 'Color', currblack(2,:));
%      set(leg_line(2*2-1),'Color',currblack(2,:));
end

set(leg,'Position',[.9 .5 .1 .1])
    
set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','dimensionality_comparison','.pdf');
% print(filename,'-dpdf','-r300');

%%
i=1;
j=(Nd-1);

ratios=[.5 1 2];
Nr=length(ratios);

twoflexthresh=zeros(Nr,Nt,2);
oneflexthresh=zeros(Nr,Nt1,2);

for k=1:Nr
    r=ratios(k);
    twoflexthresh(k,:,1)=interp2(threshvals,threshvals,reshape(twomat(i,j,:,:,1),Nt,[]),threshvals,r*threshvals);
    twoflexthresh(k,:,2)=interp2(threshvals,threshvals,reshape(twomat(i,j,:,:,2),Nt,[]),threshvals,r*threshvals);
    oneflexthresh(k,:,1)=interp2(threshvals1,threshvals1,reshape(onemat(i,j,:,:,1),Nt1,[]),threshvals1,r*threshvals1);
    oneflexthresh(k,:,2)=interp2(threshvals1,threshvals1,reshape(onemat(i,j,:,:,2),Nt1,[]),threshvals1,r*threshvals1);
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

p=zeros(1,Nr);
currblue=mycolors(Nr+1,'blue');
currred=mycolors(Nr+1,'red');

subplot_tight(2,2,1:2,[.1 .07])
set(gca,'FontSize',labfontsz)
hold on

for i=1:Nr
    p(i)=plot(reshape(twoflexthresh(i,:,2),1,[]),reshape(twoflexthresh(i,:,1),1,[]),'Color',currblue(i+1,:),'LineWidth',lw);
end

ylabel('Probability of correct output','FontSize',textfontsz)
xlabel('Expected time to output','FontSize',textfontsz)

legendlabels={};
for i=1:Nr
    legendlabels(i)=cellstr(['\rho=' num2str(ratios(i))]);
end

% leg=legend(p(1:5),['d=' num2str(ratios(1))],['d=' num2str(ratios(2))],['d=' num2str(ratios(3))],['d=' num2str(ratios(4))],['d=' num2str(ratios(5))]);
leg=legend(p,legendlabels);
legend boxoff
% leg_line=findobj(leg,'type','Line');
% for i = 1:Nr
%      set(leg_line(2*i), 'Color', black(i+1,:));
%      set(leg_line(2*i-1),'Color',black(i+1,:));
% end

set(leg,'Position',[.85 .65 .1 .1])

set(gca,'xlim',[.6 1]);

subplot_tight(2,2,3:4,[.1 .07])
set(gca,'FontSize',labfontsz)
hold on

for i=1:Nr
    p(i)=plot(reshape(oneflexthresh(i,:,2),1,[]),reshape(oneflexthresh(i,:,1),1,[]),'Color',currred(i+1,:),'LineWidth',lw);
end

ylabel('Probability of correct output','FontSize',textfontsz)
xlabel('Expected time to output','FontSize',textfontsz)

legendlabels={};
for i=1:Nr
    legendlabels(i)=cellstr(['\rho=' num2str(ratios(i))]);
end

% leg=legend(p(1:5),['d=' num2str(ratios(1))],['d=' num2str(ratios(2))],['d=' num2str(ratios(3))],['d=' num2str(ratios(4))],['d=' num2str(ratios(5))]);
leg=legend(p,legendlabels);
legend boxoff
% leg_line=findobj(leg,'type','Line');
% for i = 1:Nr
%      set(leg_line(2*i), 'Color', black(i+1,:));
%      set(leg_line(2*i-1),'Color',black(i+1,:));
% end

set(leg,'Position',[.85 .25 .1 .1])

set(gca,'xlim',[.6 1]);
    
set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','2d_flexible_thresholds','.pdf');
% print(filename,'-dpdf','-r300');