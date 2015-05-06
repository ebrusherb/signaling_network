%%
axislw=1.25;
lw=3;
trilw=2;
doubletrilw=1;
markersz=3;
textfontsz=12;
labfontsz=9;
fontname='Times New Roman';
yellow=[1 1 0];
% red=[.9 0 0];
% green=[0 .7 0];
% blue=[0 0 .9];
% greens=importdata('green.txt');
% blues=importdata('blue.txt');
% reds=importdata('red.txt');
% blacks=importdata('black.txt');

load /Users/eleanorbrush/Dropbox/signaling_network/twomat.mat
load /Users/eleanorbrush/Dropbox/signaling_network/variables
load /Users/eleanorbrush/Dropbox/signaling_network/optanalysis_output.mat
load /Users/eleanorbrush/Desktop/infogrowth_output.mat
load twomat.mat
load onemat.mat
load variables.mat
load /Users/eleanorbrush/Desktop/groupsize.mat


%% figure: dimensionality comparison 

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83*2;
h=w/3;
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
    ylabel('Probability of correct decision','FontSize',textfontsz,'FontName',fontname)
    xlabel('Expected time to decision','FontSize',textfontsz,'FontName',fontname)
text(-.4, 1.03,textlabels(k),'HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz,'FontName',fontname)
set(gca,'FontName',fontname,'FontSize',labfontsz)
set(gca,'Position',[.1+.43*(k-1) .1 .38 .8])
end

domlabs=cell(1,Nd-1);
for i=1:(Nd-1)
    domlabs{i}=num2str(domvals(i+1));
end
domlabs{1}='c = .55';

leg=legend([p(Nd) p(Nd+Nd) p(Nd+(2:Nd))],horzcat('Two','One',domlabs));
legend boxoff
leg_line=findobj(leg,'type','Line');
for i = 1:(Nd-1)
    set(leg_line(2*i),'Color',currblack(Nd+1-i,:))
    v=get(leg_line(2*i),'Xdata');
    set(leg_line(2*i),'Xdata',[.4 v(2)])
%      set(leg_line(2*1), 'Color', currblack(Nd,:));
%      set(leg_line(2*1-1),'Color',currblack(Nd,:));
%      set(leg_line(2*2), 'Color', currblack(2,:));
%      set(leg_line(2*2-1),'Color',currblack(2,:));
end
v=get(leg_line(2*Nd),'Xdata');
set(leg_line(2*Nd),'xdata',[.4 v(2)]);
v=get(leg_line(2*(Nd+1)),'Xdata');
set(leg_line(2*(Nd+1)),'xdata',[.4 v(2)]);

set(leg,'Position',[.9 .55 .1 .1])
    
set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','dimensionality_comparison','.pdf');
print(filename,'-dpdf','-r300');

%% figure: effect of group size

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83*2;
h=.25*w;
set(gcf,'Units','inches');
set(gcf,'Position',[v(1) v(2) w h]);

coloffset=4;
mycols=cbrewer('seq','Reds',Nlays+coloffset);
cvec=0:(1/(Nlays-1)):1;
leglabs=cell(1,Nlays);
for i=1:Nlays
    leglabs{i}=num2str(cvec(i));
end
leglabs{1}=strcat('c_2 = 0');
textlabs={'(a)','(b)','(c)'};

ylabs={'Average accuracy','Average accuracy of top quartile','Average accuracy of bottom quartile'};

accuracy_mat=cell(3,1);
accuracy_mat{1}=accuracy;
accuracy_mat{2}=topaccuracy;
accuracy_mat{3}=bottomaccuracy;

for k=1:3
    subplot(1,3,k)
    hold on
    p=zeros(1,Nlays);
    for j=1:Nlays
        p(j)=plot(sizes,accuracy_mat{k}(:,j,1),'Color',mycols(j+coloffset,:),'LineWidth',lw);
    %     plot(sizes,accuracy(:,j,2),'--','Color',mycols(j+coloffset,:),'LineWidth',lw)
    end
    set(gca,'ylim',[.5 1])
    set(gca,'FontSize',labfontsz,'FontName',fontname)
    xlabel('Group Size','FontName',fontname,'FontSize',textfontsz)
    ylabel(ylabs{k},'FontName',fontname,'FontSize',textfontsz)
    text(-8, 1.03,textlabs{k},'HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz,'FontName',fontname)
    set(gca,'Position',[.05+.3*(k-1)  .11 .23 .815])
    switch k
        case 3
            leg=legend(p,leglabs);
            legend('boxoff')
            set(leg,'Position',[.9 .1489 .0659 .7630])
    end
end
% 
% subplot(1,3,2)
% hold on
% p=zeros(1,Nlays);
% for j=1:Nlays
%     p(j)=plot(sizes,topaccuracy(:,j,1),'Color',mycols(j+coloffset,:),'LineWidth',lw);
% %     plot(sizes,topaccuracy(:,j,2),'--','Color',mycols(j+coloffset,:),'LineWidth',lw)
% end
% set(gca,'FontSize',labfontsz,'FontName',fontname)
% xlabel('Group Size','FontName',fontname,'FontSize',textfontsz)
% ylabel('Average accuracy of top quartile','FontName',fontname,'FontSize',textfontsz)
% text(-7, 1.02,'(b)','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz,'FontName',fontname)
% 
% 
% subplot(1,3,3)
% hold on
% p=zeros(1,Nlays);
% for j=1:Nlays
%     p(j)=plot(sizes,bottomaccuracy(:,j,1),'Color',mycols(j+coloffset,:),'LineWidth',lw);
% %     plot(sizes,bottomaccuracy(:,j,2),'--','Color',mycols(j+coloffset,:),'LineWidth',lw)
% end
% set(gca,'FontSize',labfontsz,'FontName',fontname)
% xlabel('Group Size','FontName',fontname,'FontSize',textfontsz)
% ylabel('Average accuracy of bottom quartile','FontName',fontname,'FontSize',textfontsz)
% text(-7, 1.02,'(c)','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz,'FontName',fontname)
% leg=legend(p,leglabs);
% legend('boxoff')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','group_size','.pdf');
print(filename,'-dpdf','-r300');

