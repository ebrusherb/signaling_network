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
mycolormap=cbrewer('seq','Reds',64);
mycolormap=mycolormap(1:end,:);
cornerlabels={'Error rate','Personal preference','Decision time'};
alphabet=['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h'];

load /Users/eleanorbrush/Dropbox/signaling_network/twomat.mat
load /Users/eleanorbrush/Dropbox/signaling_network/variables
load /Users/eleanorbrush/Dropbox/signaling_network/optanalysis_output.mat
load /Users/eleanorbrush/Dropbox/signaling_network/infogrowth_output.mat
load twomat.mat
load onemat.mat
load variables.mat
load /Users/eleanorbrush/Dropbox/signaling_network/groupsize.mat
load /Users/eleanorbrush/Dropbox/signaling_network/homogenabilities.mat

%% coordinate set up
N=20;
its=1000;
cstep=.1;
vec=0:cstep:1;
Nc=length(vec);
endvals=zeros(1,Nc);

% X=repmat(0:cstep:1,1,11);
% Y=reshape(repmat(0:cstep:1,11,1),1,[]);

X=[];
Y=[];
layers=cell(Nc,1);
for i=0:(Nc-1)
    layers{i+1}=length(X)+(1:(Nc-i));
    endvals(i+1)=layers{i+1}(end);
    X=[X, vec(1:(Nc-i))];
    Y=[Y,vec(i+1)*ones(1,Nc-i)];   
end
Z=1-X-Y;

transformed=transform([X',Z',Y']);
Nc=length(X);

indices=1:Nt;
leak=2;
%% figure: dimensionality comparison 
load /Users/eleanorbrush/Dropbox/signaling_network/variables
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83/2;
h=w;
set(gcf,'Units','inches');
set(gcf,'Position',[v(1) v(2) w h]);

textlabels={'(a)' '(b)'};

p=zeros(1,2*Nd);

currblue=mycolors(Nd,'blue');
currred=mycolors(Nd,'red');
currblack=mycolors(Nd,'black');
currblue=cbrewer('seq','Blues',Nd);
currred=cbrewer('seq','Reds',Nd);
currblack=cbrewer('seq','Greys',Nd);

for k=2
%     subplot_tight(1,4,(1:2)+2*(k-1),[.1 .07])
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
set(gca,'FontName',fontname,'FontSize',labfontsz)
set(gca,'Position',[.15 .1417 .64 .7833])
end

domlabs=cell(1,Nd-1);
for i=1:(Nd-1)
    domlabs{i}=num2str(domvals(i+1));
end
domlabs{1}='c = 0.55';

leg=legend([p(Nd-4) p(Nd+Nd-4) p(Nd+(2:Nd))],horzcat('Two','One',domlabs));
legend boxoff
leg_line=findobj(leg,'type','Line');
for i = 1:(Nd-1)
    set(leg_line(2*i),'Color',currblack(Nd+1-i,:))
    v=get(leg_line(2*i),'Xdata');
    set(leg_line(2*i),'Xdata',[.3 v(2)])
%      set(leg_line(2*1), 'Color', currblack(Nd,:));
%      set(leg_line(2*1-1),'Color',currblack(Nd,:));
%      set(leg_line(2*2), 'Color', currblack(2,:));
%      set(leg_line(2*2-1),'Color',currblack(2,:));
end
v=get(leg_line(2*Nd),'Xdata');
set(leg_line(2*Nd),'xdata',[.3 v(2)]);
v=get(leg_line(2*(Nd+1)),'Xdata');
set(leg_line(2*(Nd+1)),'xdata',[.3 v(2)]);

set(leg,'Position',[.82 .55 .1 .1])
    
set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','dimensionality_comparison','.pdf');
print(filename,'-dpdf','-r300');

%% figure: effect of group size

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=w/3;
set(gcf,'Units','inches');
set(gcf,'Position',[v(1) v(2) w h]);

coloffset=4;
mycols=cbrewer('seq','Reds',Nlays+coloffset);
cvec=0:(1/(Nlays-1)):1;
leglabs=cell(1,Nlays);
for i=1:Nlays
    leglabs{i}=num2str(cvec(i));
end
leglabs{1}=strcat('w_2 = 0');
textlabs={'(a)','(b)','(c)'};

ylabs={'Average accuracy','Average accuracy','Average accuracy'};

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
    text(-9, 1.07,textlabs{k},'HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz,'FontName',fontname)
    set(gca,'Position',[.08+.3*(k-1)  .17 .2 .7])
    switch k
        case 3
            leg=legend(p,leglabs);
            legend('boxoff')
            set(leg,'Position',[.9 .1489 .0659 .7630])
    end
end

leg_line=findobj(leg,'type','Line');
for i = 1:Nlays
    v=get(leg_line(2*i),'Xdata');
    set(leg_line(2*i),'Xdata',[.3 v(2)])
%      set(leg_line(2*1), 'Color', currblack(Nd,:));
%      set(leg_line(2*1-1),'Color',currblack(Nd,:));
%      set(leg_line(2*2), 'Color', currblack(2,:));
%      set(leg_line(2*2-1),'Color',currblack(2,:));
end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','group_size','.pdf');
print(filename,'-dpdf','-r300');

%% figure: dimensionality  comparison with correct biased initial conditions
load /Users/eleanorbrush/Desktop/solve_pde_init_v2.mat
d=x0-y0;
f=find(threshvals>=d);
Nf1=length(f:Nt);
figure
set(gcf,'Color','w')
set(gcf,'Units','inches');
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83/3;
h=w;
% set(gcf,'Position',[v(1) v(2) w h]);

textlabels={'(a)' '(b)'};

p=zeros(1,2*Nd);

currblue=mycolors(Nd,'blue');
currred=mycolors(Nd,'red');
currblack=mycolors(Nd,'black');
currblue=cbrewer('seq','Blues',Nd);
currred=cbrewer('seq','Reds',Nd);
currblack=cbrewer('seq','Greys',Nd);


%     subplot_tight(1,4,(1:2)+2*(k-1),[.1 .07])
    set(gca,'FontSize',labfontsz)
    hold on
 
            for i=2:Nd
            p(i+Nd)=plot(reshape(diag(reshape(onemat_par(i,f:end,f:end,2),Nf1,Nf1)),1,[]),reshape(diag(reshape(onemat_par(i,f:end,f:end,1),Nf1,Nf1)),1,[]),'-','Color',currred(i,:),'LineWidth',lw);
            p(i)=plot(reshape(diag(reshape(twomat_par(i,:,:,2),Nt,Nt)),1,[]),reshape(diag(reshape(twomat_par(i,:,:,1),Nt,Nt)),1,[]),'--','Color',currblue(i,:),'LineWidth',lw);            
            end
        
    set(gca,'xlim',[0 5],'ylim',[.5 1]);
    ylabel('Probability of correct decision','FontSize',textfontsz,'FontName',fontname)
    xlabel('Expected time to decision','FontSize',textfontsz,'FontName',fontname)
set(gca,'FontName',fontname,'FontSize',labfontsz)
set(gca,'Position',[.15 .1417 .64 .7833])


domlabs=cell(1,Nd-1);
for i=1:(Nd-1)
    domlabs{i}=num2str(domvals(i+1));
end
domlabs{1}='c = 0.55';

domlabs=cell(1,Nd-1);
for i=1:(Nd-1)
    domlabs{i}=num2str(domvals(i+1));
end
domlabs{1}='c = 0.55';

leg=legend([p(Nd-4) p(Nd+Nd-4) p(Nd+(2:Nd))],horzcat('Two','One',domlabs));
legend boxoff
leg_line=findobj(leg,'type','Line');
for i = 1:(Nd-1)
    set(leg_line(2*i),'Color',currblack(Nd+1-i,:))
    v=get(leg_line(2*i),'Xdata');
    set(leg_line(2*i),'Xdata',[.3 v(2)])
%      set(leg_line(2*1), 'Color', currblack(Nd,:));
%      set(leg_line(2*1-1),'Color',currblack(Nd,:));
%      set(leg_line(2*2), 'Color', currblack(2,:));
%      set(leg_line(2*2-1),'Color',currblack(2,:));
end
v=get(leg_line(2*Nd),'Xdata');
set(leg_line(2*Nd),'xdata',[.3 v(2)]);
v=get(leg_line(2*(Nd+1)),'Xdata');
set(leg_line(2*(Nd+1)),'xdata',[.3 v(2)]);

set(leg,'Position',[.82 .55 .1 .1])
    
set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','dimensionality_comparison_correctbias','.pdf');
print(filename,'-dpdf','-r300');
%% figure: dimensionality  comparison with wrong biased initial conditions
load /Users/eleanorbrush/Desktop/solve_pde_init_v3.mat
d=x0-y0;
f=find(threshvals>=d);
Nf1=length(f:Nt);
figure
set(gcf,'Color','w')
set(gcf,'Units','inches');
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83/2;
h=w;
% set(gcf,'Position',[v(1) v(2) w h]);

textlabels={'(a)' '(b)'};

p=zeros(1,2*Nd);

currblue=mycolors(Nd,'blue');
currred=mycolors(Nd,'red');
currblack=mycolors(Nd,'black');
currblue=cbrewer('seq','Blues',Nd);
currred=cbrewer('seq','Reds',Nd);
currblack=cbrewer('seq','Greys',Nd);


%     subplot_tight(1,4,(1:2)+2*(k-1),[.1 .07])
    set(gca,'FontSize',labfontsz)
    hold on
 
            for i=2:Nd
            p(i+Nd)=plot(reshape(diag(reshape(onemat_par(i,f:end,f:end,2),Nf1,Nf1)),1,[]),reshape(diag(reshape(onemat_par(i,f:end,f:end,1),Nf1,Nf1)),1,[]),'-','Color',currred(i,:),'LineWidth',lw);
            p(i)=plot(reshape(diag(reshape(twomat_par(i,:,:,2),Nt,Nt)),1,[]),reshape(diag(reshape(twomat_par(i,:,:,1),Nt,Nt)),1,[]),'--','Color',currblue(i,:),'LineWidth',lw);            
            end
        
    set(gca,'xlim',[0 5]);
    ylabel('Probability of correct decision','FontSize',textfontsz,'FontName',fontname)
    xlabel('Expected time to decision','FontSize',textfontsz,'FontName',fontname)
set(gca,'FontName',fontname,'FontSize',labfontsz)
set(gca,'Position',[.15 .1417 .64 .7833])


domlabs=cell(1,Nd-1);
for i=1:(Nd-1)
    domlabs{i}=num2str(domvals(i+1));
end
domlabs{1}='c = 0.55';

domlabs=cell(1,Nd-1);
for i=1:(Nd-1)
    domlabs{i}=num2str(domvals(i+1));
end
domlabs{1}='c = 0.55';

leg=legend([p(Nd-4) p(Nd+Nd-4) p(Nd+(2:Nd))],horzcat('Two','One',domlabs));
legend boxoff
leg_line=findobj(leg,'type','Line');
for i = 1:(Nd-1)
    set(leg_line(2*i),'Color',currblack(Nd+1-i,:))
    v=get(leg_line(2*i),'Xdata');
    set(leg_line(2*i),'Xdata',[.3 v(2)])
%      set(leg_line(2*1), 'Color', currblack(Nd,:));
%      set(leg_line(2*1-1),'Color',currblack(Nd,:));
%      set(leg_line(2*2), 'Color', currblack(2,:));
%      set(leg_line(2*2-1),'Color',currblack(2,:));
end
v=get(leg_line(2*Nd),'Xdata');
set(leg_line(2*Nd),'xdata',[.3 v(2)]);
v=get(leg_line(2*(Nd+1)),'Xdata');
set(leg_line(2*(Nd+1)),'xdata',[.3 v(2)]);

set(leg,'Position',[.82 .55 .1 .1])
    
set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','dimensionality_comparison_wrongbias','.pdf');
print(filename,'-dpdf','-r300');

%% extra skewness heatmaps
numfun_toplot=[2 4 6];
mainlabel='Skewness';

figure
% set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=1/3*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)
marg=[.1 .063];
xpos=.3;
exstofind=[layers{1}(1) layers{5}(1) layers{7}(1)];
colormap(mycolormap);
width=.3;

for i=1:3
  subplot(1,3,i)  
    numfun=numfun_toplot(i);
    [mybar,yl]=triimage_subplot(skewvals(numfun,:),transformed,cornerlabels,mainlabel,'off',exstofind,doubletrilw,{'b','c','d'},'k');
    v=get(mybar,'Position');
    set(mybar,'Position',[.25+width*(i-1) v(2)+.03 .01 v(4)])
    v=get(yl,'Position');
    set(yl,'Position',[21 v(2) v(3)]);
    v=get(gca,'Position');
    set(gca,'Position',[.05+width*(i-1),v(2) .18 v(4)]);
    v=get(gca,'Position');
    annotation('textbox',[v(1)-.03 .93 .01 .04],'String',strcat('(',alphabet(i),')'),'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','supp_skewness','.pdf');
print(filename,'-dpdf','-r300');
