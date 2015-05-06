%% parameters
% the old version of this file was called figures_extra 
axislw=1.25;
lw=3;
trilw=2;
doubletrilw=1;
markersz=3;
textfontsz=8;
labfontsz=6;
fontname='Helvetica';
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
cornerlabels={'Error rate','Decision preference','Decision time'};
funtypes={'number of signals','entropy','number of signalers','eigenvector centrality'};
alphabet=['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h'];
fourcols=cbrewer('qual','Set1',4);
fourcols=fourcols([3 1 2 4],:);
markertypes={'-o','-d','-v','-s'};
funstoplot=[1 2 4 6];
colorperm=[1 2 4 3];

load /Users/eleanorbrush/Dropbox/signaling_network/twomat.mat
load /Users/eleanorbrush/Dropbox/signaling_network/variables
load /Users/eleanorbrush/Dropbox/signaling_network/optanalysis_output.mat
load /Users/eleanorbrush/Dropbox/signaling_network/infogrowth_output.mat
load twomat.mat
load onemat.mat
load variables.mat
load /Users/eleanorbrush/Dropbox/signaling_network/groupsize.mat
load /Users/eleanorbrush/Dropbox/signaling_network/homogenabilities.mat

% coordinate set up
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
%% figure: dimensionality comparison with biased initial conditions 
load variables.mat 

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
w=6.83;
h=1/3*w;
set(gcf,'Units','inches');
set(gcf,'Position',[v(1) v(2) w h]);
yanno=.95;


subplot(1,3,1)
set(gca,'FontName',fontname,'FontSize',labfontsz)
p=zeros(1,2*Nd);

currblue=cbrewer('seq','Blues',Nd);
currred=cbrewer('seq','Reds',Nd);
currblack=cbrewer('seq','Greys',Nd);

for k=2
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
  
end


[x,y]=data2norm(-1.5,.5);
annotation('textbox',[x yanno .05 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(1,3,2)
load solve_pde_init_v2.mat
twomat_correct=twomat_par;
onemat_correct=onemat_par;
d=x0-y0;
f=find(threshvals>=d);
Nf1=length(f:Nt);

p=zeros(1,2*Nd);

currblue=cbrewer('seq','Blues',Nd);
currred=cbrewer('seq','Reds',Nd);
currblack=cbrewer('seq','Greys',Nd);

set(gca,'FontSize',labfontsz)
hold on

for i=2:Nd
p(i+Nd)=plot(reshape(diag(reshape(onemat_correct(i,f:end,f:end,2),Nf1,Nf1)),1,[]),reshape(diag(reshape(onemat_correct(i,f:end,f:end,1),Nf1,Nf1)),1,[]),'-','Color',currred(i,:),'LineWidth',lw);
p(i)=plot(reshape(diag(reshape(twomat_correct(i,:,:,2),Nt,Nt)),1,[]),reshape(diag(reshape(twomat_correct(i,:,:,1),Nt,Nt)),1,[]),'--','Color',currblue(i,:),'LineWidth',lw);            
end

set(gca,'xlim',[0 5],'ylim',[.5 1]);
ylabel('Probability of correct decision','FontSize',textfontsz,'FontName',fontname)
xlabel('Expected time to decision','FontSize',textfontsz,'FontName',fontname)
set(gca,'FontName',fontname,'FontSize',labfontsz)
% set(gca,'Position',[.15 .1417 .64 .7833])
    
[x,y]=data2norm(-1.5,.5);
annotation('textbox',[x yanno .05 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')


subplot(1,3,3)
load solve_pde_init_v3.mat
twomat_wrong=twomat_par;
onemat_wrong=onemat_par;
d=x0-y0;
f=find(threshvals>=d);
Nf1=length(f:Nt);


p=zeros(1,2*Nd);

currblue=cbrewer('seq','Blues',Nd);
currred=cbrewer('seq','Reds',Nd);
currblack=cbrewer('seq','Greys',Nd);

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
% set(gca,'Position',[.15 .1417 .64 .7833])

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

set(leg,'Position',[.9 .55 .08 .1])

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);


[x,y]=data2norm(-1.5,.5);
annotation('textbox',[x yanno .05 .04],'String','(c)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')


filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','dimensionality_comparison','.pdf');
print(filename,'-dpdf','-r300');

filename='/Users/eleanorbrush/Desktop/test.pdf';
print(filename,'-dpdf','-r300');
%% figure: effect of group size
load variables.mat 
load twomat.mat 

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=w/3;
set(gcf,'Units','inches');
set(gcf,'Position',[v(1) v(2) w h]);

lw=2;
markersz=3;

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
        p(j)=plot(sizes(4:end),accuracy_mat{k}(4:end,j,1),'-o','Color',mycols(j+coloffset,:),'LineWidth',lw,'MarkerFaceColor',mycols(j+coloffset,:),'MarkerSize',markersz);
    end
    set(gca,'ylim',[.5 1])
    set(gca,'FontSize',labfontsz,'FontName',fontname)
    xlabel('Group Size','FontName',fontname,'FontSize',textfontsz)
    ylab=ylabel(ylabs{k},'FontName',fontname,'FontSize',textfontsz);
    v=get(ylab,'Position');
    set(ylab,'Position',[v(1) 875 v(3)]);
    text(-9, 1.02,textlabs{k},'HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz,'FontName',fontname)
    set(gca,'Position',[.08+.3*(k-1)  .17 .2 .7])
    set(gca,'ylim',[.6 1],'xlim',[10 101])
    set(gca,'tickdir','out')
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

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','group_size','.pdf');
print(filename,'-dpdf','-r300');

filename='/Users/eleanorbrush/Desktop/test.pdf';
print(filename,'-dpdf','-r300');
%% mutual information vs accuracy figure
load /Users/eleanorbrush/Dropbox/signaling_network/optanalysis_output.mat

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
w=3.4;
h=w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
set(gca,'FontSize',labfontsz,'FontName',fontname)
lw=1;

hold on
p=zeros(1,length(funstoplot));
for i=1:length(funstoplot)
    p(i)=plot((accuracy),mutvals(funstoplot(i),:),'LineWidth',lw,'Color',fourcols(colorperm(i),:));
end
xlabel('Average pairwise accuracy','FontSize',textfontsz);
ylabel('Mutual information','FontSize',textfontsz);

[leg,hobj]=legend(p,{'unweighted in-degree','weighted in-degree','entropy','eigenvector centrality'},'FontSize',textfontsz);
legend boxoff

v=get(leg,'Position');
set(leg,'Position',[.2 .7 v(3) v(4)])
textobj=findobj(hobj,'type','line');
for i=[1 3 5 7]
    set(textobj(i),'xdata',[.16 .3])
end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','mutinfo_vs_accuracy','.pdf');
print(filename,'-dpdf','-r300');


filename='/Users/eleanorbrush/Desktop/test.pdf';
print(filename,'-dpdf','-r300');
%% extra skewness heatmaps
numfun_toplot=[1 2 4 6];
mainlabel='Skewness';

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=1/3*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
set(gca,'FontSize',labfontsz,'FontName',fontname)
marg=[.1 .063];
xpos=.3;
exstofind=[layers{1}(1) layers{5}(1) layers{7}(1)];
colormap(mycolormap);
width=.23;

for i=1:4
  subplot(1,4,i)  
    numfun=numfun_toplot(i);
    [mybar,yl]=triimage_subplot(skewvals(numfun,:),transformed,cornerlabels,mainlabel,'off',exstofind,doubletrilw,{'b','c','d'},'k');
    v=get(mybar,'Position');
    set(mybar,'Position',[.18+width*(i-1) v(2)+.08 .01 v(4)-.1])
    v=get(yl,'Position');
    set(yl,'Position',[18 v(2) v(3)]);
    v=get(gca,'Position');
    set(gca,'Position',[.01+width*(i-1),v(2) .18 v(4)]);
    v=get(gca,'Position');
    annotation('textbox',[v(1)-.01 .93 .01 .04],'String',strcat('(',alphabet(i),')'),'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','multi_skewness','.pdf');
print(filename,'-dpdf','-r300');


filename='/Users/eleanorbrush/Desktop/test.pdf';
print(filename,'-dpdf','-r300');
%% skewness histograms

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.5*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)

numfun=1;

colormap(mycolormap);
m=min(skewvals(numfun,:));
M=max(skewvals(numfun,:));
xoffset=.25;
yoffset=.2;
subplot(3,4,[1:3 5:7 9:11])

[mybar,yl]=triimage_returnhandles(skewvals(numfun,:),transformed,cornerlabels,mainlabel,'on',exstofind,trilw,{'b','c','d'},zeros(4,3),m,M,1,xoffset,yoffset);
set(gca,'xlim',[0-2/11 2+2/11],'ylim',[0-1.7321/11 1.7321+1.7321/11])
v=get(gca,'Position');
scale=1;
set(gca,'Position',[.05 v(2)-.04 scale*v(3) scale*v(4)])
set(mybar,'Position',[.53 .15 .02 .7])
v=get(yl,'Position');
set(yl,'Position',[10.5 v(2) v(3)])

annotation('textbox',[.05 .95 .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

toplot=[1:3];
for I=1:3
    i=toplot(I);
    subplot(3,3,3*(I-1)+3)
%     hold on
    hist(expowervecs(:,numfun,i));
    histobj = findobj(gca,'Type','patch');
    set(histobj,'FaceColor',fourcols(1,:));
    set(gca,'FontName',fontname,'FontSize',labfontsz)
    xlim=get(gca,'xlim');
    ylim=get(gca,'Ylim');
    switch I
        case 1
            ylab=ylabel('Frequency','FontSize',textfontsz,'FontName',fontname);
            xlab=xlabel('Number of signalers','FontSize',textfontsz,'FontName',fontname);
            v=get(xlab,'Position');
            set(xlab,'Position',[v(1) ylim(1)-diff(ylim)*.18 v(3)]);
    end
    box off
    [x,y]=data2norm(xlim(1),ylim(2));
    [x2,y2]=data2norm(xlim(2),ylim(1));
    annotation('textbox',[x-.3*(x2-x) y+.01 .05 .04],'String',['(' alphabet(I+1) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','skewness_histograms','.pdf');
print(filename,'-dpdf','-r300');

filename='/Users/eleanorbrush/Desktop/test.pdf';
print(filename,'-dpdf','-r300');
%% line plot of skewness as a fx of waiting costs
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.5*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
set(gca,'FontSize',labfontsz,'FontName',fontname)
marg=[.1 .063];
xpos=.3;

legobj=zeros(4,1);
hold on
for i=1:4
    legobj(i)=plot(vec,skewvals(funstoplot(i),endvals),markertypes{i},'Color',fourcols(colorperm(i),:),'LineWidth',lw,'Markersize',5);
%     plot(vec,skewvals(funstoplot(i),endvals),'+','Color',fourcols(i,:),'LineWidth',lw,'Markersize',markersz)
end
xlabel('Decision time','FontName',fontname,'FontSize',textfontsz);
ylabel('Skewness','FontName',fontname,'FontSize',textfontsz);
set(gca,'FontSize',labfontsz,'FontName',fontname)
set(gca,'TickDir','out','ticklength',[.02 .025])
% set(gca,'Position',[.07 .25 .17 .7]);
[leg,hobj]=legend(legobj,funtypes{3},funtypes{1},funtypes{2},funtypes{4});
legend('boxoff');
set(leg,'Position',[.6 .7 .3 .3])
textobj = findobj(hobj, 'type', 'line');
for i=[1 3 5 7]
    set(textobj(i),'XData',[.15 .25])
end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','skewness_max','.pdf');
print(filename,'-dpdf','-r300');
