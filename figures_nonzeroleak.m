%% plotting parameters
axislw=1.25;
lw=3;
markersz=3;
textfontsz=12;
labfontsz=10;
fontname='Times New Roman';
yellow=[1 1 0];
red=[.9 0 0];
green=[0 .7 0];
blue=[0 0 .9];
% mycolormap=cbrewer('seq', 'YlOrRd', 100);
% mycolormap=cbrewer('seq','BuPu',100);
mycolormap=cbrewer('seq','Reds',64);
mycolormap=mycolormap(1:end,:);
% fourcols=[yellow;red;green;blue];
fourcols=cbrewer('qual','Set1',4);
fourcols=fourcols([3 1 2 4],:);
alphabet=['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h'];
powerfuns={'D','R','Delta','H','Pi''C','probs'};
fourpowerfuns={'D','R','H','C'};
funtypes={'number of signals','entropy','number of signalers','eigenvector centrality'};
cornerlabels={'Error rate','Personal preference','Decision time'};
markertypes={'-o','-d','-v','-s'};
funstoplot=[1 2 4 6];
colorperm=[1 2 4 3];
numpowerfuns=length(powerfuns);
load /Users/eleanorbrush/Dropbox/signaling_network/variables
load /Users/eleanorbrush/Dropbox/signaling_network/optanalysis_output.mat
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

%% find how nash equilibria depend on c1,c2,c3

nasheq_ind=cell(Nc,Nd);
giveup=zeros(Nc,Nd);
nashaccuracy=zeros(Nc,Nd);

for i=1:Nc
    c1=X(i);
    c2=Y(i);
    c3=Z(i);
    
    perf=zeros(2,Nd,Nt,Nt); %individual, leak, dominance, thresholds

    perf(1,2:Nd,:,:)=c1*(1-twomat2(leak,2:Nd,:,:,1))+c2*(twomat2(leak,2:Nd,:,:,2))+c3*(1-twomat2(leak,2:Nd,:,:,1));
    perf(2,2:Nd,:,:)=c1*(1-twomat2(leak,2:Nd,:,:,1))+c2*(twomat2(leak,2:Nd,:,:,2))+c3*(twomat2(leak,2:Nd,:,:,1));

    perf(1,1,:,:)=c2*(twomat2(leak,1,:,:,2))+c3*(1-twomat2(leak,1,:,:,1));
    perf(2,1,:,:)=c2*(twomat2(leak,1,:,:,2))+c3*(twomat2(leak,1,:,:,1));

        for k=2:Nd
            T1_best=zeros(1,Nt); 
            T2_best=zeros(1,Nt);
            for u=1:Nt
            [~,m]=min(perf(1,k,:,u));
            T1_best(u)=m;
            [~,m]=min(perf(2,k,u,:));
            T2_best(u)=m;
            end
            pp1=interp1(indices,reshape(T1_best,1,[]),'linear','pp');
            pp2=interp1(indices,reshape(T2_best,1,[]),'linear','pp');
            v1=ppval(pp2,indices);
            v2=ppval(pp1,v1);
            T1=indices(abs((indices)-v2)<0.01);
            T2=ppval(pp2,T1);
            nasheq_ind{i,k}=[T1;T2];
            nashaccuracy(i,k)=twomat2(leak,k,T1,T2,1);
        end
    
end

worstaccuracy=nashaccuracy(:,2);

%% find properties of group learning

expowervecs=zeros(N,7,Nc);
skewvals=zeros(7,Nc);
mutvals=zeros(7 ,Nc);
accuracy=zeros(1,Nc);
time=zeros(1,Nc);
avgthresh=zeros(N,Nc);
stdthresh=zeros(N,Nc);
timemats=zeros(N,N,Nc);
probmats=zeros(N,N,Nc);

bottomquart=16:20;
bottomaccuracy=zeros(1,Nc);

topquart=1:5;
topaccuracy=zeros(1,Nc);

for i=1:Nc
    c1=X(i);
    c2=Y(i);
    c3=Z(i);
    t=group_props_parallelized(c1,c2,c3,its,N,'unif',Nd,Nt,twomat2,domvals,threshvals,leak);
    expowervecs(:,:,i)=t{1};
    mutvals(:,i)=t{2};
    skewvals(:,i)=t{3};
    accuracy(i)=t{4};
    time(i)=t{5};
    avgthresh(:,i)=t{7};
    stdthresh(:,i)=t{8};
    timemats(:,:,i)=t{9};
    probmats(:,:,i)=t{10};
    totalsigmat=t{11};
    bottomaccuracy(i)=sum(sum(triu(totalsigmat(bottomquart,bottomquart),1)))/sum(sum(totalsigmat(bottomquart,bottomquart)));
    topaccuracy(i)=sum(sum(triu(totalsigmat(topquart,topquart),1)))/sum(sum(totalsigmat(topquart,topquart)));
end
%% example power distribution
exstofind=[layers{1}(1) layers{2}(1) layers{5}(1) layers{9}(1)];
Nex=length(exstofind);
expowervecs=zeros(N,7,Nex);
avgskew=zeros(7,Nex);

for i=1:Nex
    I=exstofind(i);
    c1=X(I);
    c2=Y(I);
    c3=Z(I);
    t=group_props_parallelized(c1,c2,c3,its,N,'unif',Nd,Nt,twomat2,domvals,threshvals,leak);
    expowervecs(:,:,i)=t{1};
    avgskew(:,i)=t{3};
end
%% save group props output
filename='/Users/eleanorbrush/Desktop/optanalysis_output.mat';
% save(filename,'nasheq_ind','nashaccuracy','worstaccuracy','expowervec','skewvals','mutvals','accuracy','time','bottomaccuracy','topaccuracy','avgthresh','stdthresh','timemats');
%% scale vectors for heatmaps

scalednashacc=nashaccuracy-repmat(min(nashaccuracy,[],1),Nc,1);
scalednashacc=scalednashacc./repmat(max(scalednashacc,[],1),Nc,1);

scaledworst=worstaccuracy-min(worstaccuracy);
scaledworst=scaledworst/max(scaledworst);
%%
scaledskewvals=skewvals-(repmat(min(skewvals,[],2),1,Nc));
scaledskewvals=scaledskewvals./repmat(max(scaledskewvals,[],2),1,Nc);

scaledmutvals=mutvals-(repmat(min(mutvals,[],2),1,Nc));
scaledmutvals=scaledmutvals./repmat(max(scaledmutvals,[],2),1,Nc);

[~,order]=sort(mutvals(1:6,:));
bestfun=order(6,:);
bestfun(bestfun==3)=2;
bestfun(bestfun==5)=2;
bestfun(bestfun==6)=3;
scaledbestfun=bestfun-min(min(bestfun));
scaledbestfun=scaledbestfun/max(max(scaledbestfun));

scaledaccuracy=accuracy-min(accuracy);
scaledaccuracy=scaledaccuracy/max(scaledaccuracy);

scaledbottomacc=bottomaccuracy-min(bottomaccuracy);
scaledbottomacc=scaledbottomacc/max(scaledbottomacc);

scaledtopacc=topaccuracy-min(topaccuracy);
scaledtopacc=scaledtopacc/max(scaledtopacc);

scaledtime=time-min(time);
scaledtime=scaledtime/max(scaledtime);
%% pairwise accuracy triangle heatmap figure
interestingi=[layers{1}(end) layers{3}(end) Nc 1 layers{3}(1)];
mainlabel='Accuracy';

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.75*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)
marg=[.1 .063];
xpos=.3;

colormap(mycolormap);
triimage(worstaccuracy,transformed,cornerlabels,mainlabel,'on',interestingi,2,{'b','c','d','e','f','g'},'k',.5,1)
annotation('textbox',[.3 .93 .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w*.55 0 w*2 h]);

filename=strcat('/Users/eleanorbrush/Desktop/','pairwise_accuracy','.pdf');
% print(filename,'-dpdf','-r300');


%% nasheq figure

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.35*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)
marg=[.14 .1];
xpos=.35;

toplot=ceil(length(interestingi)/2)*2;
for I=1:length(interestingi)
    i=interestingi(I);
    v=zeros(2,Nd);
    for j=2:Nd
        Neq=nasheq_ind{i,j};
        v(1,j)=threshvals(Neq(1));
        v(2,j)=threshvals(Neq(2));
    end
    subplot_tight(2,toplot/2,I,marg)
%     subplot(2,toplot/2,I)
    hold on
    plot(domvals(2:Nd),v(1,2:Nd),'-','Color','blue','LineWidth',lw)
    plot(domvals(2:Nd),v(2,2:Nd),'--r','LineWidth',lw)
    set(gca,'YLim',[0 max(threshvals)+threshvals(1)]);
    set(gca,'FontSize',labfontsz,'FontName',fontname)
%     titlelab=['w2=',num2str(cvals(i))];
%     title(titlelab,'FontSize',textfontsz);
    switch I
        case 1
            ylab=ylabel('Threshold','FontSize',textfontsz,'FontName',fontname);
            v=get(ylab,'Position');
            set(get(gca,'ylabel'),'Position',[xpos v(2) v(3)])
            leg=legend('Stronger','Weaker' );
            legend boxoff
            set(leg,'Position',[.13 .67 .1 .1],'FontName',fontname)
            legendshrink(.5,'best',leg)
%         case 1
            xlab=xlabel('c','FontSize',textfontsz);
            v=get(xlab,'Position');
            set(xlab,'Position',v+[-.07 .5 0]);
    end
    switch i
        case Nc
            v=get(gca,'Position');
            set(gca,'Position',[v(1) .4 v(3) v(4)]);
    end
    box on
    ylim=get(gca,'Ylim');
    [x,y]=data2norm(.6,ylim(2));
    switch ismember(I,1:3)
        case 1
            annotation('textbox',[x-.1 y+.12 .05 .04],'String',['(' alphabet(I+1) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
        case 0
            annotation('textbox',[x-.1 y+.04 .05 .04],'String',['(' alphabet(I+1) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
    end
end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','nasheq_thresholds','.pdf');
% print(filename,'-dpdf','-r300');


%% group accuracy triangle heatmap figure
interestingi=[11 21 1 12];
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.5*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)
marg=[.1 .05];
xpos=.3;

vals=[accuracy;bottomaccuracy];
% mainlabel='Average accuracy';
colormap(mycolormap);
triimage_double(vals,transformed,cornerlabels,'Accuracy','on',interestingi,2,{'c','d','e','f','g','h'},'k',.5,1)

annotation('textbox',[.05 .93 .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
annotation('textbox',[.45 .93 .01 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','groupaccuracy','.pdf');
% print(filename,'-dpdf','-r300');

%% group eq thresholds figure
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.35*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)
marg=[.14 .063];
xpos=-3;


lineplot=[1 2 5 6];
timeplot=[3 4 7 8];

toplot=length(interestingi);
for I=1:toplot
    i=interestingi(I);
    subplot_tight(2,toplot,lineplot(I),marg)
    hold on
    plot(1:N,avgthresh(:,i),'-','Color','k','LineWidth',lw)
%     plot(1:N,avgthresh(:,i)-stdthresh(:,i),'-','Color','r','LineWidth',lw)
%     plot(1:N,avgthresh(:,i)+stdthresh(:,i),'-','Color','r','LineWidth',lw)
    set(gca,'YLim',[0 max(threshvals)+threshvals(1)]);
    set(gca,'FontSize',labfontsz,'FontName',fontname)
%     titlelab=['w2=',num2str(cvals(i))];
%     title(titlelab,'FontSize',textfontsz);
    switch I
        case 1
            ylabel('Threshold','FontSize',textfontsz,'FontName',fontname)
            v=get(get(gca,'ylabel'),'Position');
            set(get(gca,'ylabel'),'Position',[xpos v(2) v(3)])
%         case 4
            xlab=xlabel('Order of ability','FontSize',textfontsz);
            v=get(xlab,'Position');
            set(xlab,'Position',[v(1)+.3 -.55 v(3)])      
    end
    box on
    switch ismember(I,[1 2])
        case 1
            ylim=get(gca,'Ylim');
            [x,y]=data2norm(0,ylim(2));
            annotation('textbox',[x-.04 .95 .05 .04],'String',['(' alphabet(I) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
        case 0
            ylim=get(gca,'Ylim');
            [x,y]=data2norm(0,ylim(2));
            annotation('textbox',[x-.04 .5 .05 .04],'String',['(' alphabet(I) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
    end
end

maxtime=max(max(max(timemats(:,:,interestingi))));

for I=1:toplot
    i=interestingi(I);
    
    subplot_tight(2,toplot,timeplot(I),marg)
    imagesc(reshape(timemats(:,:,i),N,N))
    colormap(mycolormap);
    caxis manual
    caxis([ 0 maxtime])
    set(gca,'FontSize',labfontsz,'FontName',fontname)
    switch I
        case toplot
            c=colorbar;
            v=get(c,'position');
            set(c,'position',[.954 v(2)+.015 v(3)-.005 v(4)-.01],'FontSize',labfontsz,'FontName',fontname);
            yl=get(c,'YLabel');
            v=get(yl,'Position');
            set(yl,'String','Time','FontSize',labfontsz,'Rotation',270,'Position',v+[2 0 0],'FontName',fontname)
    end
    switch ismember(I,[1 2])
        case 1
            ylim=get(gca,'Ylim');
            [x,y]=data2norm(0,ylim(2));
            annotation('textbox',[x-.04 .95 .05 .04],'String',['(' alphabet(I+2) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
        case 0
            ylim=get(gca,'Ylim');
            [x,y]=data2norm(0,ylim(2));
            annotation('textbox',[x-.04 .5 .05 .04],'String',['(' alphabet(I+2) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
    end
end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','groupeq_thresholds','.pdf');
% print(filename,'-dpdf','-r300');


%% mutual information triangle heatmap figure
numfun=6;
mainlabel='Mutual information';

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.75*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)
marg=[.1 .063];
xpos=.3;

colormap(mycolormap);
triimage(mutvals(numfun,:),transformed,cornerlabels,mainlabel,'off')

annotation('textbox',[.3 .93 .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w*.55 0 w*2 h]);

filename=strcat('/Users/eleanorbrush/Desktop/','mutinfo','.pdf');
print(filename,'-dpdf','-r300');

%% skewness triangle heatmap figure
numfun=1;
mainlabel='Skewness';

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.75*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)
marg=[.1 .063];
xpos=.3;

colormap(mycolormap);
triimage(skewvals(numfun,:),transformed,cornerlabels,mainlabel,'on',exstofind,2,{'b','c','d'},'k')

% annotation('textbox',[.3 .93 .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w*.55 0 w*2 h]);

filename=strcat('/Users/eleanorbrush/Desktop/','skewness','.pdf');
print(filename,'-dpdf','-r300');

%% skewness histograms

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.35*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)
marg=[.14 .1];
xpos=-2;
numfun=1;

toplot=[1 3 4];
for I=1:3
    i=toplot(I);
    subplot_tight(1,3,I,marg)
%     hold on
    hist(expowervecs(:,numfun,i));
    histobj = findobj(gca,'Type','patch');
    set(histobj,'FaceColor',fourcols(1,:));
    set(gca,'FontName',fontname,'FontSize',labfontsz)
    v=get(gca,'Position');
    set(gca,'Position',[v(1) .2 v(3) v(4)]);
    xlim=get(gca,'xlim');
    ylim=get(gca,'Ylim');
    switch I
        case 1
            ylab=ylabel('Frequency','FontSize',textfontsz,'FontName',fontname);
            v=get(ylab,'Position');
%             set(get(gca,'ylabel'),'Position',[xpos v(2) v(3)])
%         case 2
            xlab=xlabel('Number of signalers','FontSize',textfontsz,'FontName',fontname);
            v=get(xlab,'Position');
            set(xlab,'Position',[v(1) ylim(1)-diff(ylim)*.12 v(3)]);
    end
    box off
    [x,y]=data2norm(xlim(1),ylim(2));
    [x2,y2]=data2norm(xlim(2),ylim(1));
    annotation('textbox',[x-.3*(x2-x) y+.01 .05 .04],'String',['(' alphabet(I+1) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','skewness_histograms','.pdf');
print(filename,'-dpdf','-r300');
%% skewness of four functions as a fx of waiting costs
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.5*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)
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

filename=strcat('/Users/eleanorbrush/Desktop/','skewness_details','.pdf');
print(filename,'-dpdf','-r300');
%% best fun triangle heatmap figure
yellow=[1 1 0];
red=[.9 0 0];
green=[0 .7 0];
blue=[0 0 .9];

mainlabel='Most informative';
funtypes2=funtypes;
funtypes2{4}='eigenvector';
funtypes2{5}='centrality';
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=2/3*6.83;
h=.75*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)
marg=[.1 .063];
xpos=.3;

triimage_discrete(scaledbestfun,bestfun,transformed,cornerlabels,funtypes2,fourcols)

% annotation('textbox',[.3 .93 .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w*.55 0 w*2 h]);

filename=strcat('/Users/eleanorbrush/Desktop/','mostinformative','.pdf');
print(filename,'-dpdf','-r300');