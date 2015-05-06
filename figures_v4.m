%% plotting parameters
axislw=1.25;
lw=3;
trilw=1;
doubletrilw=1;
markersz=3;
textfontsz=8;
labfontsz=6;
fontname='Helvetica';
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
fourcols_light=brighten(fourcols,.75);
alphabet=['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h'];
romans={'i', 'ii', 'iii', 'iv'}';
powerfuns={'D','R','Delta','H','Pi''C','probs'};
fourpowerfuns={'D','R','H','C'};
funtypes={'number of signals','entropy','number of signalers','eigenvector centrality'};
cornerlabels={'Error rate','Decision preference','Decision time'};
cornerlabels_short={'ER','DP','DT'};
markertypes={'-o','-d','-v','-s'};
funstoplot=[1 2 4 6];
colorperm=[1 2 4 3];
numpowerfuns=length(powerfuns);
load /Users/eleanorbrush/Dropbox/signaling_network/twomat.mat
load /Users/eleanorbrush/Dropbox/signaling_network/variables
load /Users/eleanorbrush/Dropbox/signaling_network/optanalysis_output.mat
load /Users/eleanorbrush/Dropbox/signaling_network/infogrowth_output.mat

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

%% find how nash equilibria depend on c1,c2,c3

nasheq_ind=cell(Nc,Nd);
giveup=zeros(Nc,Nd);
nashaccuracy=zeros(Nc,Nd);

for i=1:Nc
    c1=X(i);
    c2=Y(i);
    c3=Z(i);
    
    perf=zeros(2,Nd,Nt,Nt); %individual, leak, dominance, thresholds

    perf(1,2:Nd,:,:)=c1*(1-twomat(leak,2:Nd,:,:,1))+c2*(twomat(leak,2:Nd,:,:,2))+c3*(1-twomat(leak,2:Nd,:,:,1));
    perf(2,2:Nd,:,:)=c1*(1-twomat(leak,2:Nd,:,:,1))+c2*(twomat(leak,2:Nd,:,:,2))+c3*(twomat(leak,2:Nd,:,:,1));

    perf(1,1,:,:)=c2*(twomat(leak,1,:,:,2))+c3*(1-twomat(leak,1,:,:,1));
    perf(2,1,:,:)=c2*(twomat(leak,1,:,:,2))+c3*(twomat(leak,1,:,:,1));

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
            nashaccuracy(i,k)=twomat(leak,k,T1,T2,1);
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
    t=group_props_parallelized(c1,c2,c3,its,N,'unif',Nd,Nt,twomat,domvals,threshvals,leak);
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
exstofind=[layers{1}(1) layers{5}(1) layers{7}(1)];
Nex=length(exstofind);
expowervecs=zeros(N,7,Nex);
avgskew=zeros(7,Nex);

for i=1:Nex
    I=exstofind(i);
    c1=X(I);
    c2=Y(I);
    c3=Z(I);
    t=group_props_parallelized(c1,c2,c3,its,N,'unif',Nd,Nt,twomat,domvals,threshvals,leak);
    expowervecs(:,:,i)=t{1};
    avgskew(:,i)=t{3};
end
%% info over time
layerstoplot=[1 3 5];
Nlays=length(layerstoplot);
Ntime=50;
mut_time=zeros(Nlays,6,Ntime+1);
for I=1:Nlays
    i=layers{layerstoplot(I)}(1);
    c1=X(i);
    c2=Y(i);
    c3=Z(i);
    reshapedtime=reshape(timemats(:,:,i),N,N);
    tmin=min(min(reshapedtime(reshapedtime~=0)));
    tmax=max(max(reshapedtime));
    step=(tmax-tmin)/(Ntime-1);
    timevec=(tmin-step):step:tmax;
    t=group_props_time(c1,c2,c3,its,N,'unif',Nd,Nt,twomat,domvals,threshvals,leak,timevec);
    mut_time(I,:,:)=t;
end
%% save outputs
filename='/Users/eleanorbrush/Desktop/optanalysis_output.mat';
% save(filename,'nasheq_ind','nashaccuracy','worstaccuracy','expowervec','skewvals','mutvals','accuracy','time','bottomaccuracy','topaccuracy','avgthresh','stdthresh','timemats');

filename='/Users/eleanorbrush/Desktop/infogrowth_output.mat';
% save(filename,'mut_time','layerstoplot','Nlays','Ntime');

filename='/Users/eleanorbrush/Dropbox/signaling_network/expowervecs.mat';
save(filename,'expowervecs','avgskew','exstofind');

%% scale vectors for heatmaps

scalednashacc=nashaccuracy-repmat(min(nashaccuracy,[],1),Nc,1);
scalednashacc=scalednashacc./repmat(max(scalednashacc,[],1),Nc,1);

scaledworst=worstaccuracy-min(worstaccuracy);
scaledworst=scaledworst/max(scaledworst);

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
%% pairwise nash eq figure
interestingi=[layers{1}(end) 1 layers{3}(1) Nc ];
mainlabel='Accuracy';

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.25*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
set(gca,'FontSize',labfontsz)
xpos=.43;
X=.22;
Y=.5;
height=.33;
width=.18;
xmarg=.1;
ymarg=.12;
width=.17;
height=.3;

boxes=[1 2 5 6];
toplot=ceil(length(interestingi)/2)*2;
for I=1:length(interestingi)
% for I=1:2
    i=interestingi(I);
    v=zeros(2,Nd);
    for j=2:Nd
        Neq=nasheq_ind{i,j};
        v(1,j)=threshvals(Neq(1));
        v(2,j)=threshvals(Neq(2));
    end
    subplot(2,4,boxes(I))
    hold on
    plot(domvals(2:Nd),v(1,2:Nd),'-','Color',fourcols(3,:),'LineWidth',lw)
    plot(domvals(2:Nd),v(2,2:Nd),'--','Color',fourcols(2,:),'LineWidth',lw)
    set(gca,'YLim',[0 max(threshvals)+threshvals(1)]);
    set(gca,'FontSize',labfontsz,'FontName',fontname)
    switch I
        case 1
            ylab=ylabel('Threshold','FontSize',textfontsz,'FontName',fontname);
            v=get(ylab,'Position');
            set(get(gca,'ylabel'),'Position',[xpos v(2)-.1 v(3)])
            leg=legend('Stronger','Weaker' );
            legend boxoff
            set(leg,'Position',[.49 .4 .1 .1],'FontName',fontname)
            legendshrink(.5,'best',leg)
            xlab=xlabel('Difficulty, c','FontSize',textfontsz);
            v=get(xlab,'Position');
            set(xlab,'Position',[.75 v(2) v(3)]);
    end
    x=1-mod(I,2);
    y=(sign(2.5-I)+1)/2;
    set(gca,'Position',[xmarg+x*X ymarg+y*Y width height]);
%     box on
    set(gca,'tickdir','out')
    set(gca,'xlim',[.5 1])
    ylim=get(gca,'Ylim');
    [x,y]=data2norm(.6,ylim(2));
    annotation('textbox',[x-.07 y+.07 .05 .04],'String',['(' alphabet(I) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
end

subplot(2,4,[3 4 7 8])
colormap(mycolormap);
xoffset=.25;
yoffset=.2;
[mybar,yl]=triimage_returnhandles(worstaccuracy,transformed,cornerlabels,mainlabel,'on',interestingi,trilw,{'a','b','c','d'},'k',.5,1,1,xoffset,yoffset);

v=get(gca,'Position');
set(gca,'Position',[.55 .1 .4 .8])
set(mybar,'Position',[.9 .1 .02 .8])
set(yl,'Position',[9 .4792 1.0001])

annotation('textbox',[.6 .93 .01 .04],'String','(e)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w*.05 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','Figure1','.pdf');
print(filename,'-dpdf','-r600');


filename=strcat('/Users/eleanorbrush/Desktop/test.pdf');
print(filename,'-dpdf','-r600');



%% group nash eq figure
interestingi=[11 1 21 12];

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
w=6.83;
h=.25*w;
set(gcf,'Units','inches');
set(gcf,'Position',[v(1) v(2) w h]);
set(gca,'FontSize',labfontsz)

xpos=-3.5;
yanno=.985;

lw=1;
trilw=2;
markersz=2;

AxesH=subplot(1,4,1);
hold on
toplot=length(interestingi);
p=zeros(1,toplot);
for I=1:toplot
    i=interestingi(I);
    fill([1:N N:-1:1],[avgthresh(:,i)-stdthresh(:,i); (avgthresh(end:-1:1,i)+stdthresh(end:-1:1,i))],fourcols_light(I,:),'EdgeColor','none')
    p(I)=plot(1:N,avgthresh(:,i),'-','Color',fourcols(I,:),'LineWidth',lw);
    plot(1:N,avgthresh(:,i),'o','Color',fourcols(I,:),'LineWidth',lw,'MarkerSize',markersz,'MarkerFaceColor',fourcols(I,:));    
end

set(gca,'YLim',[.4 2.1]);
set(gca,'FontSize',labfontsz,'FontName',fontname)
set(gca,'tickdir','out','xtick',[1 10 20])
ylabel('Threshold','FontSize',textfontsz,'FontName',fontname)
v=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',[xpos v(2) v(3)])
xlab=xlabel('Order of ability','FontSize',textfontsz);
v=get(xlab,'Position');
set(xlab,'Position',[v(1)+.3 v(2) v(3)])      

box off
%     axis on

[leg,hobj]=legend(p,romans{1:4});
legend boxoff
set(leg,'Position',[.02 .25 .09 .33])
textobj=findobj(hobj,'type','line');
for i=[1 3 5 7]
    set(textobj(i),'xdata',[.65 .73])
end

v=get(gca,'Position');
set(gca,'Position',[.06 .17 v(3) .65])

xlim=get(gca,'xlim');
ylim=get(gca,'Ylim');
[x,y]=data2norm(-1,ylim(2));
annotation('textbox',[x-.05 yanno .05 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')


subplot(1,4,2)
colormap(mycolormap);
xoffset=.05;
yoffset=.25;
triimage_returnhandles(accuracy,transformed,cornerlabels_short,mainlabel,'on',interestingi,trilw,{'i','ii','iii','iv'},fourcols,.5,1,0,xoffset,yoffset);
set(gca,'Position',[.24 .02 .2 .95])
annotation('textbox',[.24 yanno .05 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(1,4,3)

colormap(mycolormap);
[mybar,yl]=triimage_returnhandles(bottomaccuracy,transformed,cornerlabels_short,mainlabel,'on',interestingi,trilw,{'i','ii','iii','iv'},fourcols,.5,1,1,xoffset,yoffset);
set(gca,'Position',[.43 .02 .2 .95])

set(mybar,'Position',[.64 .19 .01 .65])
v=get(yl,'Position');
set(yl,'Position',[14 v(2:3)])

annotation('textbox',[.43 yanno .05 .04],'String','(c)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

maxtime=max(max(max(timemats(:,:,interestingi))));
xstart=.75;
width=.08;
height=.32;
timeplot=[7 8 15 16];
for I=[2 1 4 3]
    i=interestingi(I);
    
    subplot(2,8,timeplot(I))
    imagesc(reshape(timemats(:,:,i),N,N))
    colormap(mycolormap);
    caxis manual
    caxis([ 0 maxtime])
    set(gca,'FontSize',labfontsz,'FontName',fontname)
    set(gca,'ydir','normal')
    set(gca,'xtick',[1 10 20])
    set(gca,'ytick',[1 10 20])
    switch ismember(I,[2])
        case 1
            c=colorbar;
            v=get(c,'position');
            set(c,'position',[.97 v(2)+.015 v(3)-.005 v(4)-.035],'FontSize',labfontsz,'FontName',fontname);
            yl=get(c,'YLabel');
            v=get(yl,'Position');
            set(yl,'String','Time','FontSize',labfontsz,'Rotation',270,'Position',v+[2 0 0],'FontName',fontname)
    end
    mytitle=title(romans{I});
    v=get(mytitle,'Position');
    set(mytitle,'Position',[v(1) v(2)-.5 v(3)])
    switch ismember(I,[1 2])
        case 1
            ylim=get(gca,'Ylim');
            [x,y]=data2norm(0.5,ylim(2));
            annotation('textbox',[x-mod(I,2)*.03 yanno .05 .04],'String',['(' alphabet(I+3) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
        case 0
            ylim=get(gca,'Ylim');
            [x,y]=data2norm(0.5,ylim(2));
            annotation('textbox',[x-mod(I,2)*.03 .53 .05 .04],'String',['(' alphabet(I+3) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
    end
    v=get(gca,'Position');
    x=1-mod(I,2);
    set(gca,'Position',[xstart+x*(width+.04) v(2) width height])
end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','Figure2','.pdf');
print(filename,'-dpdf','-r300');



filename=strcat('/Users/eleanorbrush/Desktop/test.pdf');
print(filename,'-dpdf','-r800');



%% skewness triangle heatmap figure
numfun=1;
mainlabel='Skewness';

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=3.4;
h=.75*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)
marg=[.1 .063];
xpos=.3;
exstofind=[layers{1}(1) layers{5}(1) layers{7}(1)];
colormap(mycolormap);
m=min(skewvals(numfun,:));
M=max(skewvals(numfun,:));
xoffset=.25;
yoffset=.2;
[mybar,yl]=triimage_returnhandles(skewvals(numfun,:),transformed,cornerlabels,mainlabel,'off',exstofind,trilw,{'i','ii','iii','iv'},fourcols,m,M,1,xoffset,yoffset);
% annotation('textbox',[.3 .93 .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

v=get(gca,'Position');
set(gca,'Position',[.1 v(2) .7 .8])
set(mybar,'Position',[.65 .13 .02 .77])
v=get(yl,'Position');
set(yl,'Position',[v(1)+2 v(2) v(3)])

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w*.55 0 w*2 h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','Figure4','.pdf');
% print(filename,'-dpdf','-r300');

filename=strcat('/Users/eleanorbrush/Desktop/test.pdf');
% print(filename,'-dpdf','-r300');


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
%% combined best fun triangle and time plots
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
w=3.4;
h=.75*w;
set(gcf,'Units','inches');
set(gcf,'Position',[v(1) v(2) w h]);
% set(gca,'FontSize',labfontsz)
marg=[.1 .063];
xpos=.3;

% lw=
trilw=1;

layerstoplot=[1 3 5];
Nlays=length(layerstoplot);

highlight='on';
tohighlight=[layers{layerstoplot(1)}(1) layers{layerstoplot(2)}(1) layers{layerstoplot(3)}(1)];
pointlabels={'b','c','d'};

subplot(3,4,[1:3 5:7 9:11])

triimage_discrete2(scaledbestfun,bestfun,transformed,cornerlabels,funtypes2,fourcols,highlight,tohighlight,trilw,pointlabels)
set(gca,'xlim',[0-2/11 2+2/11],'ylim',[0-1.7321/11 1.7321+1.7321/11])
v=get(gca,'Position');
ratio=v(4)/v(3);
scale=1;
set(gca,'Position',[.05 v(2)-.04 scale*v(3) scale*v(4)])


% linesty={'-.','--','-'};
linesty={'-','-','-'};
p=zeros(1,7);
for now=1:Nlays

    subplot(3,3,3*(now-1)+3)
    hold on
    j=layerstoplot(now);
    k=layers{j}(1);
    reshapedtime=reshape(timemats(:,:,k),N,N);
    tmin=min(min(reshapedtime(reshapedtime~=0)));
    tmax=max(max(reshapedtime));
    step=(tmax-tmin)/(Ntime-1);
    timevec=(tmin-step):step:tmax;
    
    fourcolsnow=rgb2hsv(fourcols);
%     fourcolsnow(:,2)=fourcolsnow(:,2)*(1-.4*(now-1));
    fourcolsnow=hsv2rgb(fourcolsnow);
    for i=1:4
    I=funstoplot(i);
    p(i+3)=plot(log(timevec),[0; (col(mut_time(now,I,2:end)))],'-','Color',fourcolsnow(colorperm(i),:),'LineStyle',linesty{now},'LineWidth',lw);
    end
    [~,which]=max(mut_time(now,funstoplot,end));
    plot(log(timevec),[0; (col(mut_time(now,funstoplot(which),2:end)))],'-','Color',fourcols(colorperm(which),:),'LineStyle',linesty{now},'LineWidth',lw)
    set(gca,'xlim',[log(sigfig(tmin)-3*step) log(sigfig(tmax)+step)])
%     v=get(gca,'Position');
%     set(gca,'Position',[.2 v(2:end)]);
    set(gca,'FontName',fontname,'FontSize',labfontsz)
%     timeticks=[.25 .5 1:.5:5];
    switch now
        case 3
        timeticks=[.25 .5 1  2];
        case 2
        timeticks=[.25 1  2];
        case 1
        timeticks=[2:4];
    end
    set(gca,'xtick',log(timeticks),'xticklabel',(timeticks))
    v=get(gca,'Position');
    set(gca,'Position',[.72 v(2:4)])
    ylim=get(gca,'Ylim');
    xlim=get(gca,'xlim');
    [x,y]=data2norm(xlim(1),ylim(2));
    annotation('textbox',[x-.1 y+.04 .05 .04],'String',['(' alphabet(now+1) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
end

xlabel('Time','FontName',fontname,'FontSize',labfontsz)
ylab=ylabel('Mutual information','FontName',fontname,'FontSize',labfontsz);
v=get(ylab,'Position');
set(ylab,'Position',[log(timeticks(1))-.85 v(2)-.12 v(3)]);
[leg,hobj]=legend(p([7 5 4 6]),{'eigenvector centrality','weighted degree','unweighted degree','entropy'});

legend('boxoff')
set(leg,'Position',[.03 .63 .1 .1],'FontName',fontname)
textobj = findobj(hobj, 'type', 'line');
for i=[1 3 5 7]
%     set(textobj(i),'XData',[.25 .29])
end

annotation('textbox',[.05 .95 .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','Figure3','.pdf');
print(filename,'-dpdf','-r300');

filename=strcat('/Users/eleanorbrush/Desktop/test.pdf');
print(filename,'-dpdf','-r300');