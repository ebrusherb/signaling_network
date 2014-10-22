%% plotting parameters
axislw=1.25;
lw=3;
textfontsz=12;
labfontsz=10;
fontname='Times New Roman';
yellow=[1 1 0];
red=[.9 0 0];
green=[0 .7 0];
blue=[0 0 .9];
alphabet=['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h'];
powerfuns={'D','R','Delta','H','Pi''C','probs'};
numpowerfuns=length(powerfuns);
%% coordinate set up
N=20;
its=1000;
cstep=.1;
vec=0:cstep:1;
Nc=length(vec);

% X=repmat(0:cstep:1,1,11);
% Y=reshape(repmat(0:cstep:1,11,1),1,[]);

X=[];
Y=[];
layers=cell(Nc,1);
for i=0:(Nc-1)
    layers{i+1}=length(X)+(1:(Nc-i));
    X=[X, vec(1:(Nc-i))];
    Y=[Y,vec(i+1)*ones(1,Nc-i)];   
end
Z=1-X-Y;

transformed=transform([X',Z',Y']);
Nc=length(X);

%% find how nash equilibria depend on c1,c2,c3
indices=1:Nt;
leak=1;

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
    t=group_props_parallelized(c1,c2,c3,its,N,'unif',Nd,Nt,twomat,domvals,threshvals);
    expowervecs(:,:,i)=t{1};
    mutvals(:,i)=t{2};
    skewvals(:,i)=t{3};
    accuracy(i)=t{4};
    time(i)=t{5};
    avgthresh(:,i)=t{7};
    timemats(:,:,i)=t{8};
    probmats(:,:,i)=t{9};
    totalsigmat=t{10};
    bottomaccuracy(i)=sum(sum(triu(totalsigmat(bottomquart,bottomquart),1)))/sum(sum(totalsigmat(bottomquart,bottomquart)));
    topaccuracy(i)=sum(sum(triu(totalsigmat(topquart,topquart),1)))/sum(sum(totalsigmat(topquart,topquart)));
end
%% save group props output
filename='triangle_heatmap_output.mat';
save(filename,'skewvals','mutvals','accuracy','time','avgthresh');
%% scale vectors for heatmaps
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

scalednashacc=nashaccuracy-repmat(min(nashaccuracy,[],1),Nc,1);
scalednashacc=scalednashacc./repmat(max(scalednashacc,[],1),Nc,1);

scaledworst=worstaccuracy-min(worstaccuracy);
scaledworst=scaledworst/max(scaledworst);
%% pairwise accuracy triangle heatmap figure
interestingi=[11 1 21 12 39 66];
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
cornerlabels={'ER','SP','DT'};

triimage(scaledworst,worstaccuracy,transformed,cornerlabels,mainlabel,'on',interestingi,2,{'b','c','d','e','f','g'},yellow)
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
marg=[.14 .063];
xpos=.35;

toplot=length(interestingi);
for I=1:toplot
    i=interestingi(I);
    v=zeros(2,Nd);
    for j=2:Nd
        Neq=nasheq_ind{i,j};
        v(1,j)=threshvals(Neq(1));
        v(2,j)=threshvals(Neq(2));
    end
    subplot_tight(2,toplot,I,marg)
    hold on
    plot(domvals(2:Nd),v(1,2:Nd),'Color','blue','LineWidth',lw)
    plot(domvals(2:Nd),v(2,2:Nd),'--r','LineWidth',lw)
    set(gca,'YLim',[0 max(threshvals)+threshvals(1)]);
    set(gca,'FontSize',labfontsz,'FontName',fontname)
%     titlelab=['w2=',num2str(cvals(i))];
%     title(titlelab,'FontSize',textfontsz);
    switch I
        case 1
            ylabel('Threshold','FontSize',textfontsz,'FontName',fontname)
            v=get(get(gca,'ylabel'),'Position');
            set(get(gca,'ylabel'),'Position',[xpos v(2) v(3)])
            leg=legend('Stronger','Weaker');
            legend boxoff
            set(leg,'Position',[.08 .68 .1 .1],'FontName',fontname)
            legendshrink(.5,'best',leg)
%         case 4
%             xlabel('Probability stronger animal wins','FontSize',textfontsz)
    end
    box on
    ylim=get(gca,'Ylim');
    [x,y]=data2norm(.6,ylim(2));
    annotation('textbox',[x-.05 y+.1 .05 .04],'String',['(' alphabet(I+1) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
end

for I=1:toplot
    i=interestingi(I);
    vprob=zeros(1,Nd);
    vtime=zeros(1,Nd);
    for j=2:Nd
        Neq=nasheq_ind{i,j};
        vprob(j)=twomat(leak,j,Neq(1),Neq(2),1);
        vtime(j)=twomat(leak,j,Neq(1),Neq(2),2);
    end
    
    subplot_tight(2,toplot,I+toplot,marg)
    plot(domvals(2:Nd),vtime(1,2:Nd),'k','LineWidth',lw)
%     set(gca,'ylim',[0 4])
    set(gca,'FontSize',labfontsz,'FontName',fontname)
    switch I
        case 1
            ylabel('Time','FontSize',textfontsz,'FontName',fontname)
            v=get(get(gca,'ylabel'),'Position');
            set(get(gca,'ylabel'),'Position',[xpos v(2) v(3)])
        case 3
            xlab=xlabel('c','FontSize',textfontsz,'FontName',fontname);
            v=get(xlab,'Position');
            set(xlab,'Position',[v(1) .226 v(3)]);
    end

end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','nasheq_thresholds','.pdf');
% print(filename,'-dpdf','-r300');


%% group accuracy triangle heatmap figure

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
cornerlabels={'ER','SP','DT'};

% mainlabel='Average accuracy';
% triimage(scaledaccuracy,acc3uracy,transformed,cornerlabels,mainlabel,'off')
triimage_double(scaledvals,vals,transformed,cornerlabels,'Accuracy','on',interestingi,2,{'c','d','e','f','g','h'},yellow)

annotation('textbox',[.05 .93 .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
annotation('textbox',[.45 .93 .01 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','groupaccuracy_heatmap','.pdf');
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

toplot=length(interestingi);
for I=1:toplot
    i=interestingi(I);
    subplot_tight(2,toplot,I,marg)
    hold on
    plot(1:N,avgthresh(:,i),'Color','k','LineWidth',lw)
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
%             xlabel('Probability stronger animal wins','FontSize',textfontsz)
    end
    box on
    ylim=get(gca,'Ylim');
    [x,y]=data2norm(0,ylim(2));
    annotation('textbox',[x-.05 .95 .05 .04],'String',['(' alphabet(I+2) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
end

maxtime=max(max(max(timemats(:,:,interestingi))));

for I=1:toplot
    i=interestingi(I);
    
    subplot_tight(2,toplot,I+toplot,marg)
%     v=get(axes,'position');
    imagesc(reshape(timemats(:,:,i),N,N))
%     set(axes,'position',[.05 v(2:end)])
    caxis manual
    caxis([ 0 maxtime])
    set(gca,'FontSize',labfontsz,'FontName',fontname)
    switch I
%         case 1
%             ylabel('Time','FontSize',textfontsz,'FontName',fontname)
%             v=get(get(gca,'ylabel'),'Position');
%             set(get(gca,'ylabel'),'Position',[xpos v(2) v(3)])
        case 3
            xlab=xlabel('True order of ability','FontSize',textfontsz,'FontName',fontname);
            v=get(xlab,'Position');
            set(xlab,'Position',v+[0 -2.5 0])
        case toplot
            c=colorbar;
            v=get(c,'position');
            set(c,'position',[.948 v(2)+.015 v(3) v(4)-.01],'FontSize',labfontsz);
            yl=get(c,'YLabel');
            v=get(yl,'Position');
            set(yl,'String','Time','FontSize',textfontsz,'Rotation',270,'Position',v+[3 0 0],'FontName',fontname)
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
cornerlabels={'ER','SP','DT'};

triimage(scaledmutvals(numfun,:),mutvals(numfun,:),transformed,cornerlabels,mainlabel,'off')

annotation('textbox',[.3 .93 .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w*.55 0 w*2 h]);

filename=strcat('/Users/eleanorbrush/Desktop/','mutinfo_heatmap','.pdf');
% print(filename,'-dpdf','-r300');

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
cornerlabels={'ER','SP','DT'};

triimage(scaledskewvals(numfun,:),skewvals(numfun,:),transformed,cornerlabels,mainlabel,'off')

% annotation('textbox',[.3 .93 .01 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w*.55 0 w*2 h]);

filename=strcat('/Users/eleanorbrush/Desktop/','skewness_heatmap','.pdf');
print(filename,'-dpdf','-r300');

%% best fun triangle heatmap figure
yellow=[1 1 0];
red=[.9 0 0];
green=[0 .7 0];
blue=[0 0 .9];

mainlabel='Most informative';

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
cornerlabels={'ER','SP','DT'};

funtypes={'number of signals','entropy','number of signalers','eigenvector centrality'};
fourcols=[yellow;red;green;blue];
triimage_discrete(scaledbestfun,bestfun,transformed,cornerlabels,funtypes,fourcols)

% annotation('textbox',[.3 .93 .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w*.55 0 w*2 h]);

filename=strcat('/Users/eleanorbrush/Desktop/','mostinformative_heatmap','.pdf');
print(filename,'-dpdf','-r300');