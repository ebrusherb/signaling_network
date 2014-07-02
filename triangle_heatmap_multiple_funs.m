%%
axislw=1.25;
lw=3;
labfontsz=10 ;
textfontsz=12;
alpha=.7;
ticklength=[.02,.0250];
markersz=5;
markercol='white';
fontname='Times New Roman';

%%
its=1000;
N=20;

cstep=.1;
vec=0:cstep:1;
l=length(vec);

% X=repmat(0:cstep:1,1,11);
% Y=reshape(repmat(0:cstep:1,11,1),1,[]);

X=[];
Y=[];
for i=0:(l-1)
    X=[X, vec(1:(end-i))];
    Y=[Y,vec(i+1)*ones(1,l-i)];
end
Z=1-X-Y;

transformed=transform([X',Z',Y']);

%%
l=length(X);
skewvals=zeros(6,l);
mutvals=zeros(6 ,l);
accuracy=zeros(1,l);

for i=1:l
    c1=X(i);
    c2=Y(i);
    c3=Z(i);
    t=group_props_parallelized(c1,c2,c3,its,N,'unif',Nd,Nt,twomat,domvals,threshvals);
    mutvals(:,i)=t{1};
    skewvals(:,i)=t{2};
    accuracy(i)=t{3};
end
%%
scaledskewvals=skewvals-(repmat(min(skewvals')',1,l));
scaledskewvals=scaledskewvals./repmat(max(scaledskewvals')',1,l);

scaledmutvals=mutvals-(repmat(min(mutvals')',1,l));
scaledmutvals=scaledmutvals./repmat(max(scaledmutvals')',1,l);

scaledaccuracy=accuracy-min(accuracy);
scaledaccuracy=scaledaccuracy/max(scaledaccuracy);
%%
textfontsz=20;
labfontsz=15;
figure
set(gcf,'Color','w')
marg=[.0,.055];
set(gcf,'PaperUnits','inches')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=w*ratio;

c=colormap;
spaced=linspace(0,1,size(c,1));
l=size(transformed,1);

xoffset=.25;
yoffset=.2;

numfun=6;

hold on
axis off
axis equal

for i=1:l
    p=transformed(i,:);
    mutnow=scaledmutvals(numfun,i);
    ind=sum(spaced<=mutnow);
    if ind==size(c,1)
        mutcol=c(end,:);
    else
        mutcol=c(ind,:)+(mutnow-spaced(ind))/(spaced(ind+1)-spaced(ind))*(c(ind+1,:)-c(ind,:));
    end
    fill([p(1),p(1)+cstep,p(1)+cstep,p(1),p(1)-cstep,p(1)-cstep],[p(2)-2/sqrt(3)*cstep,p(2)-1/sqrt(3)*cstep,p(2)+1/sqrt(3)*cstep,p(2)+2/sqrt(3)*cstep,p(2)+1/sqrt(3)*cstep,p(2)-1/sqrt(3)*cstep],mutcol,'EdgeColor','none');
end

text(0-xoffset,-yoffset,'Accuracy','Color','k','HorizontalAlignment','left','FontSize',textfontsz)
text(2+xoffset,-yoffset,'Dominance','Color','k','HorizontalAlignment','right','FontSize',textfontsz)
text(1,sqrt(3)+yoffset,'Time','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)
% text(0,sqrt(3),'MUTUAL','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)
% text(0,sqrt(3)-yoffset,'INFORMATION','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)

    
mybar=colorbar;
tickvec=0:.2:1;
ticklabvec=tickvec*(max(mutvals(numfun,:))-min(mutvals(numfun,:)))+min(mutvals(numfun,:));
ticklabvec=round(ticklabvec*10)/10;
set(mybar,'YTick',tickvec,'YTickLabel',ticklabvec,'FontSize',labfontsz)
yl=get(mybar,'YLabel');
v=get(yl,'Position');
set(yl,'String','Mutual information','FontSize',textfontsz,'Rotation',270,'Position',v+[3 0 0])

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w/10 -h/50 1.1*w 1.05*h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','mutinfo_heatmap','.pdf');
% print(filename,'-dpdf','-r300');

%%
textfontsz=20;
labfontsz=15;
figure
set(gcf,'Color','w')
marg=[.0,.055];
set(gcf,'PaperUnits','inches')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=w*ratio;

c=colormap;
spaced=linspace(0,1,size(c,1));
l=size(transformed,1);

xoffset=.25;
yoffset=.2;

numfun=6;

hold on
axis off
axis equal

for i=1:l
    p=transformed(i,:);
    skewnow=scaledskewvals(numfun,i);
    ind=sum(spaced<=skewnow);
    if ind==size(c,1)
        skewcol=c(end,:);
    else
        skewcol=c(ind,:)+(skewnow-spaced(ind))/(spaced(ind+1)-spaced(ind))*(c(ind+1,:)-c(ind,:));
    end
    fill([p(1),p(1)+cstep,p(1)+cstep,p(1),p(1)-cstep,p(1)-cstep],[p(2)-2/sqrt(3)*cstep,p(2)-1/sqrt(3)*cstep,p(2)+1/sqrt(3)*cstep,p(2)+2/sqrt(3)*cstep,p(2)+1/sqrt(3)*cstep,p(2)-1/sqrt(3)*cstep],skewcol,'EdgeColor','none');
end

text(0-xoffset,-yoffset,'Accuracy','Color','k','HorizontalAlignment','left','FontSize',textfontsz)
text(2+xoffset,-yoffset,'Dominance','Color','k','HorizontalAlignment','right','FontSize',textfontsz)
text(1,sqrt(3)+yoffset,'Time','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)
% text(0,sqrt(3),'MUTUAL','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)
% text(0,sqrt(3)-yoffset,'INFORMATION','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)

    
mybar=colorbar;
tickvec=0:.2:1;
ticklabvec=tickvec*(max(skewvals(numfun,:))-min(skewvals(numfun,:)))+min(skewvals(numfun,:));
ticklabvec=round(ticklabvec*10)/10;
set(mybar,'YTick',tickvec,'YTickLabel',ticklabvec,'FontSize',labfontsz)
yl=get(mybar,'YLabel');
v=get(yl,'Position');
set(yl,'String','Skewness','FontSize',textfontsz,'Rotation',270,'Position',v+[3 0 0])

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w/10 -h/50 1.1*w 1.05*h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','skewness_heatmap','.pdf');
% print(filename,'-dpdf','-r300');

%%
textfontsz=20;
labfontsz=15;
figure
set(gcf,'Color','w')
marg=[.0,.055];
set(gcf,'PaperUnits','inches')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=w*ratio;

c=colormap;
spaced=linspace(0,1,size(c,1));
l=size(transformed,1);

xoffset=.25;
yoffset=.2;

hold on
axis off
axis equal

for i=1:l
    p=transformed(i,:);
    accnow=scaledaccuracy(i);
    ind=sum(spaced<=accnow);
    if ind==size(c,1)
        acccol=c(end,:);
    else
        acccol=c(ind,:)+(accnow-spaced(ind))/(spaced(ind+1)-spaced(ind))*(c(ind+1,:)-c(ind,:));
    end
    fill([p(1),p(1)+cstep,p(1)+cstep,p(1),p(1)-cstep,p(1)-cstep],[p(2)-2/sqrt(3)*cstep,p(2)-1/sqrt(3)*cstep,p(2)+1/sqrt(3)*cstep,p(2)+2/sqrt(3)*cstep,p(2)+1/sqrt(3)*cstep,p(2)-1/sqrt(3)*cstep],acccol,'EdgeColor','none');
end

text(0-xoffset,-yoffset,'Accuracy','Color','k','HorizontalAlignment','left','FontSize',textfontsz)
text(2+xoffset,-yoffset,'Dominance','Color','k','HorizontalAlignment','right','FontSize',textfontsz)
text(1,sqrt(3)+yoffset,'Time','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)
% text(0,sqrt(3),'MUTUAL','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)
% text(0,sqrt(3)-yoffset,'INFORMATION','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)

    
mybar=colorbar;
tickvec=0:.2:1;
ticklabvec=tickvec*(max(accuracy)-min(accuracy))+min(accuracy);
ticklabvec=round(ticklabvec*10)/10;
set(mybar,'YTick',tickvec,'YTickLabel',ticklabvec,'FontSize',labfontsz)
yl=get(mybar,'YLabel');
v=get(yl,'Position');
set(yl,'String','Accuracy','FontSize',textfontsz,'Rotation',270,'Position',v+[3 0 0])

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w/10 -h/50 1.1*w 1.05*h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','accuracy_heatmap','.pdf');
% print(filename,'-dpdf','-r300');
