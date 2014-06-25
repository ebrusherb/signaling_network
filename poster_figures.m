
%%
textfontsz=20;
labfontsz=15;
lw=5;
figure
set(gcf,'Color','w')

set(gcf,'PaperUnits','centimeters');
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=17.5;
h=.65*w;

xoffset=.01;

currblue=mycolors(Nd,'blue');
currred=mycolors(Nd,'red');
currblack=mycolors(Nd,'black');


set(gca,'FontSize',labfontsz)
hold on
k=1;
toplot=3:2:Nd;
p=zeros(1,length(toplot));
leglabels=cell(1,length(toplot));

for I=1:length(toplot)
    i=toplot(I);
    p(I)=plot(reshape(diag(reshape(twomat(k,i,:,:,2),Nt,[])),1,[]),reshape(diag(reshape(twomat(k,i,:,:,1),Nt,[])),1,[]),'Color',currblue(i,:),'LineWidth',lw);
    leglabels{I}=['c=' num2str(domvals(i))];
end

set(gca,'ylim',[.5 1]);
yl=ylabel('Accuracy','FontSize',textfontsz);
xlabel('Time to decision','FontSize',textfontsz)
v=get(yl,'Position');
set(yl,'Position',v-[xoffset 0 0]);

leg=legend(p(end:-1:1),leglabels{end:-1:1});
legend boxoff
set(leg,'Position',[.85 .75 .1 .1])
    
set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w/50 0 .9*w 1.035*h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','2d_accuracy_poster','.pdf');
% print(filename,'-dpdf','-r300');
%%
textfontsz=20;
labfontsz=15;
figure
set(gcf,'Color','w')
marg=[.0,.055];
set(gcf,'PaperUnits','centimeters')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=17.5;
h=w*ratio;

c=colormap(jet);
spaced=linspace(0,1,size(c,1));
l=size(transformed,1);

xoffset=.25;
yoffset=.2;

hold on
axis off
axis equal

for i=1:l    
    p=transformed(i,:);
    skewnow=scaledskewvals(i);
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
% text(0,sqrt(3),'SKEWNESS','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)
    
mybar=colorbar;
tickvec=0:.2:1;
ticklabvec=tickvec*(max(skewvals)-min(skewvals))+min(skewvals);
ticklabvec=round(ticklabvec*10)/10;
ticklabvec(1)=0;
set(mybar,'YTick',tickvec,'YTickLabel',ticklabvec,'FontSize',labfontsz)
yl=get(mybar,'YLabel');
v=get(yl,'Position');
set(yl,'String','Skewness','FontSize',textfontsz,'Rotation',270,'Position',v+[3 0 0])

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w/10 -h/50 1.1*w 1.05*h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','skewness_heatmap_poster','.pdf');
% print(filename,'-dpdf','-r300');

%%
textfontsz=20;
labfontsz=15;
figure
set(gcf,'Color','w')
marg=[.0,.055];
set(gcf,'PaperUnits','centimeters')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=17.5;
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
    mutnow=scaledmutvals(i);
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
ticklabvec=tickvec*(max(mutvals)-min(mutvals))+min(mutvals);
ticklabvec=round(ticklabvec*10)/10;
set(mybar,'YTick',tickvec,'YTickLabel',ticklabvec,'FontSize',labfontsz)
yl=get(mybar,'YLabel');
v=get(yl,'Position');
set(yl,'String','Mutual information','FontSize',textfontsz,'Rotation',270,'Position',v+[3 0 0])

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w/10 -h/50 1.1*w 1.05*h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','mutinfo_heatmap_poster','.pdf');
% print(filename,'-dpdf','-r300');

%%
%TWO ROWS
textfontsz=15;
labfontsz=10;
lw=2;
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=11;
h=w/4;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
marg=[.12 .063];
xoffsetperc=.08;
yoffsetperc=.14;

l=1;

ivals=[1 2 11];

for I=1:length(ivals)
    i=ivals(I);
    v=zeros(2,Nd);
    for j=2:Nd
        N=nasheq_ind2{i,l,j};
        v(1,j)=threshvals(N(1,:));
        v(2,j)=threshvals(N(2,:));
    end
    subplot_tight(2,length(ivals),I,marg)
    hold on
    plot(domvals(2:Nd),v(1,2:Nd),'Color','blue','LineWidth',lw)
    plot(domvals(2:Nd),v(2,2:Nd),'--r','LineWidth',lw)
    set(gca,'YLim',[0 max(threshvals)+threshvals(1)]);
    set(gca,'XLim',[min(domvals) max(domvals)]);
    set(gca,'FontSize',labfontsz)
    switch I
        case 1
            yl=ylabel('Threshold','FontSize',textfontsz);
            v=get(yl,'Position');
            Xvec=get(gca,'xlim');
            xdiff=diff(Xvec);
            xpos=Xvec(1)-xoffsetperc*xdiff;
            set(yl,'Position',[xpos v(2) v(3)])
            
            leg=legend('Stronger','Weaker');
            legend boxoff
            set(leg,'Position',[.1 .63 .1 .1])
            legendshrink(.5,'best',leg)
    end
    box on
end

for I=1:length(ivals)
   i=ivals(I);
    vprob=zeros(1,Nd);
    vtime=zeros(1,Nd);
    for j=2:Nd
        N=nasheq_ind2{i,l,j};
        vprob(j)=twomat(l,j,N(1,:),N(2,:),1);
        vtime(j)=twomat(l,j,N(1,:),N(2,:),2);
    end
    subplot_tight(2,length(ivals),I+length(ivals),marg)
    plot(domvals(2:Nd),vtime(1,2:Nd),'k','LineWidth',lw)
%     set(gca,'ylim',[0 4])
    set(gca,'XLim',[min(domvals) max(domvals)]);
    set(gca,'FontSize',labfontsz)
    yvec=get(gca,'ylim');
    ydiff=diff(yvec);
    xl=xlabel('c','FontSize',textfontsz);
    v=get(xl,'Position');
    ypos=yvec(1)-yoffsetperc*ydiff;
    set(xl,'Position',[v(1) ypos v(2)])
    switch I
        case 1
            yl=ylabel('Time','FontSize',textfontsz);
            v=get(yl,'Position');
            Xvec=get(gca,'xlim');
            xdiff=diff(Xvec);
            xpos=Xvec(1)-xoffsetperc*xdiff;
            set(yl,'Position',[xpos v(2) v(3)])
    end

end

annotation('arrow',[.2 .8],[.915 .915])
annotation('textbox',[.25 0.95 .5 .072],'string','Increasing cost of waiting','FontSize',textfontsz,'EdgeColor','none','HorizontalAlignment','center')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','nasheq_thresholds_poster','.pdf');
print(filename,'-dpdf','-r300');

%%
%ONE ROW
textfontsz=15;
labfontsz=10;
lw=2;
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=11;
h=w/6;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
marg=[.2 .063];
xoffsetperc=.08;
yoffsetperc=.12;

l=1;

ivals=[1 2 11];

for I=1:length(ivals)
    i=ivals(I);
    v=zeros(2,Nd);
    for j=2:Nd
        N=nasheq_ind2{i,l,j};
        v(1,j)=threshvals(N(1,:));
        v(2,j)=threshvals(N(2,:));
    end
    subplot_tight(1,length(ivals),I,marg)
    hold on
    plot(domvals(2:Nd),v(1,2:Nd),'Color','blue','LineWidth',lw)
    plot(domvals(2:Nd),v(2,2:Nd),'--r','LineWidth',lw)
    set(gca,'YLim',[0 max(threshvals)+threshvals(1)]);
    set(gca,'XLim',[min(domvals) max(domvals)]);
    set(gca,'FontSize',labfontsz)
    yvec=get(gca,'ylim');
    ydiff=diff(yvec);
    xl=xlabel('c','FontSize',textfontsz);
    v=get(xl,'Position');
    ypos=yvec(1)-yoffsetperc*ydiff;
    set(xl,'Position',[v(1) ypos v(2)])
    switch I
        case 1
            yl=ylabel('Threshold','FontSize',textfontsz);
            v=get(yl,'Position');
            Xvec=get(gca,'xlim');
            xdiff=diff(Xvec);
            xpos=Xvec(1)-xoffsetperc*xdiff;
            set(yl,'Position',[xpos v(2) v(3)])
            
            leg=legend('Stronger','Weaker');
            legend boxoff
            set(leg,'Position',[.1 .4 .1 .1])
            legendshrink(.5,'best',leg)
    end
    box on
end

annotation('arrow',[.2 .8],[.87 .87])
annotation('textbox',[.25 0.95 .5 .072],'string','Increasing waiting cost','FontSize',textfontsz,'EdgeColor','none','HorizontalAlignment','center')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','nasheq_thresholds_poster','.pdf');
print(filename,'-dpdf','-r300');
%%
textfontsz=15;
labfontsz=10;
lw=2;
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=11;
h=w/4;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
marg=[.12 .063];
xoffsetperc=.08;
yoffsetperc=.17;

l=1;

nrows=2;

ivals=[1 2 4];

N=size(timemats(:,:,1),1);

c=colormap;
spaced=linspace(0,1,size(c,1));

for I=1:length(ivals)
    i=ivals(I);
    v=avgthresh(:,i);
    subplot_tight(nrows,length(ivals),I,marg)
    hold on
    plot(v,'Color','k','LineWidth',lw)
    set(gca,'XLim',[1 N],'YLim',[0 max(threshvals)+threshvals(1)]);
    set(gca,'FontSize',labfontsz)
    set(gca,'YLim',[0 max(threshvals)+threshvals(1)]);
    set(gca,'FontSize',labfontsz)
    switch I
        case 1
            yl=ylabel('Threshold','FontSize',textfontsz);
            v=get(yl,'Position');
            Xvec=get(gca,'xlim');
            xdiff=diff(Xvec);
            xpos=Xvec(1)-xoffsetperc*xdiff;
            set(yl,'Position',[xpos v(2) v(3)])
    end
    box on
end

clims=[0 ceil(max(max(max(timemats))))];

for I=1:length(ivals)
    i=ivals(I);
    subplot_tight(nrows,length(ivals),[I+length(ivals)] ,marg)
%     imagesc(timemats(:,:,i))
    set(gca,'xlim',[.5 20.5],'ylim',[.5 20.5])
    hold on
    
    scaledtime=timemats(:,:,i);
    scaledtime(~scaledtime)=nan;
    scaledtime=timemats(:,:,i)-min(min(scaledtime));
    scaledtime=scaledtime/max(max(scaledtime));
    for k=1:N
        for l=[1:(k-1) (k+1):N]
            timenow=scaledtime(k,l);
            ind=sum(spaced<=timenow);
            if ind==size(c,1)
                timecol=c(end,:);
            else
                timecol=c(ind,:)+(timenow-spaced(ind))/(spaced(ind+1)-spaced(ind))*(c(ind+1,:)-c(ind,:));
            end
            fill([k-.5 k+.5 k+.5 k-.5],[l-.5 l-.5 l+.5 l+.5],timecol,'EdgeColor','none')
        end
    end
    set(gca,'FontSize',labfontsz)
    yvec=get(gca,'ylim');
            ydiff=diff(yvec);
            xl=xlabel('True order','FontSize',textfontsz);
            v=get(xl,'Position');
            ypos=yvec(1)-yoffsetperc*ydiff;
            set(xl,'Position',[v(1) ypos v(2)])
    switch I
        case 1
            yl=ylabel('Time','FontSize',textfontsz);
            v=get(yl,'Position');
            Xvec=get(gca,'xlim');
            xdiff=diff(Xvec);
            xpos=Xvec(1)-xoffsetperc*xdiff;
            set(yl,'Position',[xpos v(2) v(3)])
        case 2
%         case 3
%             mybar=colorbar;
%             set(mybar,'Position',[.95 .0976 .0229 .3523])
    end
    
    box on
    
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end

annotation('arrow',[.2 .8],[.915 .915])
annotation('textbox',[.25 0.95 .5 .072],'string','Increasing waiting cost','FontSize',textfontsz,'EdgeColor','none','HorizontalAlignment','center')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','groupeq_thresholds_poster','.pdf');
print(filename,'-dpdf','-r300');

