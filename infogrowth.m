N=20;
its=100;
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

indices=1:Nt;
leak=2;

% timevec=.1:.1:5;
% Ntime=length(timevec);
%%
mut_time=zeros(Nc,6,Ntime);

for i=1:Nc
    c1=X(i);
    c2=Y(i);
    c3=Z(i);
    tmax=max(max(timemats(:,:,i)));
    timevec=tmax/50:tmax/50:tmax;
    t=group_props_time(c1,c2,c3,its,N,'unif',Nd,Nt,twomat2,domvals,threshvals,leak,timevec);
    mut_time(i,:,:)=t;
end
%%
filename='/Users/eleanorbrush/Desktop/infogrowth_output.mat';
% save(filename,'mut_time');
%%
% diag=1:(N+1):N^2;
% offdiag=setdiff(1:(N^2),diag);
layerstoplot=[1 3 5];
Nlays=length(layerstoplot);
Ntime=100;
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
    t=group_props_time(c1,c2,c3,its,N,'unif',Nd,Nt,twomat2,domvals,threshvals,leak,timevec);
    mut_time(I,:,:)=t;
end
%% horizontal time plots look like crap
% figure;
% set(gcf,'Color','w')
% v=get(gcf,'Position');
% ratio=v(4)/v(3);
% w=6.83;
% h=.25*w;
% set(gcf,'Units','inches');
% set(gcf,'Position',[.5 1 w h]);
% % set(gca,'FontSize',labfontsz)
% marg=[.17 .09];
% xpos=.3;
% 
% layerstoplot=[1 3 5 7];
% for now=1:length(layerstoplot)
%     subplot_tight(1,4,now,marg)
%     j=layerstoplot(now);
%     k=layers{j}(1);
%     tmax=max(max(timemats(:,:,k)));
%     timevec=tmax/50:tmax/50:tmax;
% 
%     hold on
%     for i=1:4
%         I=funstoplot(i);
%         f=find(col(mut_time(k,I,:))>0,1,'first');
%         if f<=3
%             f=4;
%         end
%         plot(timevec(f-3:end),col(mut_time(k,I,f-3:end)),'-','Color',fourcols(colorperm(i),:),'LineWidth',1)
%     end
%     [~,w]=max(mut_time(k,funstoplot,end));
%     plot(timevec(f-3:end),col(mut_time(k,funstoplot(w),f-3:end)),'-','Color',fourcols(colorperm(w),:),'LineWidth',1)
%     v=get(gca,'Position');
%     set(gca,'Position',[(now-1)*.2 v(2)+.05 .2 v(4)]);
%     set(gca,'ylim',[0 1.25*max(mut_time(k,I,:))],'xlim',[timevec(f-3) 1.05*timevec(end)]);
%     set(gca,'FontSize',labfontsz,'FontName',fontname)
%     set(gca,'TickDir','out','ticklength',[.02 .025])
%     switch now
%     case 1
%     xlabel('Decision time','FontName',fontname,'FontSize',textfontsz);
%     ylabel('Skewness','FontName',fontname,'FontSize',textfontsz);
% 
%     case 4
%     [leg,hobj]=legend(funtypes);
%     legend('boxoff');
%     set(leg,'Position',[.75 .3 .3 .3])
%     textobj = findobj(hobj, 'type', 'line');
%     for i=[1 3 5 7]
%         set(textobj(i),'XData',[.22 .27])
%     end
%     end
% end
% 
% 
% set(gcf,'PaperSize',[w h]);
% set(gcf,'PaperPosition',[0 0 w h]);
% 
% filename=strcat('/Users/eleanorbrush/Desktop/','infogrowth','.pdf');
% print(filename,'-dpdf','-r300');

%% vertical time plots
figure;
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=1/3*6.83;
h=2/3*.75*6.83;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)
marg=[.07 .15];
xpos=.3;


layerstoplot=[1 3 7];
% linesty={'-.','--','-'};
linesty={'-','-','-'};
p=zeros(1,7);
for now=Nlays:-1:1
% for now=Nlays:Nlays 
    subplot_tight(3,1,now,marg)
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
    p(i+3)=plot((log(timevec)),[0; (col(mut_time(now,I,2:end)))],'-','Color',fourcolsnow(colorperm(i),:),'LineStyle',linesty{now},'LineWidth',lw);
    [~,w]=max(mut_time(now,funstoplot,end));
    plot((log(timevec)),(col(mut_time(now,funstoplot(w),:))),'-','Color',fourcolsnow(colorperm(w),:),'LineStyle',linesty{now},'LineWidth',lw);
    end
    set(gca,'ylim',[0 1.6])
    set(gca,'xlim',[log(sigfig(tmin)-3*step) log(sigfig(tmax)+step)])
    v=get(gca,'Position');
    set(gca,'Position',[.2 v(2:end)]);
    set(gca,'FontName',fontname,'FontSize',labfontsz)
%     timeticks=[.25 .5 1:.5:5];
    switch now
        case 3
        timeticks=.1:.01:.5;
        case 2
        timeticks=[.25 .5 1 1.5 2];
        case 1
        timeticks=[2:.5:4.5];
    end
    set(gca,'xtick',log(timeticks),'xticklabel',(timeticks))

end
xlabel('Time','FontName',fontname,'FontSize',textfontsz)
ylab=ylabel('Mutual information','FontName',fontname,'FontSize',textfontsz);
v=get(ylab,'Position');
set(ylab,'Position',[.55 v(2) v(3)]);
[leg,hobj]=legend(p(4:end),{'number of signals','number of signalers','entropy','eigenvector centrality'});
legend('boxoff')
set(leg,'Position',[.52 .73 .1 .1],'FontName',fontname)
textobj = findobj(hobj, 'type', 'line');
for i=[1 3 5 7]
    set(textobj(i),'XData',[.17 .27])
end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','infogrowth','.pdf');
print(filename,'-dpdf','-r300');
%% combined triangle and time plots
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
w=6.83;
h=.75*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)
marg=[.1 .063];
xpos=.3;

highlight='on';
tohighlight=[layers{1}(1) layers{3}(1) layers{5}(1)];
pointlabels={'b','c','d'};

subplot(3,4,[1:3 5:7 9:11])
% plot(1:10)
triimage_discrete2(scaledbestfun,bestfun,transformed,cornerlabels,funtypes2,fourcols,highlight,tohighlight,lw,pointlabels)
set(gca,'xlim',[0-2/11 2+2/11],'ylim',[0-1.7321/11 1.7321+1.7321/11])
v=get(gca,'Position');
ratio=v(4)/v(3);
scale=1;
set(gca,'Position',[.05 v(2) scale*v(3) scale*v(4)])

layerstoplot=[1 3 7];
% linesty={'-.','--','-'};
linesty={'-','-','-'};
p=zeros(1,7);
for now=1:Nlays
% for now=Nlays:Nlays

    subplot(3,3,3*(now-1)+3)
%     plot(1:10)
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
    p(i+3)=plot((log(timevec)),[0; (col(mut_time(now,I,2:end)))],'-','Color',fourcolsnow(colorperm(i),:),'LineStyle',linesty{now},'LineWidth',lw);
    end

    set(gca,'xlim',[log(sigfig(tmin)-3*step) log(sigfig(tmax)+step)])
%     v=get(gca,'Position');
%     set(gca,'Position',[.2 v(2:end)]);
    set(gca,'FontName',fontname,'FontSize',labfontsz)
%     timeticks=[.25 .5 1:.5:5];
    switch now
        case 3
        timeticks=.1:.01:.5;
        case 2
        timeticks=[.25 .5 1 1.5 2];
        case 1
        timeticks=[2:.5:4.5];
    end
    set(gca,'xtick',log(timeticks),'xticklabel',(timeticks))
    v=get(gca,'Position');
    set(gca,'Position',[.72 v(2:4)])
    ylim=get(gca,'Ylim');
    xlim=get(gca,'xlim');
    [x,y]=data2norm(xlim(1),ylim(2));
    annotation('textbox',[x-.1 y .05 .04],'String',['(' alphabet(now+1) ')'],'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
end
xlabel('Time','FontName',fontname,'FontSize',textfontsz)
ylab=ylabel('Mutual information','FontName',fontname,'FontSize',textfontsz);
v=get(ylab,'Position');
set(ylab,'Position',[log(.228) v(2)-.09 v(3)]);
[leg,hobj]=legend(p(4:end),{'number of signals','number of signalers','entropy','eigenvector centrality'});
legend('boxoff')
set(leg,'Position',[.05 .73 .1 .1],'FontName',fontname)
textobj = findobj(hobj, 'type', 'line');
for i=[1 3 5 7]
    set(textobj(i),'XData',[.17 .27])
end
annotation('textbox',[.1 .93 .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','info','.pdf');
print(filename,'-dpdf','-r300');