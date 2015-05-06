yellow=[1 1 0];
% red=[.9 0 0];
% green=[0 .7 0];
% blue=[0 0 .9];
% greens=importdata('green.txt');
% blues=importdata('blue.txt');
% reds=importdata('red.txt');
% blacks=importdata('black.txt');
mycolormap=cbrewer('seq','Reds',30);
alphabet=['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h'];
circlw=3;
fontname='Times';
textfontsz=12;
labfontsz=10;
circsz=13;
smallwidth=.05;
bigwidth=.3;
headlength=5;
headwidth=5;
textheight_up=.5;
textheight_down=.01;
topheight=.64;
rad=1;
circx=rad*cos(0:(2*pi/(N)):(N-1)/N*2*pi);
circy=rad*sin(0:(2*pi/(N)):(N-1)/N*2*pi);
pointrad=.11;
pointcircx=pointrad*cos(0:(2*pi/100):2*pi);
pointcircy=pointrad*sin(0:(2*pi/100):2*pi);


N=5;

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=w/2;
% set(gcf,'Units','inches');
% set(gcf,'Position',[v(1) v(2) w h]);
% set(gcf,'PaperSize',[w h]);
% set(gcf,'PaperPosition',[0 0 w h]);

subplot(2,4,1)
delx=.1; %#ok<*NASGU>
dely=.03;
hold on
for i=1:N
    colnow=mycolormap(5*(N+1-i),:);
    plot(1,N-i,'o','MarkerSize',circsz,'LineWidth',circlw,'Color',colnow,'MarkerFaceColor',colnow)
    text(1,N-i,['a_{',num2str(i),'}'],'FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','middle')
end
set(gca,'ylim',[0 N])

axis off
v=get(gca,'Position');
set(gca,'Position',[.05 topheight smallwidth v(4)])

v=get(gca,'Position');
annotation('textbox',[.0 textheight_up .2 .05],'string','1. Draw values.','EdgeColor','none','FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','bottom')

subplot(2,4,2)
hold on

for i=1:N
    colnow=mycolormap(5*(N+1-i),:);
    plot(circx(i),circy(i),'o','MarkerSize',circsz,'LineWidth',circlw,'Color',colnow,'MarkerFaceColor',colnow)
    text(circx(i),circy(i),['a_{',num2str(i),'}'],'FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','middle')
end

axis square
axis off
v=get(gca,'Position');
set(gca,'Position',[.05+smallwidth+.08 topheight bigwidth v(4)])

for i=1:N
    for j=(i+1):N
        normx=zeros(2,1);
        normy=zeros(2,1);
        [normx(1),normy(1)]=data2norm(circx(i),circy(i));
        [normx(2),normy(2)]=data2norm(circx(j),circy(j));
        distnow=norm([circx(i) circy(i)]-[circx(j) circy(j)]);
        scale=(distnow-pointrad*2)/distnow;
        normx=scalevec(normx,scale);
        normy=scalevec(normy,scale);
        annotation('doublearrow',normx,normy,'Head1Width',headwidth,'Head2Width',headwidth,'Head1Length',headlength,'Head2Length',headlength)
        if i==2 && j==3
            delx=.03;
            dely=.05;
            annotation('textbox',[mean(normx)-delx mean(normy)-dely .1 .1],'string','c_{23}','EdgeColor','none','FontName',fontname,'FontSize',textfontsz)
        end
    end
end

v=get(gca,'Position');
annotation('textbox',[v(1) textheight_up .25 .05],'string','2. Difficulty of pairwise decisions depends on difference in values.','EdgeColor','none','FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','bottom')

subplot(2,4,3)

delx=.1;
dely=.03;
hold on
for i=1:N
    colnow=mycolormap(5*(N+1-i),:);
    plot(1,N-i,'o','MarkerSize',circsz,'LineWidth',circlw,'Color',colnow,'MarkerFaceColor',colnow)
    text(1,N-i,['a_{',num2str(i),'}'],'FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','middle')
    text(2.5,N-i,['T_{',num2str(i),'}'],'FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','middle')
end
set(gca,'ylim',[0 N])

axis off
v=get(gca,'Position');
set(gca,'Position',[.55 topheight smallwidth v(4)])
v=get(gca,'Position');
annotation('textbox',[.5 textheight_up .2 .05],'string','3. Nash thresholds.','EdgeColor','none','FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','bottom')

subplot(2,4,4)

hold on

for i=1:N
    colnow=mycolormap(5*(N+1-i),:);
    plot(circx(i),circy(i),'o','MarkerSize',circsz,'LineWidth',circlw,'Color',colnow,'MarkerFaceColor',colnow)
    text(circx(i),circy(i),['a_{',num2str(i),'}'],'FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','middle')
end

axis square
axis off
v=get(gca,'Position');
set(gca,'Position',[.3+smallwidth*2+bigwidth topheight bigwidth v(4)])

for i=1:N
    for j=(i+1):N
        normx=zeros(2,1);
        normy=zeros(2,1);
        [normx(1),normy(1)]=data2norm(circx(i),circy(i));
        [normx(2),normy(2)]=data2norm(circx(j),circy(j));
        distnow=norm([circx(i) circy(i)]-[circx(j) circy(j)]);
        scale=(distnow-pointrad*2)/distnow;
        normx=scalevec(normx,scale);
        normy=scalevec(normy,scale);
        annotation('doublearrow',normx,normy,'Head1Width',headwidth,'Head2Width',headwidth,'Head1Length',headlength,'Head2Length',headlength)
        if i==2 && j==3
            delx=.07;
            dely=.04;
            annotation('textbox',[mean(normx)-delx mean(normy)-dely .1 .1],'string','p_{23}, DT_{23}','EdgeColor','none','FontName',fontname,'FontSize',textfontsz)
        end
    end
end

v=get(gca,'Position');

annotation('textbox',[v(1) textheight_up .3 .05],'string','4. Probability of edge forming each direction and expected time for edge to form.','EdgeColor','none','FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','bottom')

subplot(2,4,6)

hold on

for i=1:N
    colnow=mycolormap(5*(N+1-i),:);
    plot(circx(i),circy(i),'o','MarkerSize',circsz,'LineWidth',circlw,'Color',colnow,'MarkerFaceColor',colnow)
    text(circx(i),circy(i),['a_{',num2str(i),'}'],'FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','middle')
end

axis square
axis off

subplot(2,4,7)

delx=.1;
dely=.03;
hold on
for i=1:N
    colnow=mycolormap(5*(N+1-i),:);
    plot(1,N-i,'o','MarkerSize',circsz,'LineWidth',circlw,'Color',colnow,'MarkerFaceColor',colnow)
    text(1,N-i,['a_{',num2str(i),'}'],'FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','middle')
    text(2.5,N-i,['s_{',num2str(i),'}'],'FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','middle')
end
set(gca,'ylim',[0 N])

axis off
v=get(gca,'Position');
set(gca,'Position',[.62 v(2) smallwidth v(4)])
v=get(gca,'Position');
annotation('textbox',[.58 textheight_down .2 .05],'string','6. Consensus scores.','EdgeColor','none','FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','bottom')

subplot(2,4,6)
set(gca,'Position',[.05+smallwidth+.15 .07 bigwidth v(4)])
for i=1:N
    for j=(i+1):N
        normx=zeros(2,1);
        normy=zeros(2,1);
        [normx(1),normy(1)]=data2norm(circx(i),circy(i));
        [normx(2),normy(2)]=data2norm(circx(j),circy(j));
        distnow=norm([circx(i) circy(i)]-[circx(j) circy(j)]);
        scale=(distnow-pointrad*2)/distnow;
        normx=scalevec(normx,scale);
        normy=scalevec(normy,scale);
        annotation('arrow',normx,normy,'HeadWidth',headwidth,'HeadLength',headlength)
        if i==2 && j==3
            delx=.03;
            dely=.05;
            annotation('textbox',[mean(normx)-delx mean(normy)-dely .1 .1],'string','w_{23}','EdgeColor','none','FontName',fontname,'FontSize',textfontsz)
        end
    end
end

v=get(gca,'Position');
annotation('textbox',[v(1) textheight_down .3 .05],'string','5. Decision network forms.','EdgeColor','none','FontName',fontname,'FontSize',textfontsz,'HorizontalAlignment','center','VerticalAlignment','bottom')


% set(gcf,'PaperSize',[w h]);
% set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','diagram','.pdf');
print(filename,'-dpdf','-r300');
