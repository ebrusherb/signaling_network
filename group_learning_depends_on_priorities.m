its=10;
N=20;

c2vals=[0 .05 .3 .5];
c2vals=[c2vals 1];
Nc2=length(c2vals);

c1=0;
c3=1;

avgthresh=zeros(N,Nc2);
indmutinfo=zeros(N,Nc2);
groupmutinfo=zeros(1,Nc2);
groupskew=zeros(1,Nc2);
threshtohist=zeros(N,Nc2);
probmats=zeros(N,N,Nc2);
timemats=zeros(N,N,Nc2);


% toreturn={avgthresh,threshtohist,probmat,timemat,indmutinfos,groupmutinfo,meanskewness};

for i=1:(Nc2-1)
    c2=c2vals(i);
    c3=1-c2;
    t=group_props(c1,c2,c3,its,N,'unif');
    avgthresh(:,i)=t{1};
    threshtohist(:,i)=t{2};
    probmats(:,:,i)=t{3};
    timemats(:,:,i)=t{4};
    indmutinfo(:,i)=t{5};
    groupmutinfo(i)=t{6};
    groupskew(i)=t{7};
end

i=Nc2;
c2=c2vals(i);
t=group_props(c1,c2,0,10,N,'unif');
avgthresh(:,i)=t{1};
threshtohist(:,i)=t{2};
probmats(:,:,i)=t{3};
timemats(:,:,i)=t{4};
indmutinfo(:,i)=t{5};
groupmutinfo(i)=t{6};
groupskew(i)=t{7};
%%
plot(indmutinfo(:,1))
hold on
plot(indmutinfo(:,2))
plot(indmutinfo(:,3))
plot(indmutinfo(:,4))
%%
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
nrows=3;
h=nrows/4*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);

marg=[.1 .065];
xpos=-10;



ivals=[1 2 3 4];

for I=1:length(ivals)
    i=ivals(I);
    v=avgthresh(:,i);
    subplot_tight(nrows,length(ivals),I,marg)
    hold on
    plot(v,'Color','k','LineWidth',lw)
    set(gca,'YLim',[0 max(threshvals)+threshvals(1)]);
    set(gca,'FontSize',labfontsz)
    titlelab=['w2=',num2str(c2vals(i))];
    title(titlelab,'FontSize',textfontsz);
    switch I
        case 1
            ylabel('Threshold','FontSize',textfontsz)
            v=get(get(gca,'ylabel'),'Position');
            set(get(gca,'ylabel'),'Position',[xpos v(2) v(3)])
%         case 2
%             xlabel('True order of fighting ability','FontSize',textfontsz)
    end
    box on
end

for I=1:length(ivals)
    i=ivals(I);
    v=threshtohist(:,i);
    subplot_tight(nrows,length(ivals),I+length(ivals),marg)
    hold on
%     hist(v,0:.25:2)
%     set(gca,'XLim',[0 max(threshvals)+threshvals(1)],'YLim',[0 N]);
%     set(gca,'FontSize',labfontsz)
    plot(indmutinfo(:,i),'k','LineWidth',lw)
    set(gca,'YLim',[0 max(max(indmutinfo))+.05])
    set(gca,'FontSIze',labfontsz)
    switch I
            case 1
                ylabel('Mutual Info','FontSize',textfontsz)
                v=get(get(gca,'ylabel'),'Position');
                set(get(gca,'ylabel'),'Position',[xpos v(2) v(3)])
%             case 2
%                 xlabel('True order of fighting ability','FontSize',textfontsz)
        end
    box on
end

% for I=1:length(ivals)
%     i=ivals(I);
%     v=avgthresh(:,i);
%     subplot_tight(nrows,length(ivals),I+2*length(ivals),marg)
%     imagesc(probmats(:,:,i),[0 1])
%     set(gca,'FontSize',labfontsz)
%     switch I
%         case 1
% %             ylabel('Order of fighting ability','FontSize',textfontsz)
% %             v=get(get(gca,'ylabel'),'Position');
% %             set(get(gca,'ylabel'),'Position',[xpos+.5 v(2) v(3)])
% %         case 3
% %             xlabel('True order of fighting ability','FontSize',textfontsz)
% %         case 3
% %             mybar=colorbar;
% %             set(mybar,'Position',[.95 .0976 .0229 .3523])
%     end
%     box on
% end

clims=[0 ceil(max(max(max(timemats))))];

for I=1:length(ivals)
    i=ivals(I);
    v=avgthresh(:,i);
    subplot_tight(nrows,length(ivals),I+2*length(ivals),marg)
    imagesc(timemats(:,:,i))
    set(gca,'FontSize',labfontsz)
    switch I
        case 1
            ylabel('Time','FontSize',textfontsz)
            v=get(get(gca,'ylabel'),'Position');
            set(get(gca,'ylabel'),'Position',[xpos+.5 v(2) v(3)])
%         case 3
%             xlabel('True order of fighting ability','FontSize',textfontsz)
%         case 3
%             mybar=colorbar;
%             set(mybar,'Position',[.95 .0976 .0229 .3523])
    end
    box on
end

annotation('textbox',[.25 0 .5 .072],'string','True order of fighting ability','FontSize',textfontsz,'EdgeColor','none','HorizontalAlignment','center')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','groupeq_thresholds','.pdf');
% print(filename,'-dpdf','-r300');


%%
textfontsz=20;
labfontsz=15;
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
nrows=1;
h=nrows/2.5*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);

marg=[.2 .1];
xpos=-10;

ivals=[1 2 4];

for I=1:length(ivals)
    i=ivals(I);
    v=avgthresh(:,i);
    subplot_tight(nrows,length(ivals),I,marg)
    hold on
    plot(v,'Color','k','LineWidth',lw)
    set(gca,'YLim',[0 max(threshvals)+threshvals(1)]);
    set(gca,'FontSize',labfontsz)
    switch I
        case 1
            ylabel('Threshold','FontSize',textfontsz)
            v=get(get(gca,'ylabel'),'Position');
            set(get(gca,'ylabel'),'Position',[xpos v(2) v(3)])
%         case 2
%             xlabel('True order of fighting ability','FontSize',textfontsz)
    end
    box on
end

annotation('textbox',[.15 0.06 .7 .072],'string','True order of fighting ability','FontSize',textfontsz,'EdgeColor','none','HorizontalAlignment','center')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','groupeq_thresholds','.pdf');
print(filename,'-dpdf','-r300');

%%

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=w/2;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
set(gca,'FontSize',labfontsz)
marg=[.15 .1];

subplot_tight(1,2,1,marg)
plot(c2vals(1:(end-1)),groupskew(1:(end-1)),'Color','k','LineWidth',lw)
% xlabel('Cost of fighting','FontSize',textfontsz)
ylabel('Expected skewness','FontSize',textfontsz)

subplot_tight(1,2,2,marg)
plot(c2vals(1:(end-1)),groupmutinfo(1:(end-1)),'Color','k','LineWidth',lw)
% xlabel('Cost of fighting','FontSize',textfontsz)
ylabel('Mutual information','FontSize',textfontsz)

annotation('textbox',[.25 0 .5 .072],'string','Cost of fighting','FontSize',textfontsz,'EdgeColor','none','HorizontalAlignment','center')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','group_learning','.pdf');
print(filename,'-dpdf','-r300');
