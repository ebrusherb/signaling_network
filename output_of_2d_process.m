green=importdata('green.txt');
blue=importdata('blue.txt');
red=importdata('red.txt');
black=importdata('black.txt');
%%
axislw=1.25;
lw=3;
labfontsz=15 ;
textfontsz=15;
alpha=.7;
ticklength=[.02,.0250];
markersz=5;
markercol='white';
%%
k=2;
i=3;

f2=interp1(reshape(diag(reshape(twomat(k,i,:,:,2),Nt,[])),1,[]),reshape(diag(reshape(twomat(k,i,:,:,1),Nt,[])),1,[]),'linear','pp');
f1=interp1(reshape(diag(reshape(onemat(k,i,:,:,2),Nt1,[])),1,[]),reshape(diag(reshape(onemat(k,i,:,:,1),Nt1,[])),1,[]),'linear','pp');
%%

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=w/2;
set(gcf,'Units','inches');
set(gcf,'Position',[v(1) v(2) w h]);
marg=[.12 .09];

textlabels={'a' 'b'};

p=zeros(1,Nd);

currblue=mycolors(Nd,'blue');
currred=mycolors(Nd,'red');
currblack=mycolors(Nd,'black');

for k=1:2
    subplot_tight(1,4,(1:2)+2*(k-1),marg)
    set(gca,'FontSize',labfontsz)
    hold on
    switch k
        case 2
            for i=3:2:Nd
            p(i)=plot(reshape(diag(reshape(twomat(k,i,:,:,2),Nt,[])),1,[]),reshape(diag(reshape(twomat(k,i,:,:,1),Nt,[])),1,[]),'Color',currblue(i,:),'LineWidth',lw);
%             p(i+Nd)=plot(reshape(diag(reshape(onemat(k,i,:,:,2),Nt1,[])),1,[]),reshape(diag(reshape(onemat(k,i,:,:,1),Nt1,[])),1,[]),'Color',currred(i,:),'LineWidth',lw);
            end
        otherwise 
            for i=2:2:Nd
            plot(reshape(diag(reshape(twomat(k,i,:,:,2),Nt,[])),1,[]),reshape(diag(reshape(twomat(k,i,:,:,1),Nt,[])),1,[]),'Color',currblue(i,:),'LineWidth',lw);
%             plot(reshape(diag(reshape(onemat(k,i,:,:,2),Nt1,[])),1,[]),reshape(diag(reshape(onemat(k,i,:,:,1),Nt1,[])),1,[]),'Color',currred(i,:),'LineWidth',lw);
            end
    end
    set(gca,'xlim',[0 5],'ylim',[.5 1]);
    ylabel('Probability of correct output','FontSize',textfontsz)
    xlabel('Expected time to output','FontSize',textfontsz)
    xlab=get(gca,'xlabel');
    v=get(xlab,'Position');
    set(xlab,'Position',[v(1) .46 v(3)]);
    text(-.4, 1.03,textlabels(k),'HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)
end

leg=legend([p([3:2:Nd])],['d=' num2str(domvals(3))],['d=' num2str(domvals(5))],['d=' num2str(domvals(7))],['d=' num2str(domvals(9))],['d=' num2str(domvals(Nd))]);
legend boxoff

% 
set(leg,'Position',[.8 .75 .1 .1])
%     
set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','2d_accuracy','.pdf');
print(filename,'-dpdf','-r300');

%%
i=1;
j=(Nd-1);

ratios=[.5 1 2];
Nr=length(ratios);

twoflexthresh=zeros(Nr,Nt,2);

for k=1:Nr
    r=ratios(k);
    twoflexthresh(k,:,1)=interp2(threshvals,threshvals,reshape(twomat(i,j,:,:,1),Nt,[]),r*threshvals,threshvals);
    twoflexthresh(k,:,2)=interp2(threshvals,threshvals,reshape(twomat(i,j,:,:,2),Nt,[]),r*threshvals,threshvals);
end


%%
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=w/2;
set(gcf,'Units','inches');
set(gcf,'Position',[v(1) v(2) w h]);

p=zeros(1,Nr);
currblue=mycolors(Nr+1,'blue');
currred=mycolors(Nr+1,'red');

% subplot_tight(1,2,1:2,[.1 .07])
set(gca,'FontSize',labfontsz)
hold on

for i=1:Nr
    p(i)=plot(reshape(twoflexthresh(i,:,2),1,[]),reshape(twoflexthresh(i,:,1),1,[]),'Color',currblue(Nr+2-i,:),'LineWidth',lw);
end

ylabel('Probability of correct output','FontSize',textfontsz)
xlabel('Expected time to output','FontSize',textfontsz)

legendlabels={};
for i=1:Nr
    legendlabels(i)=cellstr(['\rho=' num2str(ratios(i))]);
end

leg=legend(p,legendlabels);
legend boxoff

set(leg,'Position',[.85 .65 .1 .1],'FontSize',textfontsz)

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','2d_flexible_thresholds','.pdf');
print(filename,'-dpdf','-r300');