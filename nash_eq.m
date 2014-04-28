
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

c1=0;
c3=1;

c2vals=0:.05:1;
Nc2=length(c2vals);

nasheq_ind2=cell(Nc2,Nl,Nd);
giveup2=zeros(Nc2,Nl,Nd);

objective_best2=cell(Nc2,Nl,Nd);

for i=1:Nc2
    c2=c2vals(i);
    perf2=zeros(2,Nl,Nd,Nt,Nt); %individual, leak, dominance, thresholds

    perf2(1,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2))+c3*(1-twomat(:,2:Nd,:,:,1));
    perf2(2,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2))+c3*(twomat(:,2:Nd,:,:,1));

    perf2(1,:,1,:,:)=c2*(twomat(:,1,:,:,2))+c3*(1-twomat(:,1,:,:,1));
    perf2(2,:,1,:,:)=c2*(twomat(:,1,:,:,2))+c3*(twomat(:,1,:,:,1));
    
    objective_perf2=zeros(Nl,Nd,Nt,Nt);
    
    objective_perf2(:,2:Nd,:,:)=max([c1 c3])*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2));
    
    objective_perf2(:,1,:,:)=c2*(twomat(:,1,:,:,2));
    
%     for j=1:Nl
    for j=1:2
        for k=1:Nd
            T1_best=zeros(1,Nt); 
            T2_best=zeros(1,Nt);
            for u=1:Nt
            [~,m]=min(perf2(1,j,k,:,u));
            T1_best(u)=m;
            [~,m]=min(perf2(2,j,k,u,:));
            T2_best(u)=m;
            end
            pp1=interp1(indices,reshape(T1_best,1,[]),'linear','pp');
            pp2=interp1(indices,reshape(T2_best,1,[]),'linear','pp');
            v1=ppval(pp2,indices);
            v2=ppval(pp1,v1);
            T1=indices(abs((indices)-v2)<0.01);
            T2=ppval(pp2,T1);
            nasheq_ind2{i,j,k}=[T1;T2];
            
            if size(nasheq_ind2{i,j,k},2)==1
                giveup2(i,j,k)=(sum(nasheq_ind2{i,j,k}==[Nt;1])==2);
            end
            
            [minval1 where1]=min(reshape(objective_perf2(j,k,:,:),Nt,[]));
            [minval2 where2]=min(minval1);
            objective_best2{i,j,k}=[where1(where2); where2];
        end
    end
end


%%
whentogiveup=zeros(Nc2,Nl);

for i=1:Nc2
    for j=1:2
        whentogiveup(i,j)=domvals(sum(giveup2(i,j,:)==0));
    end
end

%%
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
xpos=.4;

l=1;

ivals=[1 2 5 11];

for I=1:length(ivals)
    i=ivals(I);
    v=zeros(2,Nd);
    for j=2:Nd
        N=nasheq_ind2{i,l,j};
        v(1,j)=threshvals(N(1,:));
        v(2,j)=threshvals(N(2,:));
    end
    subplot_tight(3,length(ivals),I,marg)
    hold on
    plot(domvals(2:Nd),v(1,2:Nd),'Color','blue','LineWidth',lw)
    plot(domvals(2:Nd),v(2,2:Nd),'--r','LineWidth',lw)
    set(gca,'YLim',[0 max(threshvals)+threshvals(1)]);
    set(gca,'FontSize',labfontsz)
    switch I
        case 1
            ylabel('Threshold','FontSize',textfontsz)
            v=get(get(gca,'ylabel'),'Position');
            set(get(gca,'ylabel'),'Position',[xpos v(2) v(3)])
            leg=legend('Stronger','Weaker');
            legend boxoff
            set(leg,'Position',[.1 .72 .1 .1])
            legendshrink(.5,'best',leg)
%         case 4
%             xlabel('Probability stronger animal wins','FontSize',textfontsz)
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
    subplot_tight(3,length(ivals),I+length(ivals),marg)    
    plot(domvals(2:Nd),vprob(1,2:Nd),'k','LineWidth',lw)
    set(gca,'ylim',[.5 1])
    set(gca,'FontSize',labfontsz)
    
    switch I
        case 1
            ylabel('Accuracy','FontSize',textfontsz)
            v=get(get(gca,'ylabel'),'Position');
            set(get(gca,'ylabel'),'Position',[xpos v(2) v(3)])
%         case 4
%             xlabel('Probability stronger animal wins','FontSize',textfontsz)
    end
   
    subplot_tight(3,length(ivals),I+length(ivals)*2,marg)
    plot(domvals(2:Nd),vtime(1,2:Nd),'k','LineWidth',lw)
%     set(gca,'ylim',[0 4])
    set(gca,'FontSize',labfontsz)
    switch I
        case 1
            ylabel('Time','FontSize',textfontsz)
            v=get(get(gca,'ylabel'),'Position');
            set(get(gca,'ylabel'),'Position',[xpos v(2) v(3)])
%         case 3
%             xlabel('Probability stronger animal wins','FontSize',textfontsz)
    end

end

annotation('textbox',[.25 0 .5 .072],'string','Probability stronger animal wins','FontSize',textfontsz,'EdgeColor','none','HorizontalAlignment','center')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','nasheq_thresholds','.pdf');
print(filename,'-dpdf','-r300');
