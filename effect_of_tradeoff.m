function effect_of_tradeoff(c1,c2,c3)

Nc=evalin('base','Nc');
Nl=evalin('base','Nl');
Nd=evalin('base','Nd');
Nt=evalin('base','Nt');
twomat=evalin('base','twomat');
threshvals=evalin('base','threshvals');
indices=1:Nt;
domvals=evalin('base','domvals');
blue=evalin('base','blue');
red=evalin('base','red');
markersz=evalin('base','markersz');
lw=evalin('base','lw');

perf=zeros(Nc,2,Nl,Nd,Nt,Nt); %criteria, individual, leak, dominance, thresholds

perf(1,1,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2));
perf(1,2,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2));
perf(2,1,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2))+c3*(1-twomat(:,2:Nd,:,:,1));
perf(2,2,:,2:Nd,:,:)=c1*(1-twomat(:,2:Nd,:,:,1))+c2*(twomat(:,2:Nd,:,:,2))+c3*(twomat(:,2:Nd,:,:,1));

perf(1,1,:,1,:,:)=c2*(twomat(:,1,:,:,2));
perf(1,2,:,1,:,:)=c2*(twomat(:,1,:,:,2));
perf(2,1,:,1,:,:)=c2*(twomat(:,1,:,:,2))+c3*(1-twomat(:,1,:,:,1));
perf(2,2,:,1,:,:)=c2*(twomat(:,1,:,:,2))+c3*(twomat(:,1,:,:,1));

nasheq_thresh=cell(Nc,Nl,Nd);
nasheq_ind=cell(Nc,Nl,Nd);

T1_best=zeros(Nc,Nl,Nd,Nt); %criteria, leak, dominance, threshold
T2_best=zeros(Nc,Nl,Nd,Nt);

for i=1:Nc
    for j=1:Nl
        for k=1:Nd
            for u=1:Nt
            [~,m]=min(perf(i,1,j,k,:,u));
            T1_best(i,j,k,u)=m;
            [~,m]=min(perf(i,2,j,k,u,:));
            T2_best(i,j,k,u)=m;
            end
            pp1=interp1(threshvals,reshape(threshvals(T1_best(i,j,k,:)),1,[]),'linear','pp');
            pp2=interp1(threshvals,reshape(threshvals(T2_best(i,j,k,:)),1,[]),'linear','pp');
            v1=ppval(pp2,threshvals);
            v2=ppval(pp1,v1);
            T1=threshvals(abs(threshvals-v2)<0.01);
            T2=ppval(pp2,T1);
            nasheq_thresh{i,j,k}=[T1;T2];
            pp1=interp1(indices,reshape(T1_best(i,j,k,:),1,[]),'linear','pp');
            pp2=interp1(indices,reshape(T2_best(i,j,k,:),1,[]),'linear','pp');
            v1=ppval(pp2,indices);
            v2=ppval(pp1,v1);
            T1=indices(abs((indices)-v2)<0.01);
            T2=ppval(pp2,T1);
            nasheq_ind{i,j,k}=[T1;T2];
        end
    end
end

figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=5*2;
h=w*2;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)

p=zeros(1,Nl*2);

for i=1:Nl
    for j=1:Nd
        subplot(2,2,i)
        hold on
        L=size(nasheq_thresh{2,i,j},2);
        d=domvals(j);
        num=plot(d*ones(1,L),nasheq_thresh{2,i,j}(1,:),'LineStyle','o','Color','blue','LineWidth',lw,'MarkerSize',markersz+5);
        p(2*i)=num(1);
        num=plot(d*ones(1,L),nasheq_thresh{2,i,j}(2,:),'LineStyle','o','Color','red','LineWidth',lw,'MarkerSize',markersz);
        p(2*i-1)=num(1);
    end
    set(gca,'YLim',[0 max(threshvals)]);
end

for i=1:Nl
    subplot(2,2,i+2)
    d=domvals(1);
    N=nasheq_ind{2,i,1};
    [hax hline1 hline2]=plotyy(d*ones(1,L),twomat(i,1,N(1,1),N(2,1),1),d*ones(1,L),twomat(i,1,N(1,1),N(2,1),2));
    set(hline1,'LineStyle','o','Color','red','LineWidth',lw,'MarkerSize',markersz);
    set(hline2,'LineStyle','o','Color','blue','LineWidth',lw,'MarkerSize',markersz);

    axes(hax(1)); hold on
    for j=1:Nd
        N=nasheq_ind{2,i,j};
        L=size(N,2);
        d=domvals(j);
        for k=1:L
            plot(d*ones(1,L),twomat(i,j,N(1,k),N(2,k),1),'LineStyle','o','Color','red','MarkerSize',markersz,'LineWidth',lw)
        end
    end

    axes(hax(2)); hold on
    for j=2:Nd
        N=nasheq_ind{2,i,j};
        L=size(N,2);
        d=domvals(j);
        for k=1:L
            plot(d*ones(1,L),twomat(i,j,N(1,k),N(2,k),2),'LineStyle','o','Color','blue','MarkerSize',markersz,'LineWidth',lw)
        end
    end
    
    set(hax(1),'xlim',[.5 1],'ylim',[0 1],'ytick',0:.1:1)
    set(hax(2),'xlim',[.5 1],'ylim',[0 10],'ytick',0:1:10)
    set(hax,{'ycolor'},{'r';'b'})
end
end