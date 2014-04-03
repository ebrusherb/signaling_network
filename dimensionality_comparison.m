deltax1=0.01;
deltax2=0.025 ;
b=1;

lvals=0:.5:2;
Nl=length(lvals);

domvals=.5:.1:1;
Nd=length(domvals);

threshvals1=.1:.1:8;
Nt1=length(threshvals1);

threshvals2=.5:.5:2;
Nt2=length(threshvals2);

%%

oneoutput=zeros(Nl,Nd,Nt1,2);
twooutput=zeros(Nl,Nd,Nt2,2);

for k=1:Nl
    l=lvals(k);
    for i=1:Nd
        d=domvals(i);
        for j=1:Nt1
            t=threshvals1(j);
            [x,y,~]=solve_pde_1d(-t,d,t,deltax1,l,b);
            oneoutput(k,i,j,:)=[x y];
        end
        for j=1:Nt2
            t=threshvals2(j);
            [x,y,~]=solve_pde_2d(-t,-t,d,5,5,deltax2,deltax2,l,b);
            twooutput(k,i,j,:)=[x y];
        end
    end
end
%%
filename=strcat('two_v_one_output','.mat');
save(filename,'twooutput','oneoutput','domvals','lvals','threshvals1','threshvals2');
%%
axislw=1.25;
lw=3;
labfontsz=15 ;
textfontsz=15;
alpha=.7;
ticklength=[.02,.0250];
markersz=5;
markercol='white';
rho=4.2/5.6;
w=8.5;
h=w*rho;

for k=Nl:-1:1
    figure
    set(gcf,'Color','w')
    v=get(gcf,'Position');
    ratio=v(4)/v(3);
    w=6.83*2;
    h=w*ratio;
    hold on
    for i=2:Nd
    plot(reshape(twooutput(k,i,:,1),1,[]),reshape(twooutput(k,i,:,2),1,[]),'blue','LineWidth',lw)
    plot(reshape(oneoutput(k,i,1:40,1),1,[]),reshape(oneoutput(k,i,1:40,2),1,[]),'red','LineWidth',lw)
    end
    xlabel('Probability of correct output','FontSize',textfontsz)
    ylabel('Expected time to output','FontSize',textfontsz)
%     leg=legend('Two dimensions','One dimension');
%     set(leg,'Position',[.6 .8 .1 .05])
end

%%
axislw=1.25;
lw=3;
labfontsz=15 ;
textfontsz=15;
alpha=.7;
ticklength=[.02,.0250];
markersz=5;
markercol='white';
rho=4.2/5.6;
w=8.5;
h=w*rho;

figure
    set(gcf,'Color','w')
    v=get(gcf,'Position');
    ratio=v(4)/v(3);
    w=6.83*2;
    h=w*ratio;
    

for j=1:2
    subplot(1,2,j)
    k=j;
    for i=2:Nd
    plot(reshape(twooutput(k,i,:,1),1,[]),reshape(twooutput(k,i,:,2),1,[]),'blue','LineWidth',lw)
    hold on
    plot(reshape(oneoutput(k,i,1:40,1),1,[]),reshape(oneoutput(k,i,1:40,2),1,[]),'red','LineWidth',lw)
    end
    xlabel('Probability of correct output','FontSize',textfontsz)
    ylabel('Expected time to output','FontSize',textfontsz)

end
%%
axislw=1.25;
lw=3;
labfontsz=15 ;
textfontsz=15;
alpha=.7;
ticklength=[.02,.0250];
markersz=5;
markercol='white';
rho=4.2/5.6;
w=8.5;
h=w*rho;

for k=Nl:-1:1
    figure
    set(gcf,'Color','w')
    v=get(gcf,'Position');
    ratio=v(4)/v(3);
    w=6.83*2;
    h=w*ratio;
    hold on
    for i=2:Nd
        d=domvals(i);
        if d==1
            Hfight=0;
        else
            Hfight=-(d*log(d)+(1-d)*log(1-d));
        end
        twop=reshape(twooutput(k,i,:,1),1,[]);
        Hsig_two=-(twop.*log(twop)+(1-twop).*log(1-twop));
        plot(Hsig_two,reshape(twooutput(k,i,:,2),1,[]),'Color',(1-.1*i)*[0 0 1],'LineWidth',lw)
        onep=reshape(oneoutput(k,i,1:40,1),1,[]);
        Hsig_one=-(onep.*log(onep)+(1-onep).*log(1-onep));
        plot(Hsig_one,reshape(oneoutput(k,i,1:40,2),1,[]),'Color',(1-.1*i)*[1 0 0],'LineWidth',lw)
    end
    xlabel('Probability of correct output','FontSize',textfontsz)
    ylabel('Expected time to output','FontSize',textfontsz)
%     leg=legend('Two dimensions','One dimension');
%     set(leg,'Position',[.6 .8 .1 .05])
end


%%
x=.5:.01:1;
Nx=length(x);
v=zeros(Nl,Nd-1,Nx);

for k=1:Nl
    for i=2:Nd
        v(k,i,:)=interp1(reshape(oneoutput(k,i,1:40,1),1,[]),reshape(oneoutput(k,i,1:40,2),1,[]),.5:.01:1);
    end
end

%%
l=lvals(2);
d=domvals(end-1);

threshvals2=.1:.25:2;
Nt2=length(threshvals2);

twomat=zeros(Nt2,Nt2,2);

for i=1:Nt2
    T1=threshvals2(i);
    for j=1:Nt2
        T2=threshvals2(j);
        [x,y,~]=solve_pde_2d(-T1,-T2,d,5,5,deltax2,deltax2,l,b);
        twomat(i,j,:)=[x y];
    end
end

symm=diag(twomat(:,:,1));
[L,~]=contour(threshvals2,threshvals2,transpose(twomat(:,:,1)),.5:.05:1,'LineWidth',3);
L=dividecontours(L);
numlines=length(L);

M=cell(1,numlines);

for i=1:numlines
    go=size(L{i},2);
    M{i}=zeros(go,4);
    for j=1:go
        T1=L{i}(1,j);
        T2=L{i}(2,j);
        [x,y,~]=solve_pde_2d(-T1,-T2,d,5,5,deltax2,deltax2,l,b);
        M{i}(j,:)=[T1,T2,x,y];
    end
end

for i=1:4
    [~,v]=sort(M{i}(:,2));
    M{i}=M{i}(v,:);
end

%%
contour(threshvals2,threshvals2,transpose(twomat(:,:,1)),0:.05:1,'LineWidth',3);
hold on
plot(threshvals2,threshvals2,'k','LineWidth',3)
[half,~]=contour(threshvals2,threshvals2,transpose(twomat(:,:,1)),.5);
half=dividecontours(half);
half=half{1};
plot(half(1,:),half(2,:),'k','LineWidth',3)