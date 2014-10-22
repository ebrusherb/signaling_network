Lxvals=-1:.1:-.1;
Lyvals=-1:.1:-.1;
dvals=0:.2:1;
NLx=length(Lxvals);
NLy=length(Lyvals);
Nd=length(dvals);
Lxopt=zeros(NLy,Nd,2);

l=1;
b=1;

for i=1:NLy
    Ly=Lyvals(i);
    for j=1:Nd
    d=dvals(j);
    perfvals=zeros(1,NLx);
        for k=1:NLx
            Lx=Lxvals(k);
            %Lx,Ly,d,Ux,Uy,deltax,deltay,l,b
%             [x,y,z]=solve_pde(Lx,Ly,d,5,5,.05,.05,l,b);
%             perfvals(k)=z;

%             [x,y,~]=solve_pde(Lx,Ly,d,5,5,.05,.05,l,b);
%             if d<.5
%                 perfvals(k)=(1-x)/y;
%             else perfvals(k)=x/y;
%             end

            [x,y,z]=solve_pde(Lx,Ly,d,5,5,.05,.05,l,b);
            perfvals(k)=-y+x*(1-y);
        end
        [m,w]=max(perfvals);
        Lxopt(i,j,:)=[m,Lxvals(w)];
    end
end

%%
imagesc(Lyvals,dvals,Lxopt(:,:,2))
set(gca,'ydir','normal')
xlabel('Other Animal Threshold')
ylabel('Focal Animal Dominance')
colorbar

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','Lx_opt_strategy','.eps');
print(filename,'-depsc','-r300');

imagesc(Lyvals,dvals,Lxopt(:,:,1))
set(gca,'ydir','normal')
xlabel('Other Animal Threshold')
ylabel('Focal Animal Dominance')
colorbar

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','Lx_opt_performance','.eps');
print(filename,'-depsc','-r300');

%%
l=length(dvals);
figure
hold on

for i=0:(l-1)
    pp1=interp1(Lyvals,Lxopt(:,1+i,2),'linear','pp');
    pp2=interp1(Lyvals,Lxopt(:,l-i,2),'linear','pp');
    Tystar=Lyvals(abs(ppval(pp2,ppval(pp1,Lyvals))-Lyvals)<0.001);
    Txstar=ppval(pp1,Tystar);
    plot(dvals(i+1),-Txstar,'or')
    plot(dvals(i+1),-Tystar,'ob')
end