%%
lvals=0:.1:1;
Nl=length(lvals);

domvals=.5:.05:1;
Nd=length(domvals);

threshvals=.5:.1:2;
Nt=length(threshvals);

indices=1:Nt;

b=1;
deltax2=0.025;

deltax1=0.01;
threshvals1=.1:.1:8;
Nt1=length(threshvals1);

indices1=1:Nt1;


filename=strcat('variables','.mat');
save(filename,'Nl','Nd','Nt','Nt1','lvals','domvals','threshvals','threshvals1','deltax2','deltax1');


twomat=zeros(Nl,Nd,Nt,Nt,2); %leak, dominance, thresholds, prob, time

for k=2
    l=lvals(k);
    for u=4:Nd
        d=domvals(u);
        for i=1:Nt
            T1=threshvals(i);
            for j=1:Nt
                T2=threshvals(j);
                [x,y,~]=solve_pde_2d(-T1,-T2,d,5,5,deltax2,deltax2,l,b);
                twomat(k,u,i,j,:)=[x y];
            end
        end
    end
end

% 
% onemat=zeros(Nl,Nd,Nt1,Nt1,2);
% 
% for k=1:Nl
%     l=lvals(k);
%     for u=1:Nd
%         d=domvals(u);
%         for i=1:Nt1
%             t1=threshvals1(i);
%             for j=1:Nt1
%                 t2=threshvals1(j);
%                 [x,y,~]=solve_pde_1d(-t1,d,t2,deltax1,l,b);
%                 onemat(k,u,i,j,:)=[x y];
%             end
%         end
%     end
% end


filename=strcat('/Users/eleanorbrush/Desktop/twomat_long_run','.mat');
save(filename,'twomat','lvals','domvals','threshvals','threshvals');