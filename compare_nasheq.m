comp=zeros(Nc,Nd,2,2);

for i=1:Nc
    for j=1:Nd
        n=nasheq_ind{i,j};
        if isempty(n)
            n=[0; 0];
        end
        n2=nasheq_ind2{i,j};
        if isempty(n2)
            n2=[0;0];
        end
        comp(i,j,1,:)=n;
        comp(i,j,2,:)=n2;
    end
end
%%
j=2;
plot(col(comp(:,j,1,2)),'-o')
hold on
plot(col(comp(:,j,2,2)),'r')
hold off