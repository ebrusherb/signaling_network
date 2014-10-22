function logprobmat = nummat2logprobmat(nummat)
[Nb1,Nb2]=size(nummat);
tot=sum(sum(nummat));
logprobmat=zeros(Nb1,Nb2);
for i=1:Nb1
    px=sum(nummat(i,:))/(tot);
    for j=1:Nb2
        py=sum(nummat(:,j))/(tot);
        pxy=nummat(i,j)/(tot);
        if pxy~=0
            logprobmat(i,j)=pxy*log(pxy/px/py);
        end
    end
end