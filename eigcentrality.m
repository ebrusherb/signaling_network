function v=eigcentrality(M,weight)
N=size(M,2);
totalsum=sum(sum(M));
if totalsum==0
    v=zeros(N,1);
else
    stochM=M/totalsum;
    stochM(1:N+1:end)=1-sum(stochM,1);

    unifmat=1/N*ones(N);

    tofindeigenvectorsof=(1-weight)*stochM+weight*unifmat;
    [vecs,vals]=eig(tofindeigenvectorsof);
    f=find(abs(real(diag(vals))-1)<1e-5,1,'first');
    v=vecs(:,f);
    v=abs(v);
end


