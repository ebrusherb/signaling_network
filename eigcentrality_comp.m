function [v, v2]=eigcentrality_comp(mat,weight)

N=size(mat,2);
totalsum=sum(sum(mat));
if totalsum==0
    v=zeros(N,1);
else
    stochM=mat/totalsum;
    stochM(1:N+1:end)=1-sum(stochM,1);    

    unifmat=1/N*ones(N);

    tofindeigenvectorsof=(1-weight)*stochM+weight*unifmat;
    [vecs,vals]=eig(tofindeigenvectorsof);
    f=find(abs(real(diag(vals))-1)<1e-5,1,'first');
    v=vecs(:,f);
    v=abs(v);
    
    out=sum(mat,1);
    out(out==0)=1;
    stochM=mat./repmat(out,N,1);
    
    unifmat=1/N*ones(N);

    tofindeigenvectorsof=(1-weight)*stochM+weight*unifmat;
    [vecs,vals]=eig(tofindeigenvectorsof);
    f=find(abs(real(diag(vals))-1)<1e-5,1,'first');
    v2=vecs(:,f);
    v2=abs(v2);
end