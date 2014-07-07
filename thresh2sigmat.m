function [sigmat, timemat]=thresh2sigmat(Tvals,twomat,dmat,l)
reshapedTvals=reshape(Tvals,1,[]);
N=size(reshapedTvals,2);
probmat=zeros(N,N);
timemat=zeros(N,N);
Nd=size(twomat,2);
Nt=size(twomat,3);

%size(twomat)=[Nl,Nd,Nt,Nt,2];

for i=1:N
    for j=[1:(i-1),(i+1):N]
        if dmat(i,j)<=Nd-1
            probmat(i,j)=1-twomat(l,Nd-dmat(i,j)+1,reshapedTvals(j),reshapedTvals(i),1);
            timemat(i,j)=twomat(l,Nd-dmat(i,j)+1,reshapedTvals(j),reshapedTvals(i),2);
        else
            probmat(i,j)=twomat(l,dmat(i,j)-Nd+1,reshapedTvals(i),reshapedTvals(j),1);
            timemat(i,j)=twomat(l,dmat(i,j)-Nd+1,reshapedTvals(i),reshapedTvals(j),2);
        end
    end
end
probmat(1:N+1:end)=1;

sigmat=zeros(N,N);
    
    for i=1:N
        for j=(i+1):N
            draw=rand;
            if draw <=probmat(i,j)
                sigmat(i,j)=1;
            else sigmat(j,i)=1;
            end
        end
    end