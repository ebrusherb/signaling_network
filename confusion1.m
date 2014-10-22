perf=zeros(2,Nd,Nt,Nt); %individual, leak, dominance, thresholds

perf(1,2:Nd,:,:)=c1*(1-twomat(l,2:Nd,:,:,1))+c2*(twomat(l,2:Nd,:,:,2))+c3*(1-twomat(l,2:Nd,:,:,1));
perf(2,2:Nd,:,:)=c1*(1-twomat(l,2:Nd,:,:,1))+c2*(twomat(l,2:Nd,:,:,2))+c3*(twomat(l,2:Nd,:,:,1));

perf(1,1,:,:)=c2*(twomat(l,1,:,:,2))+c3*(1-twomat(l,1,:,:,1));
perf(2,1,:,:)=c2*(twomat(l,1,:,:,2))+c3*(twomat(l,1,:,:,1));

perfvals=zeros(N,1);
indperfvals=zeros(N,N);
for i=1:N
    opp_thresh=Tvals([1:(i-1),(i+1):N]);
    ds=dmat(i,[1:(i-1),(i+1):N]);
    perfsum=0;
            for q=1:(N-1)
                if ds(q)<=Nd-1
                    perfsum=perfsum+perf(2,Nd-ds(q)+1,opp_thresh(q),Tvals(i)); 
                else 
                    perfsum=perfsum+perf(1,ds(q)-Nd+1,Tvals(i),opp_thresh(q));
                end
            end  
    perfvals(i)=perfsum;
    for j=[1:(i-1), (i+1):N]
        ds=dmat(i,j);
        if ds<=Nd-1
            indperfvals(i,j)=perf(2,Nd-ds+1,Tvals(j),Tvals(i)); 
        else 
            indperfvals(i,j)=perf(1,ds-Nd+1,Tvals(i),Tvals(j));
        end
    end
end

[probmat, timemat, accmat]=thresh2mat(Tvals,twomat,dmat,l);

v=c1*(1-accmat)+c2*timemat+c3*(1-probmat);
v=sum(v,2);