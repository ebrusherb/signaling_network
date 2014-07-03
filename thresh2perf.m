function perfvec = thresh2perf(cvec,Tvals,twomat,dmat,l)
% Tvals=reshape(Tvals,[],1);
% c1=cvec(1);
% c2=cvec(2);
% c3=cvec(3);
% Nd=size(twomat,2);
% Nt=size(twomat,3);
% N=size(Tvals,1);
% perf=zeros(2,Nd,Nt,Nt); %individual, leak, dominance, thresholds
% 
% perf(1,2:Nd,:,:)=c1*(1-twomat(l,2:Nd,:,:,1))+c2*(twomat(l,2:Nd,:,:,2))+c3*(1-twomat(l,2:Nd,:,:,1));
% perf(2,2:Nd,:,:)=c1*(1-twomat(l,2:Nd,:,:,1))+c2*(twomat(l,2:Nd,:,:,2))+c3*(twomat(l,2:Nd,:,:,1));
% 
% perf(1,1,:,:)=c2*(twomat(l,1,:,:,2))+c3*(1-twomat(l,1,:,:,1));
% perf(2,1,:,:)=c2*(twomat(l,1,:,:,2))+c3*(twomat(l,1,:,:,1));
%         
% perfvec=zeros(1,N);
% for i=1:N
%     opp_thresh=Tvals([1:(i-1),(i+1):N]);
%     ds=dmat(i,[1:(i-1),(i+1):N]);
%     perfsum=0;
%     for q=1:(N-1)
%         if ds(q)<=Nd-1
%             perfsum=perfsum+perf(2,Nd-ds(q)+1,opp_thresh(q),Tvals(i));
%         else 
%             perfsum=perfsum+perf(1,ds(q)-Nd+1,Tvals(i),opp_thresh(q));
%         end
%     end 
%     perfvec(i)=perfsum;
% end
% perfvec=reshape(perfvec,[],1);

Tvals=reshape(Tvals,1,[]);
dmat=reshape(dmat,max(size(dmat)),[]);
c1=cvec(1);
c2=cvec(2);
c3=cvec(3);
[probmat,timemat,accmat]=thresh2mat(Tvals,twomat,dmat,l);
perfvec=c1*(1-accmat)+c2*timemat+c3*(1-probmat);
perfvec=sum(perfvec,2);
perfvec=reshape(perfvec,[],1);