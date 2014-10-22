function powerindices=vec2indices(powervec,Nb2)
%     powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    %     powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
q=quantile(powervec,0:1/Nb2:1);
diffvec=diff(q);
eps=min(diffvec(diffvec>=1e-5))/10;
q(1)=q(1)-eps;q(end)=q(end)+eps;
powerbins=q;
[~,powerindices]=histc(powervec,powerbins);