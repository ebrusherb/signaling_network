
expowervecs_lognorm=zeros(N,7,Nc);
skewvals_lognorm=zeros(7,Nc);
mutvals_lognorm=zeros(7 ,Nc);
accuracy_lognorm=zeros(1,Nc);
time_lognorm=zeros(1,Nc);
avgthresh_lognorm=zeros(N,Nc);
timemats_lognorm=zeros(N,N,Nc);
probmats_lognorm=zeros(N,N,Nc);

bottomquart=16:20;
bottomaccuracy_lognorm=zeros(1,Nc);

topquart=1:5;
topaccuracy_lognorm=zeros(1,Nc);

for i=1:Nc
    c1=X(i);
    c2=Y(i);
    c3=Z(i);
    t=group_props_parallelized(c1,c2,c3,its,N,'lognorm',Nd,Nt,twomat,domvals,threshvals);
    expowervecs_lognorm(:,:,i)=t{1};
    mutvals_lognorm(:,i)=t{2};
    skewvals_lognorm(:,i)=t{3};
    accuracy_lognorm(i)=t{4};
    time_lognorm(i)=t{5};
    avgthresh_lognorm(:,i)=t{7};
    timemats_lognorm(:,:,i)=t{8};
    probmats_lognorm(:,:,i)=t{9};
    totalsigmat=t{10};
    bottomaccuracy_lognorm(i)=sum(sum(triu(totalsigmat(bottomquart,bottomquart),1)))/sum(sum(totalsigmat(bottomquart,bottomquart)));
    topaccuracy_lognorm(i)=sum(sum(triu(totalsigmat(topquart,topquart),1)))/sum(sum(totalsigmat(topquart,topquart)));
end

%%
scaledskewvals_lognorm=skewvals_lognorm-(repmat(min(skewvals_lognorm,[],2),1,Nc));
scaledskewvals_lognorm=scaledskewvals_lognorm./repmat(max(scaledskewvals_lognorm,[],2),1,Nc);

scaledmutvals_lognorm=mutvals_lognorm-(repmat(min(mutvals_lognorm,[],2),1,Nc));
scaledmutvals_lognorm=scaledmutvals_lognorm./repmat(max(scaledmutvals_lognorm,[],2),1,Nc);

[~,order]=sort(mutvals_lognorm(1:6,:));
bestfun_lognorm=order(6,:);
bestfun_lognorm(bestfun_lognorm==3)=2;
bestfun_lognorm(bestfun_lognorm==5)=2;
bestfun_lognorm(bestfun_lognorm==6)=3;
scaledbestfun_lognorm=bestfun_lognorm-min(min(bestfun_lognorm));
scaledbestfun_lognorm=scaledbestfun_lognorm/max(max(scaledbestfun_lognorm));

scaledaccuracy_lognorm=accuracy-min(accuracy_lognorm);
scaledaccuracy_lognorm=scaledaccuracy_lognorm/max(scaledaccuracy_lognorm);

scaledbottomacc_lognorm=bottomaccuracy_lognorm-min(bottomaccuracy_lognorm);
scaledbottomacc_lognorm=scaledbottomacc_lognorm/max(scaledbottomacc_lognorm);

scaledtopacc_lognorm=topaccuracy_lognorm-min(topaccuracy_lognorm);
scaledtopacc_lognorm=scaledtopacc_lognorm/max(scaledtopacc_lognorm);

