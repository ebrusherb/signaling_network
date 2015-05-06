N=20;
q=N/4;
topquart=1:q;
bottomquart=fliplr(N+1-(1:q));

layerstocheck=[1:3];
Nlaystocheck=length(layerstocheck);
acc_tocheck=zeros(Nlaystocheck,1);
bottomacc_tocheck=zeros(Nlaystocheck,1);
topacc_tocheck=zeros(Nlaystocheck,1);

for l=1:Nlaystocheck
    i=layers{l}(1);
    c1=X(i);
    c2=Y(i);
    c3=Z(i);
    t=group_props_parallelized(c1,c2,c3,its,N,'unif',Nd,Nt,twomat,domvals,threshvals,leak);
    acc_tocheck(l)=t{4};
    totalsigmat=t{11};
    bottomacc_tocheck(l)=sum(sum(triu(totalsigmat(bottomquart,bottomquart),1)))/sum(sum(totalsigmat(bottomquart,bottomquart)));
    topacc_tocheck(l)=sum(sum(triu(totalsigmat(topquart,topquart),1)))/sum(sum(totalsigmat(topquart,topquart)));
end

'Finished'

%%
bottomquart=(N+1-(1:q));
a=sum(sum(triu(totalsigmat(bottomquart,bottomquart),1)));
b=sum(sum(totalsigmat(bottomquart,bottomquart)));
[a b a/b]