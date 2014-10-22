%%
axislw=1.25;
lw=3;
labfontsz=10 ;
textfontsz=12;
alpha=.7;
ticklength=[.02,.0250];
markersz=5;
markercol='white';
fontname='Times New Roman';
cornerlabels={'ER','SP','DT'};

%%
its=1000;
N=20;

cstep=.1;
vec=0:cstep:1;
l=length(vec);

% X=repmat(0:cstep:1,1,11);
% Y=reshape(repmat(0:cstep:1,11,1),1,[]);

X=[];
Y=[];
for i=0:(l-1)
    X=[X, vec(1:(end-i))];
    Y=[Y,vec(i+1)*ones(1,l-i)];
end
Z=1-X-Y;

transformed=transform([X',Z',Y']);

%%
l=length(X);
expowervecs=zeros(N,7,l);
skewvals=zeros(7,l);
mutvals=zeros(7 ,l);
accuracy=zeros(1,l);
time=zeros(1,l);

for i=1:l
    c1=X(i);
    c2=Y(i);
    c3=Z(i);
    t=group_props_parallelized(c1,c2,c3,its,N,'unif',Nd,Nt,twomat,domvals,threshvals);
    expowervecs(:,:,i)=t{1};
    mutvals(:,i)=t{2};
    skewvals(:,i)=t{3};
    accuracy(i)=t{4};
    time(i)=t{5};
end

%%
filename='triangle_heatmap_output.mat';
save(filename,'skewvals','mutvals','accuracy','time');
%%
scaledskewvals=skewvals-(repmat(min(skewvals')',1,l));
scaledskewvals=scaledskewvals./repmat(max(scaledskewvals')',1,l);

% scaledmutvals=mutvals-(repmat(min(mutvals')',1,l));
% scaledmutvals=scaledmutvals./repmat(max(scaledmutvals')',1,l);
scaledmutvals=mutvals-min(min(mutvals));
scaledmutvals=scaledmutvals/max(max(scaledmutvals));

[~,order]=sort(mutvals(1:6,:));
bestfun=order(6,:);
scaledbestfun=bestfun-min(min(bestfun));
scaledbestfun=scaledbestfun/max(max(scaledbestfun));

[~,order]=sort(skewvals(1:6,:));
mostskewed=order(6,:);
scaledmostskewed=mostskewed-min(min(mostskewed));
scaledmostskewed=scaledmostskewed/max(max(scaledmostskewed));

scaledaccuracy=accuracy-min(accuracy);
scaledaccuracy=scaledaccuracy/max(scaledaccuracy);

scaledtime=time-min(time);
scaledtime=scaledtime/max(scaledtime);

mutvalstopc=mutvals(1:6,:)-repmat(mean(mutvals(1:6,:)),6,1);
mutvalstopc=mutvalstopc./repmat(power(var(mutvalstopc),.5),6,1);

[coeff,score,latent]=princomp(mutvalstopc');
coeff=coeff(:,1:2);
score=score(:,1:2);
score=score';
scaledscore=score-repmat(min(score')',1,l);
scaledscore=scaledscore./repmat(max(score')',1,l);
%%
for i=1:6
numfun=i;
mainlabel='Mutual information';
triimage(scaledmutvals(numfun,:),mutvals(numfun,:),transformed,cornerlabels,mainlabel)
end
% filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','mutinfo_heatmap','.pdf');
% print(filename,'-dpdf','-r300');

%%
for i=6:6
numfun=i;
mainlabel='Skewness';
triimage(scaledskewvals(numfun,:),skewvals(numfun,:),transformed,cornerlabels,mainlabel)
end

% filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','skewness_heatmap','.pdf');
% print(filename,'-dpdf','-r300');

%%
mainlabel='Accuracy';
triimage(scaledaccuracy,accuracy,transformed,cornerlabels,mainlabel)

filename=strcat('/Users/eleanorbrush/Dropbox/signaling_network/','accuracy_heatmap','.pdf');
% print(filename,'-dpdf','-r300');
%%
mainlabel='Best power formalism';
triimage(scaledbestfun,bestfun,transformed,cornerlabels,mainlabel,'off');

%%
mainlabel='Most skewed power formalism';
triimage(scaledmostskewed,mostskewed,transformed,cornerlabels,mainlabel);

%%
triimage(scaledscore(2,:),score(2,:),transformed,cornerlabels,mainlabel);
