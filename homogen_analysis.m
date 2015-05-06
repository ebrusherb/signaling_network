meanvec=mean(reshape(finalthreshmat,20000,66),1);
varvec=var(reshape(finalthreshmat,20000,66),1);
cornerlabels={'Error rate','Personal preference','Decision time'};

%% coordinate set up
N=20;
its=1000;
cstep=.1;
vec=0:cstep:1;
Nc=length(vec);
endvals=zeros(1,Nc);

% X=repmat(0:cstep:1,1,11);
% Y=reshape(repmat(0:cstep:1,11,1),1,[]);

X=[];
Y=[];
layers=cell(Nc,1);
for i=0:(Nc-1)
    layers{i+1}=length(X)+(1:(Nc-i));
    endvals(i+1)=layers{i+1}(end);
    X=[X, vec(1:(Nc-i))];
    Y=[Y,vec(i+1)*ones(1,Nc-i)];   
end
Z=1-X-Y;

transformed=transform([X',Z',Y']);
Nc=length(X);

indices=1:Nt;
leak=2;

mainlabel='';

triimage(varvec,transformed,cornerlabels,mainlabel,'off');

%%
k=layers{2}(5);
figure
hold on

for i=1:1000
    plot(sort(finalthreshmat(:,i,k)))
end