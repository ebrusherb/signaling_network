load /Users/eleanorbrush/Desktop/groupsize.mat

coloffset=4;
mycols=cbrewer('seq','YlOrRd',Nlays+coloffset);
cvec=0:(1/(Nlays-1)):1;
leglabs=cell(1,Nlays);
for i=1:Nlays
    leglabs{i}=num2str(cvec(i));
end

subplot(1,3,1)
hold on
p=zeros(1,Nlays);
for j=1:Nlays
    p(j)=plot(sizes,accuracy(:,j,1),'Color',mycols(j+coloffset,:));
    plot(sizes,accuracy(:,j,2),'--','Color',mycols(j+coloffset,:))
end
leg=legend(p,leglabs);

subplot(1,3,2)
hold on
p=zeros(1,Nlays);
for j=1:Nlays
    p(j)=plot(sizes,topaccuracy(:,j,1),'Color',mycols(j+coloffset,:));
    plot(sizes,topaccuracy(:,j,2),'--','Color',mycols(j+coloffset,:))
end
leg=legend(p,leglabs);

subplot(1,3,3)
hold on
p=zeros(1,Nlays);
for j=1:Nlays
    p(j)=plot(sizes,bottomaccuracy(:,j,1),'Color',mycols(j+coloffset,:));
    plot(sizes,bottomaccuracy(:,j,2),'--','Color',mycols(j+coloffset,:))
end
leg=legend(p,leglabs);
