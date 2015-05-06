function triimage_discrete(scaledvals,vals,coords,cornerlabels,labels,fourcols)
n=length(unique(vals));
indices=1:n;

red=evalin('base','red');
yellow=evalin('base','yellow');
green=evalin('base','green');
blue=evalin('base','blue');
textfontsz=evalin('base','textfontsz');
labfontsz=evalin('base','labfontsz');
fontname=evalin('base','fontname');
powerfuns=evalin('base','powerfuns');
powerfuns{strcmp(powerfuns,'Delta')}='\Delta';
layers=evalin('base','layers');

% c=colormap(jet(n));
c=fourcols;

l=size(coords,1);
cstep=(coords(1,1)-coords(2,1))/2;

xoffset=.25;
yoffset=.2;

% v=get(axes,'position');
% set(axes,'position',[0.02 v(2:end)],'off');
hold on
axis off
axis equal


for i=1:l
    p=coords(i,:);
    valnow=scaledvals(i);
    valcol=c(valnow==unique(scaledvals),:);
    fill([p(1),p(1)+cstep,p(1)+cstep,p(1),p(1)-cstep,p(1)-cstep],[p(2)-2/sqrt(3)*cstep,p(2)-1/sqrt(3)*cstep,p(2)+1/sqrt(3)*cstep,p(2)+2/sqrt(3)*cstep,p(2)+1/sqrt(3)*cstep,p(2)-1/sqrt(3)*cstep],valcol,'EdgeColor','none');
end
xl=get(gca,'xlim');
set(gca,'xlim',[-.5 xl(2)]) 

text(0-xoffset,-yoffset,cellstr(cornerlabels{1}),'Color','k','HorizontalAlignment','left','FontSize',textfontsz,'FontName',fontname)
text(2+xoffset,-yoffset,cellstr(cornerlabels{2}),'Color','k','HorizontalAlignment','right','FontSize',textfontsz,'FontName',fontname)
text(1,sqrt(3)+yoffset,cellstr(cornerlabels{3}),'HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz,'FontName',fontname)

[x1,y1]=data2norm(coords(layers{10}(1),1),coords(layers{10}(1),2));
annotation('textarrow',[x1+.03 x1-.05], [y1 y1],'String',['  ' labels{1}],'HeadStyle','none','FontSize',labfontsz,'FontName',fontname)
[x2,y2]=data2norm(coords(layers{1}(1),1),coords(layers{1}(1),2));
annotation('line',[x2-.15 x1+.07],[y2 y1-.03])
[x3,y3]=data2norm(coords(layers{7}(end),1),coords(layers{7}(end),2));
annotation('textarrow',[x3-.01 x3+.04], [y3 y3],'String',[labels{2} '  ' ],'HeadStyle','none','FontSize',labfontsz,'FontName',fontname)
[x3,y3]=data2norm(coords(layers{5}(end),1),coords(layers{5}(end),2));
annotation('textarrow',[x3+.03 x3+.07], [y3 y3],'String',[labels{3} '  ' ],'HeadStyle','none','FontSize',labfontsz,'FontName',fontname)
[x3,y3]=data2norm(coords(layers{3}(end),1),coords(layers{3}(end),2));
annotation('textarrow',[x3+.04 x3+.095], [y3 y3],'String',[labels{4} '  ' ],'HeadStyle','none','FontSize',labfontsz,'FontName',fontname)
if size(labels,2)>4
    annotation('textbox',[x3-.04 y3-.04 .01 .04],'String',labels{5},'FitBoxToText','on','FontSize',labfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
end