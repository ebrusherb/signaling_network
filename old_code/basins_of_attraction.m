global f eps e t b
f=0.1;
eps=2;
e=.5;
b=-2*eps*e*(eps+f/e)/(f/e-eps);
t=(2*eps*e*(eps+f/e)+b*(f/e-eps))/(2*eps*e-b);

xvals = -t-eps-1:.1:f/e+rate(-f/e)/e+3;
yvals = -t-eps-1:.1:-f/e+rate(f/e)/e+3;

stepsx = min(xvals):2:max(xvals);
stepsy = min(yvals):2:max(yvals);


xgradient=zeros([length(yvals),length(xvals),2]);
ygradient=zeros([length(yvals),length(xvals),2]);
zgradient=zeros([length(yvals),length(xvals),2]);
for i=1:length(xvals)
    for j=1:length(yvals)
        xgradient(j,i,:) = -e*xvals(i)+f+rate(yvals(j));
        ygradient(j,i,:) = -e*yvals(j)-f+rate(xvals(i));
    end
end

[X,Y,Z]=meshgrid(xvals,yvals,1:2);
[sx,sy,sz]=meshgrid(stepsx,stepsy,1);
figure
h=streamslice(X,Y,Z,xgradient,ygradient,zgradient,[],[],[1]);
set(h,'Color','black')
set(gca,'box','off')
axis tight
hold on

plot (f/e,-f/e+rate(f/e)/e,'ro')
plot (f/e+rate(-f/e)/e,-f/e,'ro')
xhat = b(eps-t)/(2*eps*e+b)-2*eps*f/(b-2*eps*e);
yhat = b(eps-t)/(2*eps*e+b)+2*eps*f/(b-2*eps*e);
plot (xhat,yhat,'ro')

plotx = -t-eps-1:.01:timefunx(-t-eps);
xhat = timefunx(-t-eps);
yhat = -t-eps;
ploty = (yhat-(b-f)/e)/(xhat-(b+f)/e)*(plotx-(b+f)/e)+(b-f)/e;
plot(plotx,ploty,'r')

ploty = -t-eps:.1:-t-eps+4*eps*f/(b-2*eps*e);
plotx = timefunx(ploty);
plot(plotx,ploty,'r')

plotx = -t-eps:.1:-t+eps-4*eps*f/(b-2*eps*e);
ploty = plotx+4*eps*f/(b-2*eps*e);
plot(plotx,ploty,'r')

plotx = -t+eps-4*eps*f/(b-2*eps*e):.1:-t+eps;
ploty = timefuny(plotx);
plot(plotx,ploty,'r')

plotx = -t+eps:.1:max(xvals)+1;
xhat = -t+eps;
yhat = timefuny(-t+eps);
ploty = (yhat+f/e)/(xhat-f/e)*(plotx-f/e)-f/e;
plot(plotx,ploty,'r')

plot(min(xvals):(max(xvals)+1),0*(min(xvals):(max(xvals)+1)),'color',[0 0 0])
plot(0*(min(yvals):(max(yvals)+1)),(min(yvals):(max(yvals)+1)),'color',[0 0 0])

plotx=min(xvals):.1:max(xvals);
ploty=-f/e+rate(plotx)/e;
plot(plotx,ploty,'b')

ploty=min(yvals):.1:max(yvals);
plotx=f/e+rate(ploty)/e;
plot(plotx,ploty,'b')
