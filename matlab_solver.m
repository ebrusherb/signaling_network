%% geometry
maxX=1;
maxY=1;
Tx=1;
Ty=1;
width = 1; height = 1;
% define  the square by giving the 4 x-locations followed by the 4
% y-locations of the corners.
gdmTrans = [3 4 -Tx maxX maxX -Tx -Ty -Ty maxY maxY];
g = decsg(gdmTrans', 'S1', ('S1')');
% Create the triangular mesh on the square with approximately
% ten elements in each direction.
hmax = .01; % element size
[p, e, t] = initmesh(g, 'Hmax', hmax);

%%
%%parameters

l=1;
b=1;
d=.7;
a=b*(2*d-1);
c=b^2;

%%
A = @(p,t,u,time) acoeffunction(p,t,l,a);
C = @(p,t,u,time) ccoeffunction(p,t,c);
F = [0;0;0];
B = @(p,e,u,time) boundaryFile(p,e);

u = assempde(B, p, e, t, C, A, F);
n = size(p,2); % number of finite element nodes
uu = reshape(u, n, []);

pdeplot(p, e, t, 'xydata', uu(:,1), 'contour', 'on');