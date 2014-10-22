function A = acoeffunction(p,t,l,a)

N = 3; % Number of equations
% Triangle point indices
it1=t(1,:);
it2=t(2,:);
it3=t(3,:);

% Find centroids of triangles
xpts=(p(1,it1)+p(1,it2)+p(1,it3))/3;
ypts=(p(2,it1)+p(2,it2)+p(2,it3))/3;

nt = size(t,2); % Number of columns
A = zeros(N^2,nt); % Allocate a

% Now the particular functional form of f
A(4,:) = 1;
A(7,:) = 1;
A(6,:) = -l*xpts+a;
A(9,:) = -l*ypts-a;