function C = ccoeffunction(p,t,c)

N = 3; % Number of equations
% Triangle point indices
it1=t(1,:);
it2=t(2,:);
it3=t(3,:);

% Find centroids of triangles
xpts=(p(1,it1)+p(1,it2)+p(1,it3))/3;
ypts=(p(2,it1)+p(2,it2)+p(2,it3))/3;

nt = size(t,2); % Number of columns
C = zeros(N^2*4,nt); % Allocate f

% Now the particular functional form of f
C(2,:) = -ypts;
C(3,:) = ypts;
C(6,:)=xpts;
C(7,:)=-xpts;
C(22,:)=c*(xpts-ypts);
C(23,:)=c*(ypts-xpts);
C(34,:)=c*(xpts-ypts);
C(35,:)=c*(ypts-xpts);
