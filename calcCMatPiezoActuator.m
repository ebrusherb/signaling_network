function c = calcCMatPiezoActuator( p, t, c11, c12, c22, c13, c23, c33 )
%CALCCMATPIEZOACTUATOR C-matrix for piezoelectric actuator example
%   c = CALCCMATPIEZOACTUATOR( p, t, c11, c12, c22, c13, c23, c33 )
%   returns the 'c' coefficient matrix for the piezoelectric actuator
%   example given the point and element matrices along with the
%   constitutive submatrices (cij) for the PVDF material.

% Copyright 2012 The MathWorks, Inc.

numElems = size(t,2);
c=zeros(21,numElems);
%
% Although the material in both layers is PVDF, in the top layer
% the polarization direction points down (minus y direction) and in the
% bottom layer, it points up. That is, the top layer has d-coefficients
% that are the negative of those in the bottom layer.
%
% The code below examines the y-location of the centroid of each
% triangular element and assigns the correct material properties to
% element depending on whether it is in the top or bottom layer.
%
ctop = [c11(:); c12(:); c22(:); -c13(:); -c23(:); -c33(:)];
cbot = [c11(:); c12(:); c22(:);  c13(:);  c23(:); -c33(:)];
% calculate y-coordinate of triangle centers
yCenter=(p(2,t(1,:)) + p(2,t(2,:)) + p(2,t(3,:)))/3;
for i=1:numElems
  if(yCenter(i) < 0)
    c(:,i) = cbot;
  else
    c(:,i) = ctop;
  end
end

end