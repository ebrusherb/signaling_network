function [ q, g, h, r ] = boundaryFilePiezoActuator( p, e, V )
%BOUNDARYFILEPIEZOACTUATOR Boundary conditions for piezoelectric actuator example
%   [ q, g, h, r ] = BOUNDARYFILEPIEZOACTUATOR( p, e, V ) returns the
%   Neumann BC (q, g) and Dirichlet BC (h, r) matrices for the
%   piezoelectric actuator example example.
%   p is the point matrix returned from INITMESH
%   e is the edge matrix returned from INITMESH
%   V is the voltage applied to the top layer of the beam

% Copyright 2012 The MathWorks, Inc.

N = 3;
ne = size(e,2);
q = zeros(N^2, ne);
g = zeros(N, ne);
h = zeros(N^2, 2*ne);
r = zeros(N, 2*ne);
voltage = V;
for i=1:ne
  geomEdgeID = e(5,i);
  if(geomEdgeID == 1)
    % top geometry edge
    % set the voltage to V at both vertices on all element edges
    % on this geometry edge 
    h(9,i) = 1;
    h(9,i+ne) = 1;
    r(3,i) = voltage;
    r(3,i+ne) = voltage;
  elseif(geomEdgeID == 2)
    % bottom geometry edge
    % set the voltage to zero at both vertices on all element edges
    % on this geometry edge
    % (the entries in r have already been set to zero above)
    h(9,i) = 1;
    h(9,i+ne) = 1;
  elseif(geomEdgeID == 6 || geomEdgeID == 7)
    % left geometry edge
    % set the u and v displacements to zero at both vertices on all
    % element edges on this geometry edge
    % (the entries in r have already been set to zero above)
    h(1,i) = 1;
    h(1,i+ne) = 1;
    h(5,i) = 1;
    h(5,i+ne) = 1;
  end
end

end