function [ qmatrix, gmatrix, hmatrix, rmatrix ] = boundaryFile( p, e )

N = 3; % Set N = the number of equations
ne = size(e,2); % number of edges
qmatrix = zeros(N^2,ne);
gmatrix = zeros(N,ne);
hmatrix = zeros(N^2,2*ne);
rmatrix = zeros(N,2*ne);

for k = 1:ne
    hmatrix([1 5 9],k)=1;
    hmatrix([1 5 9],k+ne)=1;
    switch e(5,k)
        case {1,2}
            % Fill in hmatrix,rmatrix or qmatrix,gmatrix
            rmatrix(1,k)=1;
            rmatrix(1,k+ne)=1;
    end
end
