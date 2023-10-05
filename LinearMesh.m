function [M,K]=LinearMesh(n)
%LinearMesh.m - Generates the FE matrices of a 1-dimensional mesh to model diffusion on a line.
%
%INPUTS:
%   -n: Number of elements of the linear mesh
%
%OUTPUTS:
%   -M: Mass matrix of the mesh.
%   -K: Stiffness matrix of the mesh.

L = 1;              % length of interval
% n = 19;             % number of elements
dt=0.01;
h = L/n;            % element length
D = 0;            % diffusion coefficient
f = @(x) 0;         % external force function

% Construct matrices
M = zeros(n+1,n+1);
K = zeros(n+1,n+1);
for i = 1:n
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1) + h/6*[2 1;1 2];
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + 1/h*[1 -1;-1 1] + dt*D/h*[1 -1;-1 1];
end

% Modify the stiffness matrix for Neumann boundary conditions
K(1,1) = K(1,1) + dt*D/h;
K(1,2) = K(1,2) - dt*D/h;
K(n+1,n+1) = K(n+1,n+1) + dt*D/h;
K(n+1,n) = K(n+1,n) - dt*D/h;

end
