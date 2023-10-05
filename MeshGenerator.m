function [FEM,mesh,model,line] = MeshGenerator(x,y)
%MeshGenerator.m - Generates the FEM matrices, a triangular mesh, the
%PDE model and the edges forming the boundary of a polygon. It uses the
%PDEmodeler App from Matlab
%
%INPUTS:
%   -x: X-coordinates of the vertices of the polygon.
%   -y: Y-coordinates of the vertices of the polygon.

shp = polyshape(x,y); %% Create the poligon from the data points
tr = triangulation(shp);
model = createpde;
tnodes = tr.Points';
telements = tr.ConnectivityList';
geometryFromMesh(model,tnodes,telements);
specifyCoefficients(model,'m',0,'d',1,'c',1,'a',0,'f',0); %%Parameters for the equation: m (d^2u/dt^2) + d du/dt + nabla · (c nabla u) + a u= f.
mesh = generateMesh(model,'GeometricOrder','linear'); %%Linear Geometric order is used in order to have a less dense mesh in which calculations run faster.
FEM = assembleFEMatrices(model);
line=cell(1,length(x));
line{length(x)} = unique(findNodes(mesh,'nearest',[linspace(x(end),x(1))',linspace(y(end),y(1))']'));
for j=2:length(x)
    line{j-1}= unique(findNodes(mesh,'nearest',[linspace(x(j-1),x(j))',linspace(y(j-1),y(j))']'));
end
