%% Default Simplices
function [mesh,Mass,sink,source] = DefaultSimplex(A,Types)
%DefaultSimplex.m - Generates the metaplex structure of a simplicial
%complex. Uses the functions MeshGenerator.m and LinearMesh.m
%
%INPUTS:
%   -A: Graph of the simplicial complex, where the edges and triangles
%   are also considered as nodes.
%   -Types: List of the type of object that it is each node.
%  Types: 1=Node, 2=line, 3=triangle with sinks on the vertex, 4= triangle
%  with sinks on the edges.
%  (Example of a triangle simplex:
%  A=[0 1 1 1 0 0 0;1 0 0 0 1 1 0;1 0 0 0 0 1 1;1 0 0 0 1 0 1;0 1 0 1 0 0 0;0 1 1 0 0 0 0;0 0 1 1 0 0 0],
%  Types=[4,2,2,2,1,1,1])
%
%OUTPUTS:
%   -mesh: One dimensional cell of matrices representing the diffusion inside each
%   node in the metaplex.
%   -Mass: One dimensional cell of matrices representing the structure of each node in
%   the metaplex.
%   -Sink/Source: Two dimensional cells containing the matrices that
%   represent the sink and sources from one node to another.


N=length(A);

if length(Types)~=N
    error('Need to assign a type to each node')
end

mesh=cell(1,N);
Mass=cell(1,N);
sink=cell(N,N);
source=cell(N,N);

for i=1:N
    switch Types(i)
        case 1
            Mass{i} = 1;
            mesh{i} = 0;
            for j=1:N
                if A(i,j)==1
                    sink{i,j} = 1;
                    source{i,j} = 1;
                end
            end
        case 2
            linelen=20; %% This is the size of the lines in the triangles by default
            [M,K]=LinearMesh(linelen-1);
            Mass{i} = M/sum(M,'All'); %% Mass of the nodes over the total length of the line
            mesh{i} = K;
            c=0;
            for j=1:N
                if A(i,j)==1
                    if Types(j)==4
                        sink{i,j} = [1:linelen];
                        source{i,j} = [1:linelen];
                    else
                        sink{i,j} = 1+(linelen-1)*mod(c,2);
                        source{i,j} = 1+(linelen-1)*mod(c,2);
                    end
                    c=c+1;
                end
            end
        case 3
            x = [0, 0.5, 1];
            y = [0, sqrt(3)/2, 0];
            [FEM,~,~] = MeshGenerator(x,y);
            mesh{i}(:,:)=FEM.K;
            mesh{i}=mesh{i}/max(max(mesh{i}));
            Mass{i}(:,:)=FEM.M; %% Mass of the nodes over the total area of the figure
            Mass{i}(:,:)=Mass{i}(:,:)/sum(Mass{i}(:,:),'All');
            c=0;
            for j=1:N
                if A(i,j)==1
                    sink{i,j} = 1+mod(c,3);
                    source{i,j} = 1+mod(c,3);
                    c=c+1;
                end
            end
        case 4
            x = [0, 0.5, 1];
            y = [0, sqrt(0.75), 0];
            [FEM,~,~,line] = MeshGenerator(x,y);
            mesh{i}(:,:)=FEM.K;
            mesh{i}=mesh{i}/max(max(mesh{i}));
            Mass{i}(:,:)=FEM.M/sum(FEM.M,'All'); %% Mass of the nodes over the total area of the figure
            c=0;
            for j=1:N
                if A(i,j)==1
                    sink{i,j} = line{1+mod(c,3)};
                    source{i,j} = line{1+mod(c,3)};
                    c=c+1;
                end
            end
    end
end
end
