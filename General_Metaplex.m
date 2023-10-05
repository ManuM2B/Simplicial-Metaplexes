%% GeneralMetaplex
%This script can be used to test a diffusion process on random graphs and two triangles connected in
%which the 4th node is a metaplex connected to the 3rd and 5th node through
%a point.

clear all
% Polygon vertices of the sahpe that the nodes in the metaplex will have
x = [0, 0.5, 1];
y = [0, sqrt(3)/2, 0];

% We generate the mesh and matrices of the metaplex
[FEM,MeshModel,~] = MeshGenerator(x,y);
%
nummesh=length(MeshModel.Nodes);

%% Network creation
%Random graph
% N=2; p=1;
% A=rand(N,N);
% A=A<p;
% A=triu(A);
% A=A+A';
% A=A-diag(diag(A));

% A=[0 1 0 0; 1 0 1 0;0 1 0 1;0 0 1 0]; %Triangle graph
A = [0 1 1 0 0 0 0;1 0 1 0 0 0 0;1 1 0 1 0 0 0;0 0 1 0 1 0 0;0 0 0 1 0 1 1;0 0 0 0 1 0 1;0 0 0 0 1 1 0]; %Triangle connected with two nodes connected to two other nodes.
N=length(A);

%Degree-bias
K=diag(sum(A,2));
Ak=K^(-1)*A*K;
LapOut=diag(sum(Ak,2))-Ak;

mesh=cell(1,N);
Mass=cell(1,N);

sink=cell(N,N);
source=cell(N,N);


for i=1:N
    if i~=4
        mesh{i}(:,:)= 0;
        Mass{i}(:,:)= 1;
    end
    for j=1:N
        if A(i,j)==1
            sink{i,j} = 1;
            source{j,i} = 1;
        end
    end
end

mesh{4}(:,:)=FEM.K;
Mass{4}(:,:)=FEM.M/sum(FEM.M,'All'); %% Normalized so the total area is equal to 1

sink{4,5}=3; source{4,5}=3;
sink{5,4}=3; source{5,4}=3;
sink{4,3}=1; source{4,3}=1;
sink{3,4}=1; source{3,4}=1;


%% Metaplex structure
alpha=ones(N,N);
c=0;
for i=1:N
    s=c+1;
    c=c+size(mesh{i}(:,:),1);
    SS(s:c,s:c)=mesh{i}(:,:);
    MM(s:c,s:c)=Mass{i}(:,:);
    VV(s:c,s:c)=0;
end
a=0;
for i=1:N
    b=0;
    for j=1:N
        VV(a+source{i,j},b+sink{j,i})=VV(a+source{i,j},b+sink{j,i})-A(i,j)*alpha(i,j);
        VV(b+sink{j,i},b+sink{j,i})=VV(b+sink{j,i},b+sink{j,i})+A(i,j)*alpha(i,j);
        b=b+size(mesh{j}(:,:),1);
    end
    a=a+size(mesh{i}(:,:),1);
end

VV=VV';

%% Metaplex dynamics

nn = 10^3;    % Number of time steps used to find uht
h = 0.01;     % Time step

n = size(SS,1);
uht = zeros(n,nn);

I = @(x,y) -min(((x-0.5).^2+(y-0.5).^2)-0.1,0)/100; % Initial distribution in the centre of the triangle

u0=zeros(n,1);
ss =N;

%Different possible initial states
u0 = [0;0;0; I(MeshModel.Nodes(1,:),MeshModel.Nodes(2,:))'/sum(Mass{N}*I(MeshModel.Nodes(1,:),MeshModel.Nodes(2,:))');zeros((N-ss)*size(mesh{N},1),1); 0; 0; 0];
% u0(n)=10;

uht(:,1) = u0;

PP = zeros(nn,N);
PPfinal = zeros(nn/100,N);

bb=linspace(0,N,N+1);

%% Dynamic equation

[L,U,pp] = lu(MM(:,:)+h*SS(:,:)+h*VV(:,:),'vector');


for j=2:nn
    gggg=MM(:,:)*uht(:,j-1);
    uht(:,j) =  U\(L\(gggg(pp,:)));
    %Uncomment these lines to see the evolution of the process in the first
    %triangle
    %      trimesh(mesh.Elements,mesh.Nodes(1,:),mesh.Nodes(2,:),uht(1:nummesh,j))
    %     xlim([-0.1 1.1])
    %     ylim([-0.1 1.1])
    %     zlim([0 3])
    %     title(num2str(j))
    %     drawnow
end

%Density in each of the nodes
c=0;
for rrr=1:N
    s=c+1;
    c=c+size(mesh{rrr}(:,:),1);
    for i=1:nn
        PP(i,rrr) = (sum(Mass{rrr}*uht(s:c,i)))/sum(MM*u0)*100;
    end
end

%% Plots
figure;
plot(PP)