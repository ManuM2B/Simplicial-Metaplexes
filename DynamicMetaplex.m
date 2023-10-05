function [uht,PP] = DynamicMetaplex(A,Mass,mesh, sink, source,u0)
    %DynamicMetaplex.m - Computes a diffusion process on the metaplex with
    %subjacent network A. It can be used with the outputs of function
    %DefaultSimplex.m
    %
    %INPUTS:
    %   -A: Graph of the simplicial complex, where the edges and triangles
    %   are also considered as nodes.
    %   -Mass: One dimensional cell of matrices representing the structure of each node in
    %   the metaplex.
    %   -mesh: One dimensional cell of matrices representing the diffusion inside each
    %   node in the metaplex.
    %   -Sink/Source: Two dimensional cells containing rectangular matrices
    %   whose position and values indicate the location and the strength of
    %   the link between two nodes.
    %   -u0: Initial conditions of the system. There are some examples of
    %   possible initial conditions commented in the code.
    %
    %OUTPUTS:
    %   -uht: Matrix whose columns give the state of the function in each
    %   of the points inside all nodes at each time step.
    %   -PP: Matrix whose columns give the density of the function for each
    %   of the nodes at each time step.

%% Matrix generation    

N=length(A);
K=sum(A,2);

%Alpha is the transition strength between different nodes. It is one by
%default, but the next lines can be uncommented to 
alpha=ones(N,N);
% for i=1:N
%     for j=1:N
%         alpha(i,j)=1/K(i);
%     end
% end

%Compute the final size of the matrices when all of the meshses are joint
realsize=0;
for i=1:N
    realsize=realsize+length(mesh{i});
end
% SS=zeros(realsize);
% MM=zeros(realsize);
VV=zeros(realsize);

% SS is the block-diagonal matrix of the discretized diffusion process.
% MM is the block-diagonal matrix of the meshes in the nodes.
% c=0;
% for i=1:N
%     s=c+1;
%     c=c+size(mesh{i}(:,:),1);
%     SS(s:c,s:c)=mesh{i}(:,:);
%     MM(s:c,s:c)=Mass{i}(:,:);
% end
SS=blkdiag(mesh{:});
MM=blkdiag(Mass{:});
% MM=MM/(sum(MM,'All'));


% VV is the matrix of transitions between nodes. These fors are used to
% fill the matrix. 
a=0; %This parameter is for taking into account the size of the mesh of the row node
for i=1:N %Iteration over nodes for the rows
    b=0; %This parameter is for taking into account the size of the mesh of the column node
    for j=1:N %Iteration over nodes for the columns
        if A(i,j)~=0
            maxlen=max(length(source{i,j}),length(sink{i,j})); %Sinks and sources should have equal length
            for l=1:maxlen %Iteration over the sink and source points
                VV(a+source{i,j}(l),b+sink{j,i}(l))=VV(a+source{i,j}(l),b+sink{j,i}(l))-A(i,j)*alpha(i,j);
                VV(b+sink{j,i}(l),b+sink{j,i}(l))=VV(b+sink{j,i}(l),b+sink{j,i}(l))+A(i,j)*alpha(i,j);
            end
        end
    b=b+size(mesh{j}(:,:),1); %Update the size of the column node
    end
    a=a+size(mesh{i}(:,:),1); %Update the size of the row node
end

%% Dynamics

nn = 10^4;    % Number of time steps used to find uht
h = 10^-4;     % Time step

n = size(SS,1);
uht = zeros(n,nn);

%Initial condition
% u0=zeros(n,1);

% u0=ones(n,1);
% u0(3)=100;
% u0(273)=10;
% u0(end)=1;
% u0(6:26)=1;
% u0(end-210:end)=1;
% u0=randi(100,n,1);

uht(:,1) = u0;
% uht(:,1) = MM^-1*u0;

PP = zeros(nn,N);  

bb=linspace(0,N,N+1);
    
%Equation implementation

[L,U,pp] = lu(MM(:,:)+h*SS(:,:)+h*VV(:,:),'vector');
for j=2:nn
    gggg=MM(:,:)*uht(:,j-1);
    uht(:,j) =  U\(L\(gggg(pp,:)));
end

%Computation of the density in each node
c=0;
for rrr=1:N
    s=c+1;
    c=c+size(mesh{rrr}(:,:),1);
for i=1:nn
PP(i,rrr) = sum(Mass{rrr}*uht(s:c,i))/sum(MM*u0)*100;
end
end
figure, plot(PP)

% % This for does the same without the LU decompostion
% for j=2:nn
%    uht(:,j) = (MM+h*(SS+VV))^-1*MM*uht(:,j-1); 
% end
% 
% % Computation of the density in each node
% c=0;
% for rrr=1:N
%     s=c+1;
%     c=c+size(mesh{rrr}(:,:),1);
% for i=1:nn
% PP(i,rrr) = sum(Mass{rrr}*uht(s:c,i))/sum(MM*u0)*100;
% end
% end
% 
% figure,plot(PP)

end