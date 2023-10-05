function [PP,uht] = DynamicMetaplexBigNet(A,Mass,mesh, sink, source,u0)
    %DynamicMetaplexBigNet.m - Computes a diffusion process on the metaplex with
    %subjacent network A without storing the matrices of the dynamics.
    % It is useful when MATLAB can not handle matrices of the needed dimension.
    % It can be used with the outputs of function
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


nn = 10^4;    % Number of time steps used to find uht
h = 10^-4;     % Time step

N=length(A);
realsize=0;
for i=1:N
    realsize=realsize+length(mesh{i});
end
alpha=ones(N,N);

n = realsize;
uht = zeros(n,nn);

%Different possible initial states
% u0=zeros(n,1);
% u0([1:100])=1;

uht(:,1) = u0;

PP = zeros(nn,N);


for j=2:nn
    c=0;
    for i=1:N
        s=c+1;
        c=c+size(mesh{i}(:,:),1);
        uht(s:c,j)=(Mass{i}(:,:)- h*mesh{i}(:,:))*uht(s:c,j-1);
        t=0;
        for l=1:N
            if ~isempty(source{i,l}) && ~isempty(sink{l,i})
                uht(s-1+source{i,l},j)=uht(s-1+source{i,l},j)+h*alpha(i,l)*A(i,l)*uht(t+sink{l,i},j-1);
            end
            if ~isempty(sink{i,l})
                uht(s-1+sink{i,l},j)=uht(s-1+sink{i,l},j)-h*alpha(i,l)*A(i,l)*uht(s-1+sink{i,l},j-1);
            end
            t=t+size(mesh{l}(:,:),1);
        end
        uht(s:c,j) = Mass{i}(:,:)^-1*uht(s:c,j);
    end
end

c=0;
InitialDens=0;
for i=1:N
    s=c+1;
    c=c+size(mesh{i}(:,:),1);
    InitialDens=InitialDens+sum(Mass{i}(:,:)*u0(s:c));
end

%Density in each of the nodes
c=0;
for rrr=1:N
    s=c+1;
    c=c+size(mesh{rrr}(:,:),1);
    for i=1:nn
        PP(i,rrr) = (sum(Mass{rrr}*uht(s:c,i)))/InitialDens*100;
    end
end
figure,plot(PP)

end