load CElegansMetaplex.mat
MaxTries=10;

for Try=1
realsize=0;
for i=1:length(A)
    realsize=realsize+length(mesh{i});
end

u0=zeros(realsize,1);
% u0(randi(280,50,1))=1;
u0([1:280])=1;
% batch(@DynamicMetaplexBigNet,2,{A,Mass,mesh,sink,source,u0,Try})
DynamicMetaplexBigNet(A,Mass,mesh,sink,source,u0,Try);

% mkdir("/home/mmiranda/Desktop/ManuMetaplex/NuredddunaBigMetaplexes/MetaplexTry"+convertCharsToStrings(sprintf('%d',Try)));
% fileplace = "/home/mmiranda/Desktop/ManuMetaplex/NuredddunaBigMetaplexes/MetaplexTry"+convertCharsToStrings(sprintf('%d',Try));
% 
% c=1;
% for i = 1:length(A)
% temp= uht(c:c+length(mesh{i})-1,:);
% filename = fileplace+"/u"+convertCharsToStrings(sprintf('%d',i));
% save(filename,"temp");
% c=c+length(mesh{i});
% end

% clear uht
% clear all
end