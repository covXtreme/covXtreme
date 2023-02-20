
n=1000;
nB=100;  %number of bootsraps
X=randn(n,nB);

nRls=1e3;
nBin=10;
Bn=linspace(-3,3,nBin);

tic
for iR=1:nRls
    [~,A]=histc(X(:),Bn);
    A=reshape(A,size(X));
end
toc;

tic
for iR=1:nRls
    A2=discretize(X,Bn);
end
toc;

tic
for iR=1:nRls
    nBn=numel(Bn);    
    dX=1./(Bn(2)-Bn(1));
        
    A3=ceil((X-Bn(1)).*dX);
    
    A3(A3<1)=NaN;    
    A3(A3>length(Bn))=NaN;
end
toc


[A,A3];

