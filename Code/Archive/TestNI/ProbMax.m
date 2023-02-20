

clf;close;

mu=[1,3,0];
sig=[1,1,2];

nR=100000;
nB=3;

Z=randn(nR,nB).*sig+mu;

[~,I]=max(Z,[],2);

p=accumarray(I,I,[nB,1],@numel,0)./nR;


figure(1)
clf;
subplot(2,1,1)

hold on;
for i=1:3
    histogram(Z(:,i),linspace(-15,15,100),'displaystyle','stairs','normalization','pdf')
end

%% 
nG=100;
X=linspace(-15,15,nG)';

C=NaN(nG,nB);
P=NaN(nG,nB);

for iB=1:nB
    C(:,iB)=normcdf(X,mu(iB),sig(iB));
    P(:,iB)=normpdf(X,mu(iB),sig(iB));
end

mean(P./sum(P,2))


p2=NaN(nB,1);
for iB=1:3
    J=setdiff(1:3,iB);
    Y=prod(C(:,J),2); %max of 2,3
    pY=[diff(Y);0]; %py
    X=C(:,iB); %1
    
    p2(iB)=1-sum(X.*pY);
    
    fprintf('Bin %g: Analytical: %g, Empirical: %g\n',iB,p2(iB),p(iB))
end








