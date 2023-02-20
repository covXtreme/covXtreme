clf;

nB=10;
R=sort(randn(1000,nB)); %residuals


X=linspace(1,5,1000)';
Y=linspace(1,5,1000);

a=0.7+randn(nB,1)*0.05;
b=0.5+randn(nB,1)*0.05;
m=0+randn(nB,1)*0.05;
s=1+randn(nB,1)*0.05;
Xb=X.^shiftdim(b,-2);
Ri=(Y-shiftdim(a,-2).*X-Xb.*shiftdim(m,-2))./(shiftdim(s,-2).*Xb);

tic
for iR=1:100
Pi=NaN(size(Ri));
for iB=1:nB
    Pi(:,:,iB)=CDFInterp1(R(:,iB),Ri(:,:,iB));
end
end
toc

tic
for iR=1:100
tP=linspace(0,1,numel(R))';
I2=knnsearch(R,reshape(Ri,[],nB));
Pi2=reshape(tP(I2),size(Pi));
end
toc

tic
for iR=1:100
tP=linspace(0,1,numel(R))';
Pi3=interp1(tP,Rs,Ri);
end
toc

tic
for iR=1:100
tP=linspace(0,1,numel(R))';
Pi4 = qinterp1(tP,Rs,Ri,1);
end
toc

clf;
plot(Pi2-Pi)