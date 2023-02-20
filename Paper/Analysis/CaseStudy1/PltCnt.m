if 1;
subplot(2,2,1); plot(x(:,1), x(:,2),'k.'); title 'HTIS';
subplot(2,2,2); plot(angles, C_theta,'ko-'); title 'Ctheta';

tmp=Cnt.XY{IMth}(:,:,iBin,IGd(iRtr),iAsc);
subplot(2,2,3); plot(tmp(:,1), tmp(:,2),'ko-'); title 'ContourEstimate';

save tDat tmp x angles C_theta;
end;

if 0;
clf; clc; clear;

load tDat;

subplot(2,2,1); plot(tmp(:,1), tmp(:,2),'ko-'); title 'ContourEstimate';
end;

if 0;
clf; hold on; 
clear;

%% Set up interesting sample
n1=700*2;
n2=200*2;
n3=100*2;
n=n1+n2+n3;
Sigma1=2*[1 0.95;0.95 1];
X1=randn(n1,2)*sqrtm(Sigma1);
Sigma2=[0.5 0;0 2];
X2=randn(n2,2)*sqrtm(Sigma2);
Sigma3=0.5*[1 0;0 1];
X3=randn(n3,2)*sqrtm(Sigma3);
X=[X1+ones(n1,1)*[0 1];X2+ones(n2,1)*[-2 2];X3+ones(n3,1)*[-1 1]];

Tau=0.999999;
nGrd=40;
CntPnt=[-0.5;4];
KrnWdt=1;
[XVal,YVal] = LocalTangent(X,Tau,nGrd,CntPnt,KrnWdt);
clf; hold on;
plot(X(:,1), X(:,2),'k.');
plot(XVal(:,1), YVal(:,1),'ro-'); 
title 'ContourEstimate';
end;