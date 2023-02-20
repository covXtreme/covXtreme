clc; clear;

addpath('R:\David.Randell\Extremes\CE\Code')

%% HT Tawn paramters
a=0.7;
b=0.5;
mu=0;
sig=1;
NEP=0.9;  %HT NEP

n=1000;  %n observations
R=randn(n,1);  %residuals

%% Simulate from HT model
tic
nRls=10000;
U= rand(nRls,1)*(1-NEP)+NEP; %sample uniform value
SmlX = -log(-log(U));  %sample gumbel value

Z=randsample(R,nRls,true);  %sample residuals
SmlY=a.* SmlX + SmlX.^b.* ( mu+ sig.*Z);
toc
     
%% Numerical Integration
tic
nG=100;
XRng = [-log(-log(NEP)),-log(-log(1-1e-5))];
XGrd=linspace(XRng(1),XRng(2),nG);
FX=((exp(-exp(-XGrd)))-NEP)./(1-NEP);
fX=exp(-(XGrd +exp(-XGrd)))./(1-NEP);
fXNrm=fX./sum(fX);

YRng = [-log(-log(1e-5)),-log(-log(1-1e-5))];
YGrd=linspace(YRng(1),YRng(2),nG)';

%find P(R< (y-aX)/x^b)
Rs=sort(R);
Ri=(YGrd-a.*XGrd)./(XGrd.^b);
Pi=CDFInterp1(Rs,Ri);

%compute integral
FY=sum(Pi.*fXNrm,2);  
toc

%%
clf;
subplot(1,3,1)
plot(SmlX,SmlY,'k.')
xlabel('X_G')
ylabel('Y_G')
axis equal 
title('MC Gumbel Margin')

subplot(1,3,2)
hold on
plot(sort(SmlX),linspace(0,1,nRls),'k-','linewidth',2)
hold on
plot(XGrd,FX,'r--','linewidth',2);
xlabel('X_G')
ylabel('cumulative probability')
title('F(X_G)')

subplot(1,3,3)
plot(sort(SmlY),linspace(0,1,nRls),'k-','linewidth',2)
hold on
plot(YGrd,FY,'r--','linewidth',2)
title('F(Y_G|X_G)')
xlabel('Y_G')
ylabel('cumulative probability')
legend('MC','NI');

savePics('HTNI')