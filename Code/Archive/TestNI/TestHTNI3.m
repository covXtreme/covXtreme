clc; clear;

addpath('R:\David.Randell\Extremes\CE\Code')

rho=2; %annual rate
T=100; %T year return period
%% Marginal GP X
sigX=2; %gp scale
xiX=-0.1; %gp shape
uX=0; %threshold

%% Marginal GP Y
sigY=3; %gp scale
xiY=-0.2; %gp shape
uY=0; %threshold

%% HT Tawn paramters
a=0.7;
b=0.5;
m=0;
s=1;
HTNEP=0.5;  %HT HTNEP

n=1000;  %n observations
R=randn(n,1);  %residuals

%% Simulate from HT model
tic
% simulate marignal return value
nRls=10000;
nStr=num2cell(poissrnd(rho*T,nRls,1)); %number of storms T year period
Sml.XO=cellfun(@(n)max(gprnd(xiX,sigX,uX,n,1)),nStr);  %max storm in T-year period
Sml.XU= gpcdf(Sml.XO,xiX,sigX,uX); %sample uniform value
Sml.XG = -log(-log(Sml.XU));  %sample gumbel value

% simulate under HT model
Z=randsample(R,nRls,true);  %sample residuals
Sml.YG=a.* Sml.XG + Sml.XG.^b.* ( m+ s.*Z);
Sml.YU=exp(-exp(-Sml.YG)); %cdf
Sml.YO=gpinv(Sml.YU,xiY,sigY,uY);  %inv

toc                               
     
%% Numerical Integration
tic
nG=100;
%XRng on Original Scale

XRng = [uX,gpinv(1-1e-5,xiX,sigX,uX)];
XGrdOrg=linspace(XRng(1),XRng(2),nG);
XGrdU = gpcdf(XGrdOrg,xiX,sigX,uX);
XGrdGmb = -log(-log(XGrdU));

%Marignal return value
FX=exp(-rho.*T*(1-XGrdU));
fX=[diff(FX),0];
fXNrm=fX./sum(fX);

YRng = [0,15];
YGrdOrg=linspace(YRng(1),YRng(2),nG)';
YGrdU=gpcdf(YGrdOrg,xiY,sigY,uY);
YGrdGmb= -log(-log(YGrdU));

%find P(R< (y-aX)/x^b) cdf of normalised residuals
Rs=sort(R);
Xb=XGrdGmb.^b;
Ri=(YGrdGmb-a.*XGrdGmb-Xb.*m)./(s.*Xb);
Pi=CDFInterp1(Rs,Ri);

%compute integral
FY=sum(Pi.*fXNrm,2);  
toc

%% Plot
clf;
subplot(2,3,1)
plot(Sml.XG,Sml.YG,'k.')
xlabel('X_G')
ylabel('Y_G')
axis equal 
title('MC Gumbel Margin')

subplot(2,3,2)
hold on
plot(sort(Sml.XG),linspace(0,1,nRls),'k-','linewidth',2)
hold on
plot(XGrdGmb,FX,'r--','linewidth',2);
xlabel('X_G')
ylabel('cumulative probability')
title('F(X_G)')

subplot(2,3,3)
plot(sort(Sml.YG),linspace(0,1,nRls),'k-','linewidth',2)
hold on
plot(YGrdGmb,FY,'r--','linewidth',2)
title('F(Y_G|X_G)')
xlabel('Y_G')
ylabel('cumulative probability')
legend('MC','NI');

subplot(2,3,4)
plot(Sml.XO,Sml.YO,'k.')
xlabel('X')
ylabel('Y')
axis equal 
title('MC Original Margin')

subplot(2,3,5)
hold on
plot(sort(Sml.XO),linspace(0,1,nRls),'k-','linewidth',2)
hold on
plot(XGrdOrg,FX,'r--','linewidth',2);
xlabel('X')
ylabel('cumulative probability')
title('F(X)')

subplot(2,3,6)
plot(sort(Sml.YO),linspace(0,1,nRls),'k-','linewidth',2)
hold on
plot(YGrdOrg,FY,'r--','linewidth',2)
title('F(Y|X)')
xlabel('Y')
ylabel('cumulative probability')


savePics('HTNI')