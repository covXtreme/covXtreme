clc; clear;

addpath('R:\David.Randell\Extremes\CE\Code')

nBin=4;  %number of bins
nB=100;  %number of bootstraps

rho=[2,1,2,1]'; %annual rate
rho=rho+randn(nBin,nB)*0.1;

T=1000; %T year return period
%% Marginal GP X
sigX=[2,1.5,1.8,2.3]'; %gp scale
sigX=sigX+randn(nBin,nB)*0.1;
xiX=-0.1; %gp shape
xiX=xiX+randn(nBin,nB)*0.05;
uX=ones(nBin,1)+randn(1,nB)*0.1; %threshold

%% Marginal GP Y
sigY=[2.9,1.7,2.8,3.3]'; %gp scale
sigY=sigY+randn(nBin,nB)*0.1;
xiY=-0.2; %gp shape
xiY=xiY+randn(nBin,nB)*0.05;
uY=1.5*ones(nBin,1)+randn(1,nB)*0.1; %threshold

%% HT Tawn paramters
a=0.7+randn(1,nB)*0.01;
b=0.5+randn(1,nB)*0.01;
m=0+randn(1,nB)*0.01;
s=1+randn(1,nB)*0.01;
HTNEP=0.5;  %HT HTNEP

n=1000;  %n observations
R=randn(n,1);  %residuals

%% Simulate from HT model
tic
% simulate marignal return value
nRls=10000;
I=randi(nB,nRls,1);  %bootstrap samples to use

nStr=num2cell(poissrnd(rho(:,I)*T)); %number of storms T year period
            
Sml.XO=cellfun(@(x,s,u,n)max(gprnd(x,s,u,n,1)),num2cell(xiX(:,I)),...
    num2cell(sigX(:,I)),num2cell(uX(:,I)),nStr);  %max storm in T-year period

Sml.XU= gpcdf(Sml.XO,xiX(:,I),sigX(:,I),uX(:,I)); %sample uniform value
Sml.XG = -log(-log(Sml.XU));  %sample gumbel value

% simulate under HT model
J=randi(n,nBin,nRls);  %sample residuals
Z=R(J);

Sml.YG=a(:,I).* Sml.XG + Sml.XG.^b(:,I).* ( m(:,I)+ s(:,I).*Z);
Sml.YU=exp(-exp(-Sml.YG)); %cdf
Sml.YO=gpinv(Sml.YU,xiY(:,I),sigY(:,I),uY(:,I));  %inv

toc                               
     
%% Numerical Integration
tic
nG=100;
%XRng on Original Scale

XRng = [0,30];
XGrdOrg=linspace(XRng(1),XRng(2),nG)';
XGrdU=ReturnValue.GeneralisedPareto.gpcdf(XGrdOrg,shiftdim(xiX,-1),shiftdim(sigX,-1),shiftdim(uX,-1));
XGrdGmb = -log(-log(XGrdU));
XGrdGmb=permute(XGrdGmb,[4,2,3,1]); %transpose to find all combinations of X/Y

%Marignal return value
FX=exp(-shiftdim(rho,-1).*T.*(1-XGrdU));
fX=[diff(FX);zeros(1,nBin,nB)];
fXNrm=fX./sum(fX,1);
fXNrm=permute(fXNrm,[4,2,3,1]);

FXOmni=mean(prod(FX,2),3);
FX=mean(FX,3);

YRng = [0,30];
YGrdOrg=linspace(YRng(1),YRng(2),nG)';
YGrdU=ReturnValue.GeneralisedPareto.gpcdf(YGrdOrg,shiftdim(xiY,-1),shiftdim(sigY,-1),shiftdim(uY,-1));
YGrdGmb= -log(-log(YGrdU));

%find P(R< (y-aX)/x^b) cdf of normalised residuals
Rs=sort(R);
Xb=XGrdGmb.^shiftdim(b,-1);
Ri=(YGrdGmb-shiftdim(a,-1).*XGrdGmb-Xb.*shiftdim(m,-1))./(shiftdim(s,-1).*Xb);
Pi=CDFInterp1(Rs,Ri);

%compute integral
FY=sum(Pi.*fXNrm,4);  
FYOmni=mean(prod(FY,2),3);  
FY=mean(FY,3);
toc

%% Plot
clf;
%subplot(2,3,1)
% plot(Sml.XG,Sml.YG,'k.')
% xlabel('X_G')
% ylabel('Y_G')
% axis equal 
% title('MC Gumbel Margin')

% subplot(2,3,2)
% hold on
% plot(sort(Sml.XG,2),linspace(0,1,nRls),'k-','linewidth',2)
% hold on
% plot(XGrdGmb,FX,'r--','linewidth',2);
% xlabel('X_G')
% ylabel('cumulative probability')
% title('F(X_G)')

% subplot(2,3,3)
% plot(sort(Sml.YG,2),linspace(0,1,nRls),'k-','linewidth',2)
% hold on
% plot(YGrdGmb,FY,'r--','linewidth',2)
% title('F(Y_G|X_G)')
% xlabel('Y_G')
% ylabel('cumulative probability')
% legend('MC','NI');

subplot(1,3,1)
plot(Sml.XO,Sml.YO,'k.')
xlabel('X')
ylabel('Y')
%axis equal 
title('MC Original Margin')

subplot(1,3,2)
hold on
plot(sort(Sml.XO,2),linspace(0,1,nRls),'k-','linewidth',2)
hold on
plot(sort(max(Sml.XO),2),linspace(0,1,nRls),'k-','linewidth',2)
plot(XGrdOrg,FX,'r--','linewidth',2);
hold on
plot(XGrdOrg,FXOmni,'g--','linewidth',2);
xlabel('X')
ylabel('cumulative probability')
title('F(X)')

subplot(1,3,3)
plot(sort(Sml.YO,2),linspace(0,1,nRls),'k-','linewidth',2)
hold on
plot(sort(max(Sml.YO),2),linspace(0,1,nRls),'k-','linewidth',2)
plot(YGrdOrg,FY,'r--','linewidth',2)
plot(YGrdOrg,FYOmni,'g--','linewidth',2)
title('F(Y|X)')
xlabel('Y')
ylabel('cumulative probability')


savePics('HTNI')