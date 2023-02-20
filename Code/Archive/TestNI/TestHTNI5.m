clc; clear; 

addpath('R:\David.Randell\Extremes\CE\Code')

%AnalysisDir='R:\Public Folder\Extremes\Analysis\NorthSea\SurgeAnalysis\Brent\Rng_NonStn\Output';
AnalysisDir='R:\Public Folder\Extremes\Analysis\NorthSea\SurgeAnalysis\Anga\Rng_NonStn\Output';
load(fullfile(AnalysisDir,'MM1'),'MM')  %Load margin 1
Mrg=MM;
load(fullfile(AnalysisDir,'MM2'),'MM')  %Load margin 1
Mrg=cat(1,Mrg,MM); clear MM;
load(fullfile(AnalysisDir,'HT'),'HT')  %Load margin 1

nBin=size(Mrg(1).Scl,1);  %number of bins

rho=Mrg(1).RatExc;
nB=Mrg(1).nB;

T=10; %T year return period
%% Marginal GP X
sigX=Mrg(1).Scl; %gp scale
xiX=repmat(Mrg(1).Shp',nBin,1);
uX=Mrg(1).Thr; %threshold

%% Marginal GP Y
sigY=Mrg(2).Scl; %gp scale
xiY=repmat(Mrg(2).Shp',nBin,1);
uY=Mrg(2).Thr; %threshold

%% HT Tawn paramters
a=HT.Prm(1:nBin,:);
b=repmat(HT.Prm(nBin+1,:),nBin,1);
m=repmat(HT.Prm(nBin+1,:),nBin,1);
s=repmat(HT.Prm(nBin+1,:),nBin,1);
HTNEP=HT.NEP;  %HT HTNEP


R=randn(1000,1);  %residuals
n=length(R);  %n observations
%% Simulate from HT model
tic
% simulate marignal return value
nRls=10000;
I=randi(nB,nRls,1);  %bootstrap samples to use

nStr=num2cell(poissrnd(rho(:,I)*T)); %number of storms T year period
    
StrOn=cell2mat(nStr)>0;
Sml.XO=NaN(nBin,nRls);
tXi=num2cell(xiX(:,I));
tSig=num2cell(sigX(:,I));
tU=num2cell(uX(:,I));
Sml.XO(StrOn)=cellfun(@(x,s,u,n)max(gprnd(x,s,u,n,1)),tXi(StrOn),tSig(StrOn),tU(StrOn),nStr(StrOn));  %max storm in T-year period

Sml.XU= gpcdf(Sml.XO,xiX(:,I),sigX(:,I),uX(:,I)); %sample uniform value
Sml.XG = -log(-log(Sml.XU));  %sample gumbel value

% simulate under HT model
J=randi(n,nBin,nRls);  %sample residuals
Z=R(J);



Sml.YG=a(:,I).* Sml.XG + Sml.XG.^b(:,I).* ( m(:,I)+ s(:,I).*Z);
Sml.YU=exp(-exp(-Sml.YG)); %cdf
Sml.YO=gpinv(Sml.YU,xiY(:,I),sigY(:,I),uY(:,I));  %inv

[Sml.XOmni,J]=max(Sml.XO);

WMC=accumarray(J',J',[nBin,1],@numel);

[Sml.YOmni]=Sml.YO(sub2ind([nBin,nRls],J,(1:nRls)));

toc                               
     
%% Numerical Integration
tic
nG=400;
%XRng on Original Scale

XRng = [0,30];
XGrdOrg=linspace(XRng(1),XRng(2),nG)';
XGrdU=ReturnValue.GeneralisedPareto.gpcdf(XGrdOrg,shiftdim(xiX,-1),shiftdim(sigX,-1),shiftdim(uX,-1));
XGrdGmb = -log(-log(XGrdU));
XGrdGmb=permute(XGrdGmb,[4,2,3,1]); %transpose to find all combinations of X/Y

fpdf=ReturnValue.GeneralisedPareto.gppdf(XGrdOrg,shiftdim(xiX,-1),shiftdim(sigX,-1),shiftdim(uX,-1));

%Marignal return value
rT=shiftdim(rho,-1).*T;
FX=exp(-shiftdim(rho,-1).*T.*(1-XGrdU));
%fX=[diff(FX);zeros(1,nBin,nB)];
fX=FX.*rT.*fpdf;

fXNrm=fX./sum(fX,1);

FXOmni=prod(FX,2);


%% Prob of max bin

if 0
    W=nansum(fXNrm./FX.*FXOmni);
   % W=W./sum(W,2); %make sure sum to 1 
else
    
    W=NaN(1,nBin,nB);
    
    for iB=1:nBin
        J=setdiff(1:nBin,iB);
        Y=prod(FX(:,J,:),2); %max of 2,3
        
        pY=[diff(Y);zeros(1,1,nB)]; %py
        X=FX(:,iB,:); %1
        
        W(1,iB,:)=1-sum(X.*pY);
    end
    W=W./sum(W,2); %make sure sum to 1
end

FXOmni=mean(FXOmni,3);
FX=mean(FX,3);



%% Conditional Y
YRng = [min(Mrg(2).Y),max(Mrg(2).Y)*3];
YGrdOrg=linspace(YRng(1),YRng(2),nG)';
YGrdU=ReturnValue.GeneralisedPareto.gpcdf(YGrdOrg,shiftdim(xiY,-1),shiftdim(sigY,-1),shiftdim(uY,-1));
YGrdGmb= -log(-log(YGrdU));

%find P(R< (y-aX)/x^b) cdf of normalised residuals
Rs=sort(R);
Xb=XGrdGmb.^shiftdim(b,-1);
Ri=(YGrdGmb-shiftdim(a,-1).*XGrdGmb-Xb.*shiftdim(m,-1))./(shiftdim(s,-1).*Xb);
Pi=CDFInterp1(Rs,Ri);

%compute integral
fXNrm=permute(fXNrm,[4,2,3,1]);
FY=sum(Pi.*fXNrm,4);  
FYOmni=mean(sum(W.*FY,2),3);  
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
plot(Mrg(1).Y,Mrg(2).Y,'.','color',[1,1,1]*0.7);
hold on
plot(Sml.XO,Sml.YO,'k.')
hold on
plot(Sml.XOmni,Sml.YOmni,'go')
xlabel('X')
ylabel('Y')
%axis equal 
title('MC Original Margin')

subplot(2,3,[2,3])
hold on
plot(sort(Sml.XO,2,'MissingPlacement','first'),linspace(0,1,nRls),'k-','linewidth',2)
hold on
plot(sort(max(Sml.XO),2,'MissingPlacement','first'),linspace(0,1,nRls),'b-','linewidth',2)
plot(XGrdOrg,FX,'r--','linewidth',2);
hold on
plot(XGrdOrg,FXOmni,'g--','linewidth',2);
xlabel('X')
ylabel('cumulative probability')
title('F(X)')

subplot(2,3,[5,6])
plot(sort(Sml.YO,2,'MissingPlacement','first'),linspace(0,1,nRls),'k-','linewidth',2)
% hold on
% plot(sort(max(Sml.YO),2),linspace(0,1,nRls),'k-','linewidth',2)
hold on
plot(sort(Sml.YOmni,'MissingPlacement','first'),linspace(0,1,nRls),'b-','linewidth',2)
plot(YGrdOrg,FY,'r--','linewidth',2)
plot(YGrdOrg,FYOmni,'g--','linewidth',2)
title('F(Y|X)')
xlabel('Y')
ylabel('cumulative probability')


%savePics('HTNI')