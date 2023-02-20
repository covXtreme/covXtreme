clc; clear; clf;

addpath('R:\David.Randell\Extremes\CE\Code')

AnalysisDir='R:\David.Randell\Extremes\CE\UnitTests\Output';

load(fullfile(AnalysisDir,'MM1'),'MM')  %Load margin 1
Mrg=MM;
load(fullfile(AnalysisDir,'MM2'),'MM')  %Load margin 1
Mrg=cat(1,Mrg,MM); clear MM;
load(fullfile(AnalysisDir,'HT'),'HT')  %Load margin 1

nBin=Mrg(1).Bn.nBin;  %number of bins

rho=Mrg(1).RatExc;
nB=Mrg(1).nB;

T=100; %T year return period



%% Simulate from HT model
tic
if 0
    %% Marginal GP X
    sigX=Mrg(1).Scl; %gp scale
    xiX=repmat(Mrg(1).Shp',nBin,1);
    uX=Mrg(1).Thr; %threshold
    
    %% Marginal GP Y
    sigY=Mrg(2).Scl; %gp scale
    xiY=repmat(Mrg(2).Shp',nBin,1);
    uY=Mrg(2).Thr; %threshold
    
    %% HT Tawn paramters
    a=permute(HT.Prm(1:HT.nAlp,1,:),[1,3,2]);
    b=repmat(permute(HT.Prm(HT.nAlp+1,1,:),[1,3,2]),nBin,1);
    m=repmat(permute(HT.Prm(HT.nAlp+2,1,:),[1,3,2]),nBin,1);
    s=repmat(permute(HT.Prm(HT.nAlp+3,1,:),[1,3,2]),nBin,1);
    HTNEP=HT.NEP;  %HT HTNEP
    R=HT(1).Rsd;  %residuals
    n=length(R);  %n observations
    % simulate marignal return value
    nRls=10000;
    I=randi(nB,nRls,1);  %bootstrap samples to use
    %I=randi(1,nRls,1);  %bootstrap samples to use
    
    nStr=num2cell(poissrnd(rho(:,I)*T)); %number of storms T year period
    
    if 0
        Sml.XO=cellfun(@(x,s,u,n)max(gprnd(x,s,u,n,1)),num2cell(xiX(:,I)),...
            num2cell(sigX(:,I)),num2cell(uX(:,I)),nStr);  %max storm in T-year period
    else
        U=rand(nBin,nRls); %think about non occurence!!
        rho=Mrg(1).RatExc(:,I)+Mrg(1).RatBlw(:,I);
        Sml.XU=1+log(U)./(rho*T);    %adjust for return period  (inv of exp(-L*T*(1-C));
        Sml.XO=INV(Mrg(1),Sml.XU,I);
    end
    
    Sml.XG=INV_Standard(Mrg(1),Sml.XU);
    
    %draw random resdiuals
    Z=cell2mat(cellfun(@(x)x(randi(numel(x),nBin,1))',R(I),'uniformoutput',false))';
    
    % simulate under HT model
    Sml.YG=a(:,I).* Sml.XG + Sml.XG.^b(:,I).* ( m(:,I)+ s(:,I).*Z);
    
    %transform from standard to uniform margins using CDF
    Sml.YU=CDF_Standard(Mrg(1),Sml.YG);
    
    tP=(Sml.YU-Mrg(2).NEP(I)')./(1-Mrg(2).NEP(I)');
    tP(tP<0)=NaN;
    Sml.YO=gpinv(tP,xiY(:,I),sigY(:,I),uY(:,I));  %inv
    
    [Sml.XOmni,J]=max(Sml.XO);
    
    WMC=NaN(nBin,nB);
    for iB=1:nB
        WMC(:,iB)=accumarray(J(I==iB)',J(I==iB)',[nBin,1],@numel);
    end
    WMC=WMC./sum(WMC);
    
    
    [Sml.YOmni]=Sml.YO(sub2ind([nBin,nRls],J,(1:nRls)));
else %test HT code
    Mrg(1).RtrPrd=T;
    HT.nRtr=1;
    HT=ConditionalReturnValue(HT,Mrg);
    
    Sml.XG=HT.RV.X_Stn;
    Sml.XO=HT.RV.X;
    Sml.XOmni=HT.RV.XOmni;
    Sml.YG=squeeze(HT.RV.Y_Stn);
    Sml.YO=squeeze(HT.RV.Y);
    Sml.YOmni=squeeze(HT.RV.YOmni);
end

toc                               
     
%% Numerical Integration to compute conditional Return Value
tic
%% CDF X  X--> U
CDF_X=permute(Mrg(1).RVX,[3,2,1]);

XGrdU=CDF(Mrg(1),CDF_X);

%% PIT  U--> Gumbel
XGrdStn=INV_Standard(Mrg(1),XGrdU);
XGrdStnReg=permute(linspace(0,10,numel(Mrg(1).RVX)),[1,4,3,2]);
XGrdStn=permute(XGrdStn,[1,2,4,3]); %transpose to find all combinations of X/Y

nRtr=Mrg(1).nRtr;
RVCnd=NaN(HT.nBin+1,Mrg(1).nRVX); %preallocate

%Marginal return value exp(-rho.*T (1-C))
FX=exp(-bsxfun(@times,(Mrg(1).RatExc+Mrg(1).RatBlw).*T,(1-XGrdU)));
fX=cat(3,diff(FX,1,3),zeros(HT.nBin,HT.nB));

sF=sum(fX,3);
PrbZero=sum(sF==0,2)./HT.nB;

fXNrm=fX./sF; %if have zero rate get NaNs in sum
fXNrm(isnan(fXNrm))=0;
fXNrm=permute(fXNrm,[1,2,4,3]);

FXOmni=prod(FX);
    


iDmn=2;
%% CDF Y
CDF_Y=permute(Mrg(iDmn).RVX,[3,2,1]);

dX=diff(Mrg(iDmn).RVX(1:2));
YGrdU=CDF(Mrg(iDmn),CDF_Y);
YGrdStn=INV_Standard(Mrg(iDmn),YGrdU);
%find P(R< (y-aX)/x^b) cdf of normalised residuals

Pi=NaN(HT.nBin,HT.nB,Mrg(iDmn).nRVX,Mrg(1).nRVX);

for iB=1:HT.nB
    Rs=sort(HT.Rsd{iB}(:,iDmn-1));
    Xb=XGrdStn(:,iB,:,:).^HT.Prm(end-2,iDmn-1,iB);
    Ri=(YGrdStn(:,iB,:)-HT.Prm(1:HT.nAlp,iDmn-1,iB).*XGrdStn(:,iB,:,:)-Xb.*HT.Prm(end-1,iDmn-1,iB))./(HT.Prm(end,iDmn-1,iB).*Xb);
    Pi(:,iB,:,:)=CDFInterp1(Rs,Ri);
end

%compute integral
tFY=sum(Pi.*fXNrm,4)+PrbZero;

W=fX./FX.*FXOmni;
%W=W./nansum(W,3);

FYOmni=sum(nansum(Pi.*permute(W,[1,2,4,3]),1),4);
FYOmni=FYOmni./max(FYOmni,[],3);
FY=[tFY;FYOmni];

RVCnd=squeeze(mean(FY,2));
RVCnd(RVCnd>=1)=1; %check well defined CDF
toc



%% Plot
nRls=HT.RV.nRls;
figure(1);
clf;
subplot(1,3,1)

hold on
plot(Sml.XO,Sml.YO,'.k')
plot(Mrg(1).Y,Mrg(2).Y,'.','color',[1,1,1]*0.7);
hold on
plot(Sml.XOmni,Sml.YOmni,'og')
xlabel('X')
ylabel('Y')
%axis equal 
title('MC Original Margin')

subplot(2,3,[2,3])
hold on
plot(sort(Sml.XO,2,'MissingPlacement','first'),linspace(0,1,nRls),'-','linewidth',2)
hold on
plot(sort(max(Sml.XO),2,'MissingPlacement','first'),linspace(0,1,nRls),'k-','linewidth',2)
plot(Mrg(1).RVX,squeeze(mean(FX,2)),'r--','linewidth',2);
hold on
plot(Mrg(1).RVX,squeeze(mean(FXOmni,2)),'g--','linewidth',2);
xlabel('X')
ylabel('cumulative probability')
title('F(X)')
%xlim([5,20])

subplot(2,3,[5,6])
plot(sort(Sml.YO,2,'MissingPlacement','first'),linspace(0,1,nRls),'-','linewidth',2)
% hold on
% plot(sort(max(Sml.YO),2),linspace(0,1,nRls),'k-','linewidth',2)
hold on
plot(sort(Sml.YOmni,'MissingPlacement','first'),linspace(0,1,nRls),'k-','linewidth',2)
plot(Mrg(2).RVX,RVCnd(1:end-1,:),'r--','linewidth',2)
plot(Mrg(2).RVX,RVCnd(end,:),'g--','linewidth',2)
title('F(Y|X)')
xlabel('Y')
ylabel('cumulative probability')
%xlim([10,20])

%savePics('HTNI')

%%
%figure(2)
% clf;
% subplot(2,1,1)
% histogram(XGrdU)
% subplot(2,1,2)
% histogram(YGrdU)
% 
% XGrdU