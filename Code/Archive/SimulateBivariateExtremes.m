clear; clc;

nDat=10000;
DrcEdg=[10,60,140,340]'; %NB> when nBin=1 Loc of DrcEdg doesn't matter...stationary
nBin=length(DrcEdg);
xi=-0.1;  %shape
sig=[2,3,1,1.8]';  %scale
PrmTrue=[xi;sig];
thr=zeros(nBin,1);  %threshold

if exist('DrcX','var') && exist('Y','var') %if already have data
    MM = MarginalModel(DrcX,Y,DrcEdg,PrmTrue); 
else
    MM=MarginalModel();   %empty constructor
    MM=SimulateMarginal(MM,nDat,DrcEdg,PrmTrue,thr); %Simulate data if don't have any    
end
Tau=0.8;    %non exceedence probability
nBt=1;  %if want to bootstrap, make nBt > 1

%starting param ------
p0=[0;ones(nBin,1)*2];  

%CV settings --------
CVAllBt=0; %=1 then CV smoothness for every bootstrap resample (slow), =0 use opt. smoothness from orig dataset CV only 
nSigLam=10;  %no. smoothness penalties tried in CV (limit variation across neighbouring bins)
SigLam=logspace(-8,4,nSigLam); %try range smoothness penalties for sigma varying by bin
nCV=10; %no. cross-validation groups

%Bootstrap -----

MM = Fit(MM,Tau,nBt,CVAllBt,SigLam,p0,nCV);    %TODO merge FitMargModel and BootMargModel


% %% marginal i  
% MM{i}.Xi =     %nB x 1
% MM{i}.Sig   %nBin_i x nB
% MM{i}.Thr   % nB x 
% MM{i}.Rho=   % nB x nBin_i
% MM{i}.BinEdg   %nBin_i x 1

%% simulate some data from assymertic logistic model

n=10000;  %number of data
alp=[NaN,NaN,0.5];  %three parameters in nD=2;
theta={0,0.5,[1,0.5]};
F=rndasymlgs(alp,theta,nD,n)';  %simulated on frechet scale

G=log(F);
%% transform to uniform scale
U=exp(-(1./F));

Y=NaN(n,nD);  %response
X=NaN(n,1);  %covariates


Edg=[MM(1).DrcEdg(end)-360;MM(1).DrcEdg];     %bin edges (add final bin for wrapping)
BinSze=Edg(2:end)-Edg(1:end-1);  %bin width
%Simulate covariate with right rate
RatCdf=cumsum(MM(1).CntExc)./sum(MM(1).CntExc); %get rate cdf
[~,A]=histc(rand(n,1),[0;RatCdf]);  %bin allocation
X=mod(rand(n,1).*BinSze(A)+Edg(A),360);  %covariate value within bin        

for i=1:nD     %loop over dimension                     
    %Simulate reponse
    Y(:,i)=gpinv(U(:,i),MM(i).XiTrue,MM(i).SigTrue(A),0);
end

clf;
subplot(2,2,1)
plot(X,Y(:,1),'k.')
subplot(2,2,2)
plot(X,Y(:,2),'k.')
subplot(2,2,3)
plot(U(:,1),U(:,2),'k.')
subplot(2,2,4)
scatter(Y(:,1),Y(:,2),20,A,'filled')