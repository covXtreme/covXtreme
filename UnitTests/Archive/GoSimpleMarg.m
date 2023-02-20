clear; clc;
addpath('../Code')
%addpath('R:\Emma.Ross\Extremes\Devel_Code\BARS') %for savePics only

%Code strucure
%- MarginalModel (class) < +Basis\PPC
%- Hefferenan Tawn (class)
%
%TODO Marginal
% - remove simulate (outside fun)
%- Help to choose bins (30+ obs p.bin check) [later]
%- update plots for bootstrapped
%- add lots of comments/documentation. [later]

%TODO CEC
%- Fit HT
%- Threshold diagnostic plot  (leave this until others are done).
%- simulate then transform back to original.
%- need to include poisson rate of occurence in simulatation ?...? 
%- Extend CEC code to do similar cross validation ?...smoothness?

%% Inputs to sim data
nDat=10000;
DrcEdg=[10,60,140,340]'; %NB> when nBin=1 Loc of DrcEdg doesn't matter...stationary
nBin=length(DrcEdg);
xi=-0.1;  %shape
sig=[2,3,1,1.8]';  %scale
PrmTrue=[xi;sig];
thr=zeros(nBin,1);  %threshold

%% Consuctor/simulate
if exist('DrcX','var') && exist('Y','var') %if already have data
    MM = MarginalModel(DrcX,Y,DrcEdg,PrmTrue); 
end
    %MM=MarginalModel();   %empty constructor
    %MM=SimulateMarginal(MM,nDat,DrcEdg,PrmTrue,thr); %Simulate data if don't have any    
%end


%% Generalised Pareto: fit for Exc
Tau=0.8;    %non exceedence probability
nBt=1;  %if want to bootstrap, make nBt > 1

%starting param ------
p0=[0;ones(nBin,1)*2];  
%CV settings --------
CVAllBt=0; %=1 then CV smoothness for every bootstrap resample (slow), =0 use opt. smoothness from orig dataset CV only 
nSigLam=10;  %no. smoothness penalties tried in CV (limit variation across neighbouring bins)
SigLam=logspace(-8,4,nSigLam); %try range smoothness penalties for sigma varying by bin
nCV=10; %no. cross-validation groups

%Fit ----------------
if 0
    PlotDiagnGP(MM) %diagnostic plots !! TO FIX FOR MULTI BOOT
end

%Bootstrap -----
MM = Fit(MM,Tau,nBt,CVAllBt,SigLam,p0,nCV);    %TODO merge FitMargModel and BootMargModel
%PlotData(MM)


%% PIT: Transform to Gumbel margins using probability integral transform
MM = GumbelMargins(MM);
PlotData(MM);




