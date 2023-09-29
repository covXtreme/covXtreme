clear; clc; close all;
%% Stage4_FitH&T
% Fit Heffernan & Tawn conditional extremes model

FM=cellstr(ls('Output/MM*.mat'));
for iM=1:numel(FM) %load all marginal models
    load(sprintf('Output/MM%g',iM),'MM');
    if iM==1
        Mrg=MM;
    else
        Mrg=cat(1,Mrg,MM);
    end
    clear MM;
end

%% Fit Heffernan and Tawn model options
HTNEP=[0.7,0.85];  
%Note: HTNEP
%- conditional extremes non-exceedance threshold interval
NonStationary=[true,false,false,false];  
%Note: NonStationary
%- flag for whether alpha, beta, mu and sigma are non-stationary with covariate
%- if true, corresponding parameter assumed non-stationary, with covariate bins from marginal analysis
SampleLocalResid = true; 
%Note: SampleLocalResid
%- if true: when simulating under H&T model, resample residuals locally from same covariate sector; 
%- if false, resample residuals from any sector
%- if you have any bins with very few observations; set this to false.

%% Cross Validation defaults (Optional can specify some or all of these)
CV.CVMth=0;     %0 Only Cross Validate smoothness for original dataset (fast);
                %1 Cross Validate smoothness for every bootstrap resample (slow),
CV.nCV=10;      %number cross-validation groups
CV.nSmth=10;    %number smoothnesses tried in CV
CV.SmthLB=-3;   %lower bound (log10)  for smoothness range
CV.SmthUB=2;   %upper bound (log10)  for smoothness range

%% Fit model
HT=HeffernanTawn(Mrg,HTNEP,NonStationary,CV,SampleLocalResid); %Fit Heffernan & Tawn joint exceedance model
HT=HT.Fit(Mrg);

%% Save result
if ~exist('Output','dir')
   mkdir('Output') 
end
save('Output/HT','HT')

%% Plot results
Plot(HT,Mrg)







