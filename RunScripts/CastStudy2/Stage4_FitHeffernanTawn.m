clear; clc; close all;

%% Stage4_FitH&T

FM=cellstr(ls('Output/MM*.mat'));
for iM=1:numel(FM) %load all marginal models
    load(sprintf('Output/MM%g',iM),'MM')  %Load margin 2
    if iM==1
        Mrg=MM;
    else
        Mrg=cat(1,Mrg,MM);
    end
    clear MM;
end

%% Fit Heffernan and Tawn model options
HTNEP=[0.8,0.9];  
%Note: HTNEP
%- conditional extremes non-exceedance threshold interval
NonStationary=[true,false,false,false];  
%Note: NonStationary
%- flag for whether alpha, beta, mu and sigma are non-stationary with covariate
%- if true, corresponding parameter assumed non-stationary, with covariate bins from marginal analysis
SampleLocalResid = false; 
%Note: SampleLocalResid
%- if true: when simulating under H&T model, resample residuals locally from same covariate sector; 
%- if false, resample residuals from any sector
%- if you have any bins with very few observations; set this to false.

%% Cross Validation defaults (Optional can specify some or all of these)
CV.CVMth=0;     %0 Only Cross Validate smoothness for original dataset (fast);
                %1 Cross Validate smoothness for every bootstrap resample (slow),
CV.nCV=10;      %number cross-validation groups
CV.nSmth=10;    %number smoothnesses tried in CV
CV.SmthLB=-8;   %lower bound (log10) for smoothness range
CV.SmthUB=8;   %upper bound (log10) for smoothness range

%% Fit model
HT=HeffernanTawn(Mrg,HTNEP,NonStationary,CV,SampleLocalResid); %Fit Heffernan & Tawn joint exceedance model
HT=HT.Fit(Mrg);

%% Save result
if ~exist('Output','dir')
   mkdir('Output') 
end
save('Output/HT','HT')

Plot(HT,Mrg)







