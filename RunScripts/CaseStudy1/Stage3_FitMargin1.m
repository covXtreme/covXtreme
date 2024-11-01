clear; clc; close all;
%% Stage3_FitMarginal Generalised Pareto
% Fit Piecewise constant extreme value model for margin / dimension 'iDmn'
% 
% - Piecewise constant threshold using local bin quantile given NEP range
% - Annual rate-of-occurence using local Poisson max-likelihood estimate
% for each bin: (No.observations)/(No. years of data)
% - Generalised Pareto fitted above the threshold
%    = Shape is assumed constant acorss all covariate bins
%    = Scale is non-stationary (varies by bin) but smoothed across bins, with smoothness penalty chosen using cross validation
% - Transformation of data to standard-margins: can choose Laplace or
% Gumbel standard margins

%% Inputs
load('Output/Data','Dat')  %load data from Stage 1 
load('Output/Bin','Bn')  %load bins from Stage 2
%% Dimension to fit
%
iDmn=1;  %which dimension to fit marginal model to
%
NEP=[0.2,0.9];  %GP non exceedence probability range
nB=100;   %number bootstrap resamples
Yrs=34;  %number of years of data
RtrPrd=[10,100]; %vector of return Periods

%% Cross Validation defaults (Optional can specify some or all of these)
CV.CVMth=0;     %0 Only Cross Validate smoothness for original dataset (fast);
                %1 Cross Validate smoothness for every bootstrap resample (slow),
CV.nCV=10;      %number cross-validation groups
CV.nSmth=10;    %number smoothnesses tried in CV
CV.SmthLB=-4;   %lower bound (log10)  for smmothness range
CV.SmthUB=4;   %upper bound (log10)  for smmothness range

%% Transformation to Standard-margins 
MarginType='Laplace'; %Laplace or Gumbel

%% Fit model
MM=MarginalModel(Dat,iDmn,NEP,Bn,nB,Yrs,RtrPrd,CV,MarginType);

%% Save result
if ~exist('Output','dir')
   mkdir('Output') 
end
save(sprintf('Output/MM%g',iDmn),'MM')

