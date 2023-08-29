% Copyright 2023 covXtreme
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%      http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

clear; clc; close all;
%% Stage3_FitMarginal Generalised Pareto
%%%%%%%%%%%%%%%%%%%%%%%%%

%Fit Piecewise constant EVA Model
%
% Estimates piecewise constant
% Threshold using local bin quantile given NEP range
% Annual Rate using local Possion mle for each bin (Number of
% observations)/Number of Years of data
%
% Generalised Pareto Above the threshold
%     Shape is assumed constant acorss all bins
%     Scale is local but smoothed across bins. Smoothness can be
%     estimated using crossvalidation.
%% Inputs
load('Output/Data','Dat')  %load data from Stage 1 
load('Output/Bin','Bn')  %load bins from Stage 2
%% Dimension to fit
%
iDmn=2;  %which dimension to fit 
%
NEP=[0.5,0.7];  %GP non exceedence probability range
nB=50;   %number bootstrap resamples
Yrs=34;  %number of years of data
RtrPrd=[1,10,100,1000,10000,100000]; %vector of return Periods

%% Cross Validation defaults (Optional can specify some or all of these)
CV.CVMth=0;     %0 Only Cross Validate smoothness3 for original dataset (fast);
                %1 Cross Validate smoothness for every bootstrap resample (slow),
CV.nCV=10;      %number cross-validation groups
CV.nSmth=10;    %number smoothnesses tried in CV
CV.SmthLB=-4;   %lower bound (log10)  for smoothness range
CV.SmthUB=4;   %upper bound (log10)  for smoothness range

%% Optional Margin type Laplace or Gumbel
MarginType='Laplace';

%% Fit model
MM=  MarginalModel(Dat,iDmn,NEP,Bn,nB,Yrs,RtrPrd,CV,MarginType); %initialse PPC marginal model clas 

%% Save result
if ~exist('Output','dir')
   mkdir('Output') 
end
save(sprintf('Output/MM%g',iDmn),'MM')

