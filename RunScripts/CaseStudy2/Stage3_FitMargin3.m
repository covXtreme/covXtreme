% Copyright © [2023] Shell Global Solutions International B.V. All Rights Reserved.
%
% SPDX-License-Identifier: Apache-2.0-or-later
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
% Fit Piecewise constant extreme value model for margin / dimension 'iDmn'
% 
% - Piecewise constant threshold using local bin quantile given NEP range
% - Annual rate-of-occurence using local Posson max-likelihood estimate
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
iDmn=3;  %which dimension to fit
%
NEP=[0.7,0.85];  %GP non exceedence probability range
nBoot=100;   %number bootstrap resamples
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
MarginType='Laplace';  %Laplace or Gumbel

%% Fit model
MM=MarginalModel(Dat,iDmn,NEP,Bn,nBoot,Yrs,RtrPrd,CV,MarginType); 

%% Save result
if ~exist('Output','dir')
   mkdir('Output') 
end
save(sprintf('Output/MM%g',iDmn),'MM')

