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







