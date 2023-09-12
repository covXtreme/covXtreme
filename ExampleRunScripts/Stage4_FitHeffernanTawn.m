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

%% Stage4_FitH&T
%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% Fit heffernan and Tawn model
HTNEP=[0.7,0.9];  %Conditional Extreme NEP

NonStationary=[true,false,false,false];  %flag for whether to estimate the HT alpha parameter in a non -stationary way using the same bins as the Marginal analysis


% Cross Validation defaults (Optional can specify some or all of these)
CV.CVMth=0;     %0 Only Cross Validate smoothness for original dataset (fast);
                %1 Cross Validate smoothness for every bootstrap resample (slow),
CV.nCV=10;      %number cross-validation groups
CV.nSmth=10;    %number smoothnesses tried in CV
CV.SmthLB=-4;   %lower bound (log10)  for smmothness range
CV.SmthUB=4;   %upper bound (log10)  for smmothness range


SmpLclRsdOn=false;
%% Fit model
HT=HeffernanTawn(Mrg,HTNEP,NonStationary,CV,SmpLclRsdOn);
HT=HT.Fit(Mrg);

%% Make plots of results         
Plot(HT,Mrg)
         
%% Save result
if ~exist('Output','dir')
   mkdir('Output') 
end
save('Output/HT','HT')

            
