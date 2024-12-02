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

%% Stage6_StormTrajectorySimulation
%% Load time-series data
load('..\CNS_mo_response')

%% Parameters
RspLbl={'Hs','Tp'}; %main and associated variable labels
CvrLbl={'Direction'}; %covariate labels
Rsp=wave.hs;  %main response data
Cvr=[wave.dm];  %covariate data (e.g. direction, season, location)
IsPrdCvr=1;  %flag dictating periodicity of covariate(s). If 1, covariate data loops on 360. When have >1 covariate, this is a vector input e.g. [1,0].   
Asc=[wave.tp_sea];  %associated variable(s)   
NEP=0.7; %non-exceedence quantile level used to set peak picking threshold
IsStrTrj = 1; % Set to 1 to allow the use of storm trajectories.

% load peak picked data set
load('Output/Data','Dat')
load('Output/Bin','Bn')  %load bins from Stage 2
% Load and combine marginal models
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
%% Storm trajectory options

%% Run storm trajectory simulation
StrTrj=StormTrajectorySimulation(Rsp,Asc,Cvr,Dat,Bn,Mrg);

%% Save result
if ~exist('Output','dir')
   mkdir('Output') 
end
save('Output/StrTrj','StrTrj');

