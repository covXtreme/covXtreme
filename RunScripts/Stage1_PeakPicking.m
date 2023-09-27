% Copyright © [2023] Shell Global Solutions International B.V. All Rights Reserved.
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

%% Stage1_PeakPicking
%% Peak picking if we are reading in observed data

clear; clc; close all;
addpath('../Code');

%% Create figure directory
if ~exist('Figures','dir')
    mkdir('Figures')
end

load('CNS_mo_response')

%% creating figure directory
if ~exist('Figures','dir')
    mkdir('Figures')
end

%% Parameters
RspLbl={'maxPitch','Tp'}; %main and associated variable labels
CvrLbl={'Direction'}; %Covariate labels
%CvrLbl={'Direction'}; %Covariate labels

Rsp=vessel.maxPitch;  %main response
Cvr=wave.dm;  %covariates (e.g. direction, season, location)
IsPrdCrv=1;  %(boolean) vector flag for if covariate is periodic. 

%InlineCurr= cos(curr.dirn-wave.dm).*curr.spd; %current projected onto wave direction 
Asc=wave.tp_sea;  %associated variable(s)   wave.tp_sea  

[Rsp,Cvr,Asc] = PreProcessData(Rsp,Cvr,Asc, vessel.heading);

NEP=0.7; %non-exceedence quantile level used to get peak picking threshold 

%% Peak Picking
Dat=PeakPicking(Rsp,Cvr,Asc,IsPrdCrv,NEP,RspLbl,CvrLbl);

%% Save
if ~exist('Output','dir')
    mkdir('Output')
else
    F=cellstr(ls('Output/*.mat'));
    if numel(F)>0
       warning('Existing files in output directory') 
    end
end
save('Output\Data','Dat');
