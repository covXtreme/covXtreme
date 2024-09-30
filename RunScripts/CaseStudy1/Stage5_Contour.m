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

%% Stage5_Contour
%%%%%%%%%%%%%%%%%%%%%%%%%
%load marginal model
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
%load Heffernan and Tawn model
load('Output/HT','HT');

%% Define contour options
OptCnt=OptionsContours; %initiate OptCnt with default values; type OptionsContours for default values
%
% Optional edits to control contour estimation
% The values below have been chosen to work reasonably for Hs-Tp type conditions
%
OptionsContours.nSml=1e6;       %Number to SimuLate: number of importance samples to use
OptionsContours.nGrd=200;       %Number of GRid points for each of x and y when rectangular gridding needed
OptionsContours.nPnt=200;       %Number of PoinTs on the contour
OptionsContours.SmtWdtC=5;      %Smoothing Width for the Huseby contour C function (see Huseby et al 2015)
OptionsContours.BndWdtScl=0.05; %Band width scale for HTDns contour
%
% Contour method options are
%'Exc'    constant exceedence probability contour
%'HTDns'  constant density contour
%'HusOld' Huseby local tangent contour (original)
%'Hus'    Huseby local tangent contour (cleaned of "bow tie" effects)
%
%OptCnt.Mth={'Hus'}; %cell array of contour methods to be used 
%OptCnt.Mth={'Exc'}; %cell array of contour methods to be used 
%OptCnt.Mth={'Exc','HTDns','Hus','HusOld'};   %cell array of contour methods to be used 
OptCnt.Mth={'Exc','HTDns','Hus'};   %cell array of contour methods to be used 

%% Estimate contour
Cnt=Contour(HT,Mrg,OptCnt); %initiate contour object
Cnt=Cnt.makeContours(Mrg,HT); %estimate contours

%% Plot contours
Cnt.Plot(Mrg);

%% Save result
if ~exist('Output','dir')
   mkdir('Output') 
end
save('Output/Cnt','Cnt')


