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

%% Stage2_BinData
% Choose covariate bins

%% Inputs 
load('Output/Data','Dat')

%% Choose bins for each covariate dimension
BinEdg={[25,230,275]',[145,270]'}; %In single covariate case: input in format {[]'}. In multiple covariate case: input in format {[]',[]'}.
%Note: in single covariate case: input in format {[]'}. In multiple covariate case: input in format {[]',[]'}.
%Note: covariate bins must be suitable for all dimensions (main and associated)                             
%Note: if using non-periodic covariate (you set 'IsPrdCvr' = 0 in Stage1), and you must provide the endpoint/maximum value of the covariate as last bin-edge

%% Allocate Data to Bins and Plot
Bn=CovariateBinning(Dat.X,BinEdg,Dat.IsPrd,Dat.Y,Dat.RspLbl,Dat.CvrLbl);

%% Save bin Edges
save('Output/Bin','Bn')
