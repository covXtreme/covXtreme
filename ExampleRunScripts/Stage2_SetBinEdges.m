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

clear; clc; close all;
%% Stage2_BinData
% Stage to choose covariate (directional) bins

%% Inputs 
load('Output/Data','Dat')

%%Choose bins cell array ones for each covariate dimension
BinEdg={[0,180]'}; %NB. Must be suitable for both Margins if nMarg > 1


%% Allocate Data to Bins and Plot
Bn=CovariateBinning(Dat.X,BinEdg,Dat.IsPrd,Dat.Y,Dat.RspLbl,Dat.CvrLbl);

%% Save bin Edges
save('Output/Bin','Bn')
