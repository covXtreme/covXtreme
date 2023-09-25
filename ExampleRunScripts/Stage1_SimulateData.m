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

clear; clc;

%% Step1_SimulateData.m: CEC or Marginal

% Simulate marginal data from a single covariate (directional) extreme value model in nDmn
% dimensions.
% Data simulated as piecewise constant model defined over covariate (directional) bins.
% Model has Generalised Pareto margins with constant shape parameter and
% scale parameter which varies by covariate (directional) bin.

%% Create figure directory
if ~exist('Figures','dir')
    mkdir('Figures')
end

%% Inputs 
Mdl.nDmn= 2;  %number of dimensions/responses (currently 1D or 2D) 
Mdl.nDat=400; % 2e6;  %number of observations to simulate
Mdl.nBin=1;     %number of covariate bins (Dimension of DrcEdg, Mdl.MM(iDmn).Scl, Mdl.MM(1).Thr and Mdl.Rat must match)
Mdl.DrcEdg=[0]';%[30:30:240]'; %Location of covariate (directional) bin edges

%% Marginal Parameters
% For each margin, need to specify:
% - Generalised Pareto shape 1 x 1
% - Generalised Pareto scale nBin x 1
% - Generalised Pareto threshold nBin x 1
%Margin 1
Mdl.MM(1).Shp=-0.2;  %shape  
Mdl.MM(1).Scl=[1]';  %scale
Mdl.MM(1).Thr=zeros(Mdl.nBin,1);  %threshold
%Margin 2
if Mdl.nDmn>1
    Mdl.MM(2).Shp=-0.2;  %shape
    Mdl.MM(2).Scl=[1]';  %scale
    Mdl.MM(2).Thr=zeros(Mdl.nBin,1);  %threshold
end
% Poission rate of occurrence nBin x 1
Mdl.Rat=[10]';  %(relative) rate of occurence (common across all margins)

%% Joint parameters: specify joint dependency model if simulating 2 responses

if Mdl.nDmn>1
    % choice of joint dependency model
    Mdl.Jnt.Mth ='LGS';  % Options:
                         % - 'MVN': multivairate normal with dependency
                         % parameter rho
                         % - 'LGS': logistic with single parameter alp
                         % - 'ASL': Assymetric Logistic
    
    % Set parameters of joing dependency model                     
    switch Mdl.Jnt.Mth
        case 'MVN'  %multivariate normal N(0,Sig)
            
            %Sig = [1,rho; rho,1];
            Mdl.Jnt.Rho=0.8;%[(0.5+(0.8.*(cos(2.*pi*[30:19:320]./360)+1)/2))/1.5];%[0.15:0.05:0.9];%[0.1 0.5 0.9];  %correlation on [0,1];
        case 'LGS' %logistic with single parameter alp on [0,1]
            Mdl.Jnt.Alp=[0.1,0.9,0.6,0.4];%;%[0.15:0.05:0.9];
            %Mdl.Jnt.Alp=0.2';
            %periodic values of alpha that are  close to zero and 1 [(cos(2.*pi*[30:19:320]./360)+1)/2]
            %periodic values of alpha that are not as close to zero and 1 [(cos(2.*pi*[30:19:320]./360)+1.5)/3]
        case 'ASL' % Assymetric Logistic (misture of logistic and random)
            Mdl.Jnt.Alp=0.1;  %logistic single parameter alp on [0,1]
            Mdl.Jnt.Theta=[0,0.5]; %mixture probabilities for each margin both on [0,1]
                                 %...(proportion of 'random' data off the logistic relationship)      
    end
end



%%  Simulate
Dat=SimulateData(Mdl);
%% Save 
if ~exist('Output','dir')
   mkdir('Output') 
end
if Mdl.nDat>1e4
    save('Output\DataLarge','Dat','Mdl');
else
    save('Output\Data','Dat','Mdl');
end


%%
% XI=lhsdesign(100000,2).*max(Dat.Y);
% 
% F = ksdensity(Dat.Y,XI);
% 
% clf;
% subplot(2,1,1)
% plot(Dat.Y(:,1),Dat.Y(:,2),'k.')
% hold on
% 
% subplot(2,1,2)
% scatter(XI(:,1),XI(:,2),10,F,'filled')
% colormap('jet')

