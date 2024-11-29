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

%% Stage6_StormTrajectorySimulation
% load peak picked data set
load('Output/Data','Dat')

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

%load Heffernan and Tawn model
load('Output/HT','HT');

%% Storm trajectory options
%% Run storm trajectory simulation
StrTrj=StormTrajectorySimulation(Dat,Mrg,HT);

%% Save result
if ~exist('Output','dir')
   mkdir('Output') 
end
save('Output/StrTrj','StrTrj');



%% Using marginal model for dominant variate only
iDmn=1; %Always iDmn=1 for dominant variate

load(sprintf('Output/MM%g',iDmn),'MM');

% Simulate event set
Sml=cell(nSml,1);
for iS=1:nSml
    nStr = poissrnd(MM.nDat*RtrPrd/MM.Yrs);
    Sml{iS}=sample_MC(MM,nStr); %This is input to storm matching
end

% Diagnostic plot
load('Output/Data','Dat');
clf; hold on;
for iS=1:nSml
    tSml=sample_MC(MM,MM.nDat);
    plot(sort(Dat.Y(:,1)),sort(tSml.Org),'k.');
end
xlabel('Sorted originial sample');
ylabel('Sorted simulated sample');

%% Using HT model for all variables

% Load HT model
load('Output/HT','HT');
smlMethod="randX";

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

% Simulate event set
SmlHT=cell(nSml,1);
fprintf(1,'HT simulation: ');
for iS=1:nSml
    nStr = poissrnd(Mrg(1).nDat*RtrPrd/Mrg(1).Yrs);
    SmlHT{iS}=SimulateMC(HT,Mrg,smlMethod,nStr); %This is input to storm matching
end

% Diagnostic plot
load('Output/Data','Dat');
tSml=cell(nSml,1);
fprintf(1,'HT simulation: ');
for iS=1:nSml
    tSml{iS}=SimulateMC(HT,Mrg,smlMethod,Mrg(1).nDat); 
end
clf;
for iDmn=1:HT.nDmn
    subplot(1,HT.nDmn,iDmn); hold on;
    for iS=1:nSml
        plot(sort(Dat.Y(:,iDmn)),sort(tSml{iS}.Org(1,:,iDmn)),'k.');
    end
    xlabel('Sorted originial sample');
    ylabel('Sorted simulated sample');
    title(Dat.RspLbl(iDmn));
end %iDmn

%% Naive matching

load('Output/Data','Dat');

load('Output/Bin','Bn');

% Loop over storm peaks in Sml

% For each storm peak, identify the nearest nNgh storms in Dat.StrTrj
% This will use 
% - only dominant variate for "Kevin's" method, 
% - all variates for HT method
% Only match with a storm which has its storm peak in the same covariate bin

% Adjust matched storm trajectory based on physical constraints (do later)

% Save the full matched storm trajectory
% - bin index
% - nAsc+1 values of dominant and associated variates


%% Create storm trajectory bin allocations
fprintf(1,'Calculating bin allocations for historical trajectories\n');
for iS=1:Bn.n

    tCvr=zeros(Bn.n,1);
    tStrTrjCvr=Dat.StrTrj.Cvr{iS,:};
    tCvr(1:size(tStrTrjCvr,1))=tStrTrjCvr;
    iVrb = 0;
    tBn=BinAllocation(Bn,tCvr, iVrb);
    Dat.StrTrj.A{iS,:}=tBn.A(1:size(tStrTrjCvr,1));

end



