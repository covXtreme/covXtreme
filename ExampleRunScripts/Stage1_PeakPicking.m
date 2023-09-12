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
