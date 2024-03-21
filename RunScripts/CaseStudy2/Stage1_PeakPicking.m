clear; clc; close all;

%% Stage 1: Peak-picking
% Extract peaks data 

load('..\CNS_mo_response')

%% Parameters
RspLbl={'OTM','Hs','WS'}; %main and associated variable labels
CvrLbl={'Direction','Season'}; %Covariate labels

Rsp=jacket.otm;  %main response

tSsn = SsnDgr(t); %convert date to season (on [0,360]) 
Cvr=[wave.dm,tSsn];  %covariate data (e.g. direction, season, location)
IsPrdCrv=[1,1];   %flag dictating periodicity of covariate(s). If 1, covariate data loops on 360. When have >1 covariate, this is a vector input e.g. [1,0].   

Asc=[wave.hs,wind.spd];  %associated variable(s) 

NEP=0.7; %non-exceedence quantile level used to set peak picking threshold 

%% Peak Picking

if ~exist(fullfile(cd,'Figures'),'dir')
    mkdir('Figures')
end
    
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
