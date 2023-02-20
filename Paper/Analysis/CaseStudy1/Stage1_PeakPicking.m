clear; clc; close all;
addpath('C:\Users\Philip.Jonathan\Git\MetOcean_PPC_CE\Code');

load('Dat')

%% Parameters
RspLbl={'Hs','T2'}; %main and associated variable labels
CvrLbl={'Direction'}; %Covariate labels

Rsp=X(:,1);  %main response
Cvr=X(:,3);  %covariates (e.g. direction, season, location)
IsPrdCrv=1;  %(boolean) vector flag for if covariate is periodic. 
Asc=X(:,2);  %associated variable(s)   wave.tp_sea  

NEP=0.7; %non-exceedence quantile level used to get peak picking threshold 

%% Peak Picking
Dat=PeakPicking(Rsp,Cvr,Asc,IsPrdCrv,NEP,RspLbl,CvrLbl);

figure(1); clf;
pPltDrcSct(Dat.X,Dat.Y,RspLbl,CvrLbl);


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
