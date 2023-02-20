%% Stage5_Contour
%% Use output of Heffernan and Tawn model to estimate environmental contours

%% Clean up
clear; clc; close all; 

%% Load marginal model
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

%% Load Heffernan and Tawn model
load('Output/HT','HT');

%% Define contour options
OptCnt=OptionsContours; %initiate OptCnt with default values; type OptionsContours for default values
%
% Optional edits to control contour estimation
% The values below have been chosen to work reasonably for Hs-Tp type conditions
%
nSml=1e6;       %Number to SimuLate: number of importance samples to use
nGrd=200;       %Number of GRid points for each of x and y when rectangular gridding needed
nPnt=200;       %Number of PoinTs on the contour
SmtWdtC=5;      %Smoothing Width for the huseby contour C function (see Huseby et al 2015)
BndWdtScl=0.02; %Band with scale for HTDns contour
%
% Contour method options are
%'Exc'    constant exceedence probability contour
%'HTDns'  constant density contour
%'HusOld' Huseby local tangent contour (original)
%'Hus'    Huseby local tangent contour (cleaned of "box tie" effects)
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