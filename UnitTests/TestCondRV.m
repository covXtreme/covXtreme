clc; clear; close all;
load('Output\HT')
%%%%%%%%%%%%%%%%%%%%%%%%%
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

iRtr=1;

nRls=1e3;
Mrg(1).RtrPrd=100;
HT=HT.ConditionalReturnValue(Mrg,nRls);




