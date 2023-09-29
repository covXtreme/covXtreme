clear; clc; close all;
%% Stage2_BinData
% Choose covariate bins

%% Inputs 
load('Output/Data','Dat')

%% Choose bins for each covariate dimension
BinEdg={[0,20,60,225,270,315]'}; 
%Note: in single covariate case: input in format {[]'}. In multiple covariate case: input in format {[]',[]'}.
%Note: covariate bins must be suitable for all dimensions (main and associated)                             
%Note: if using non-periodic covariate (you set 'IsPrdCvr' = 0 in Stage1), and you must provide the endpoint/maximum value of the covariate as last bin-edge

%% Allocate Data to Bins and Plot
Bn=CovariateBinning(Dat.X,BinEdg,Dat.IsPrd,Dat.Y,Dat.RspLbl,Dat.CvrLbl);

%% Save bin Edges
save('Output/Bin','Bn')
