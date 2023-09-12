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
