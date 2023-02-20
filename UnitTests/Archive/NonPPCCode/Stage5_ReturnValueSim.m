%% Surge Analysis
% Stage 5: Return value simulations
addpath('R:/Emma.Ross/Extremes/Devel_Code/CE/Code')
addpath('R:/Emma.Ross/Extremes/CE_Surge')

%% Load data from marginal and HT analysis
cd('R:\Samta.Sam\Extremes\CE_Surge\Analysis\Anga\Mdn_Stn')

load(sprintf('Output/MM%g',1),'MM') %load and concat margin 1 and 2
Mrg=MM;
load(sprintf('Output/MM%g',2),'MM')
Mrg=cat(1,Mrg,MM); clear MM;
load('Output/HT','HT'); %load Heffernan and Tawn model

%% Simulation inputs
%inputs
nSim = 1000; %no. of repetitions/simulations
RtrPrd = 100; %return period in years  (run time = 40s for 100yr, 7 mins for 1000yr)

%% Simulate
tic
RV = ReturnValueSim(HT,Mrg,RtrPrd,nSim,'Output');
toc
%% Plot boxplots of RVs
BxpltLayout= 0; 
ReturnValuePlotsPPC(RV,BxpltLayout,'./Figures')
