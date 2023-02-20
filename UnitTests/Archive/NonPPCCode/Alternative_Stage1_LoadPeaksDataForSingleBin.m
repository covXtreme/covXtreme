%% "Stage1_LoadPeaksData" 
%- If already have peak-picked data Dat.mat, load it and copy as "Datamat"
% into "Output" folder so can run analysis straight from Stage2.

% - If want to replicate Shiraz's analysis (split data into bins first,
% then fit sep. HT model to each bin), set allScts = 0 and pick out a
% sector you want to run analysis for (iQud). This will clear all data
% outside that sector. NOTE: this use of the code is a bit of a hack, but
% did help to verify similar results to Shiraz's analysis on Surge

clear; clc;
addpath('R:/Emma.Ross/Extremes/Devel_Code/CE/Code');
allScts = 0;

%Copy data into code structure
if ~exist('Output','dir')
    mkdir('Output')
end

if ~exist('Figures','dir')
    mkdir('Figures')
end

if allScts %fit for all sectors in one go
   load('Dat.mat')  %assumes you've already saved in "Data" format
%    t=Dat.Y(:,2)>inf;
%    Dat.X=Dat.X(t==0,:);
%    Dat.Y=Dat.Y(t==0,:);
   save('./Output/Data.mat','Dat')
   %copyfile('Dat.mat','Output/Data.mat')  %assumes you've already saved in "Data" format
else %split sector, fit for iSct = ?
   DrcQud=[0,180,270]';
   iQud = 1;
   load('Dat.mat')  %assumes you've already saved in "Data" format
   A=BinAllocation(Dat.X,DrcQud,Dat.Y,true,Dat.Lbl);
   Dat.X=Dat.X(A==iQud);
   Dat.Y=Dat.Y(A==iQud,:);
   save('./Output/Data.mat','Dat')
end


% SHOULD KNOW THIS:
% NEP=0.6; %non-exceedence quantile level used to peak pick Srg data

