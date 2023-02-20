clear; clc; close all;

% addpath('R:\Emma.Ross\Extremes\Devel_Code\CE\Code');
% addpath('R:\Emma.Ross\Extremes\Devel_Code\E_Code\');

addpath('C:\PhilipGit\Metocean_CEVA_E_code');
AnlDrc = 'Z:\project\MetOcean\PPC\ResponseCasesPPC\AutoRunCases';
addpath(AnlDrc)

%TODO: Stuff to fix
% - Make coloured scatter version of the illustration of empirical density
% plot (more accurate and helps explain in a hand-wavey what the density
% and importance sampling methods are doing)

%% Cases covered in test runs:

% Factor1 and2: Bin choice and covariates [14 levels]
% - Direction: 1,2,5,10 bins (pre-set), Season: 1,3,4 bins
% - Drc (4 levels), Snn (3 levels), Ssn (3 levels -> on-periodic) , Drc-Ssn (2 levels -> both periodic),
% ...Drc-Ssn (2 levels -> direction periodic + season nonperiodic)

% Factor3: Associated Variables [3 levels]
% - [3, 3 shuffled, 4], with pre-chosen input-data as below (we've done 1 and 2D
% cases loads of times...)
% 3D = Tp,WS,Crr (idx 2,3,4) | Hs (idx 1)
% 3D shuffled = Hs,Tp,WS (idx 4,2,3) | Crr (idx 1)   %!! Hs-->Crr UPDATE
% 4D = Tp,WS,Crr,Otm  (2,3,4,5) | Hs (idx 1)

% Factor4: NEP choices [1 level]
% - Fixed to 0.7 for peak-picking; Fixed to [0.7-0.85] for Marginal
% - HT NEP: one range below the marginal ([0.6,0.7]) and one above ([0.85,0.95])

% Factor5: Marginal Whitening Disrtibution:  [1 level]
% - For now fix to Laplace only

% Factor6: Inside each test case: HT non-stationary / stationary [2 sub-levels]
% ---> only redo HT part onwards!

% Factor7: Contour types
% CntMth={'Hus','Exc'};


%% Set up data

load('Z:\project\MetOcean\PPC\ResponseCasesPPC\CNS_mo_response')
AllRspLbl={'Hs','Tp','WS','Otm','Pitch'}; %main and associated variable labels  ,'Roll','Bsr',
AllRspDat = [wave.hs,wave.tp_sea,wind.spd,curr.spd,jacket.otm,vessel.maxPitch]; %store all response data together
%vessel.maxRoll, jacket.bsr
AllCvrLbl={'Direction','Season','Day'}; %Covariate labels
AllCvrDat=[wave.dm,E_SsnDgr(t),day(datetime(t,'ConvertFrom','datenum'),'dayofyear')];  %direction, season (0 to 360), season (0 to 366)

%Factors 1 and 2: Covariate Number of Bins (hard coded bin edge choice, variable nBins)
nCvrStt = 14;
BinEdg = cell(nCvrStt,1); %structure to store bin-edges for each test case. col1 = Drc, col2 = Ssn by degrees, col3 = Ssn by month
IsPrd = cell(nCvrStt,1); %structure to store bin-edges for each test case. col1 = Drc, col2 = Ssn by degrees, col3 = Ssn by month
%-- Directional
BinEdg{1} = {0}; IsPrd{1} = 1; IdxCvr{1} = 1; %1 bin
BinEdg{2} = {[20,160]'}; IsPrd{2} = 1; IdxCvr{2} = 1; %2 bins
BinEdg{3} = {[0,20,160,270,315]'}; IsPrd{3} = 1; IdxCvr{3} = 1; %5 bins
BinEdg{4} = {[0,20,60,115,150,200,220,240,270,315]'};  IsPrd{4} = 1; IdxCvr{4} = 1; % 10 bins
%-- Seasonal
% ....NonPeriodic( on  0 to 366, i.e. ~days)
BinEdg{5} = {[0,366]'}; IsPrd{5} = 0; IdxCvr{5} = 3; %1 bin
BinEdg{6} = {[0,100,250,366]'}; IsPrd{6} = 0; IdxCvr{6} = 3; %3 bins
BinEdg{7} = {[0,70,160,270,366]'}; IsPrd{7} = 0; IdxCvr{7} = 3; %4 bins
% ....Periodic (on 0 to 360)
BinEdg{8} = {0}; IsPrd{8} = 1; IdxCvr{8} = 2;%1 bin
BinEdg{9} = {[100,250]'}; IsPrd{9} = 1; IdxCvr{9} = 2; %2 bins
BinEdg{10} = {[70,100,250,290]'}; IsPrd{10} = 1; IdxCvr{10} = 2; %4 bins
%--Directional Seasonal
% ...both periodic
BinEdg{11} = {0,0}; IsPrd{11} = [1,1]; IdxCvr{11} = [1,2]; %both have 1 bin
BinEdg{12} = {[0,20,160,270,315]',[100,250]'}; IsPrd{12} =[1,1]; IdxCvr{12} = [1,2];  %5 directional, 2 seasonal bins
% ...season non-periodic
BinEdg{13} = {0,[0,366]'}; IsPrd{13} = [1,0]; IdxCvr{13} = [1,3];   %both have 1 bin
BinEdg{14} = {[0,20,160,270,315]', [0,100,250,366]'}; IsPrd{14} = [1,0]; IdxCvr{14} = [1,3];   %5 directional, 3 seasonal bins


%Factor 3: Number of Associated Vars
nRspStt = 3;
IdxRsp = cell(3,1);  %structure to store indices of responses for each test case
IdxRsp{1} =[1,2,3,4]';   %main variable first (one we condition on), associated listed after
IdxRsp{2} = [1,4,2,3]';
IdxRsp{3} = [1,2,3,4,5]';

%Factor 4: NEP choices [2 levels for HT, 1 for Marginal]
nHTNEPRng = 2;
HTNEPRng = [[0.6,0.7];
    [0.85,0.95]];

%Store all combinations of the above 3 factor levels
nTestCase =nCvrStt * nRspStt * nHTNEPRng;
AnlStt= NaN(nTestCase,3); %matrix storing combinations of settings for each run (indices), col1 = ...
tCol1 = repmat(1:nCvrStt,(nRspStt*nHTNEPRng),1);
tCol2 = repmat(1:nRspStt,nHTNEPRng,1);  tCol2b = repmat(tCol2(:),nCvrStt,1);
tCol3 = repmat(1:nHTNEPRng,nCvrStt*nRspStt,1)';
AnlStt(:,1) = tCol1(:);
AnlStt(:,2) = tCol2b;
AnlStt(:,3) = tCol3(:);

set(0,'DefaultFigureWindowStyle','docked')
%% Set up things which are common for all test cases
%(But still user-set inputs)

%Common for all cases:
% - (Factor4) NEP choices
PksNEP = 0.7;  %peak-picking NEP
MrgNEPRng = [0.7,0.85];  %used for marginal and joint modelling
% - (Factor5): Marginal Whitening Disrtibution:  [1 level]
MarginType='Laplace';
%- (Factor7): Contour types
CntMth={'Hus','Exc','HTDns'};
% - (Factor6): Inside each test case: HT non-stationary / stationary [2 sub-levels]


% TODO: This isn't coded in, just rerun this script, changing to stationary case (RunRunAllStages) ---> only redo HT part onwards!

VssHdn=vessel.heading;
save Dsg AnlStt BinEdg AllRspLbl IdxRsp IdxCvr ... 
    AllRspDat AllCvrDat AllCvrLbl IsPrd VssHdn ...
    AnlDrc PksNEP MrgNEPRng MarginType HTNEPRng CntMth;

return;

%% Loop over and run all cases
%In terminal:
%parpool(16)

TestCompleted = false(nTestCase,1);
failures =cell(length(nTestCase),1);

diary 'PPC_Log.out'

%Specifiy which stage to run
Stage=5;
for iTest = 30
    close all;
    try
        ParRunTests(AnlStt,iTest,BinEdg,AllRspLbl,IdxRsp,IdxCvr,AllRspDat,AllCvrDat,AllCvrLbl,IsPrd,vessel.heading,AnlDrc,PksNEP,MrgNEPRng,MarginType,HTNEPRng,CntMth,nTestCase,Stage)
        TestCompleted(iTest) = true;
    catch ME
       fprintf('Case %d failed\n',iTest)
       fprintf('%s\n',getReport(ME))
       failures{iTest}=ME;
    end
    
end

%delete(gcp('nocreate'))

diary off

%fail_cell=cell2char(failures(cellfun(@length,failures)~=0))';
%uni_fail_cell=unique(fail_cell);

%cases 68, 72, 82 - histc
%casese 84 - error 

