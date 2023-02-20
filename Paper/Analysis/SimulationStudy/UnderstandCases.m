load Dsg

%% Understanding the design
%AnlStt is 84 x 3

%Column 1 has values 1:14 referring to the covariate 
%These are defined in the BinEdg cell array
%BinEdg is 14 x 1 cell array with either one or two cell entries, giving
%the covariate bin boundaries for either 1-D or 2-D covariates

%Column 2 has values 1:3 referring to the set of responses
%These are defined in the IdxRsp cell array
%IdxRsp is 3 x 1 cell array with either 4 of 5 entries, given list of
%responses, the first being the conditioning variate
%***So there are no examples conditioned on structural response?????

%Column 3 has values 1 or 2, referring to the NEP range used for the CndExt
%model
%The two ranges are defined in 2x2 numeric array HTNEPRng (1 is lower than
%2)
%PksNEP = 0.7 is the fixed range of peak-picking NEP
%MrgNEPRng = [0.7,0.85] is the fixed range used for marginal and joint modelling

%% BinEdg
% %Factors 1 and 2: Covariate Number of Bins (hard coded bin edge choice, variable nBins)
% nCvrStt = 14;
% BinEdg = cell(nCvrStt,1); %structure to store bin-edges for each test case. col1 = Drc, col2 = Ssn by degrees, col3 = Ssn by month
% IsPrd = cell(nCvrStt,1); %structure to store bin-edges for each test case. col1 = Drc, col2 = Ssn by degrees, col3 = Ssn by month
% %-- Directional
% BinEdg{1} = {0}; IsPrd{1} = 1; IdxCvr{1} = 1; %1 bin
% BinEdg{2} = {[20,160]'}; IsPrd{2} = 1; IdxCvr{2} = 1; %2 bins
% BinEdg{3} = {[0,20,160,270,315]'}; IsPrd{3} = 1; IdxCvr{3} = 1; %5 bins
% BinEdg{4} = {[0,20,60,115,150,200,220,240,270,315]'};  IsPrd{4} = 1; IdxCvr{4} = 1; % 10 bins
% %-- Seasonal
% % ....NonPeriodic( on  0 to 366, i.e. ~days)
% BinEdg{5} = {[0,366]'}; IsPrd{5} = 0; IdxCvr{5} = 3; %1 bin
% BinEdg{6} = {[0,100,250,366]'}; IsPrd{6} = 0; IdxCvr{6} = 3; %3 bins
% BinEdg{7} = {[0,70,160,270,366]'}; IsPrd{7} = 0; IdxCvr{7} = 3; %4 bins
% % ....Periodic (on 0 to 360)
% BinEdg{8} = {0}; IsPrd{8} = 1; IdxCvr{8} = 2;%1 bin
% BinEdg{9} = {[100,250]'}; IsPrd{9} = 1; IdxCvr{9} = 2; %2 bins
% BinEdg{10} = {[70,100,250,290]'}; IsPrd{10} = 1; IdxCvr{10} = 2; %4 bins
% %--Directional Seasonal
% % ...both periodic
% BinEdg{11} = {0,0}; IsPrd{11} = [1,1]; IdxCvr{11} = [1,2]; %both have 1 bin
% BinEdg{12} = {[0,20,160,270,315]',[100,250]'}; IsPrd{12} =[1,1]; IdxCvr{12} = [1,2];  %5 directional, 2 seasonal bins
% % ...season non-periodic
% BinEdg{13} = {0,[0,366]'}; IsPrd{13} = [1,0]; IdxCvr{13} = [1,3];   %both have 1 bin
% BinEdg{14} = {[0,20,160,270,315]', [0,100,250,366]'}; IsPrd{14} = [1,0]; IdxCvr{14} = [1,3];   %5 directional, 3 seasonal bins
