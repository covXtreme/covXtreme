clear; clc;
addpath('R:\Ross.Towe\Extremes\Metocean\Code\CE\Code');

%% Step1_SimulateData.m: CEC or Marginal

% Simulate marginal data from a single covariate (directional) extreme value model in nDmn
% dimensions.
% Data simulated as piecewise constant model defined over covariate (directional) bins.
% Model has Generalised Pareto margins with constant shape parameter and
% scale parameter which varies by covariate (directional) bin.

%% Inputs 
Mdl.nDmn= 2;  %number of dimensions/responses (currently 1D or 2D) 
Mdl.nDat=3000;  %number of observations to simulate
Mdl.nBin=3;     %number of covariate bins (Dimension of DrcEdg, Mdl.MM(iDmn).Scl, Mdl.MM(1).Thr and Mdl.Rat must match)
Mdl.DrcEdg=[30,100,240]'; %Location of covariate (directional) bin edges
Mdl.Jnt.Mth='MVN';



%% Marginal Parameters
% For each margin, need to specify:
% - Generalised Pareto shape 1 x 1
% - Generalised Pareto scale nBin x 1
% - Generalised Pareto threshold nBin x 1
%Margin 1
Mdl.MM(1).Shp=-0.1;  %shape  
Mdl.MM(1).Scl=[8;7;8];  %scale
Mdl.MM(1).Thr=zeros(Mdl.nBin,1);  %threshold
%Margin 2
if Mdl.nDmn>1
    Mdl.MM(2).Shp=-0.3;  %shape
    Mdl.MM(2).Scl=[0.1;0.2;0.3];  %scale
    Mdl.MM(2).Thr=zeros(Mdl.nBin,1);  %threshold
end
% Poission rate of occurrence nBin x 1
Mdl.Rat=[1,1,1]';  %(relative) rate of occurence (common across all margins)

%% Joint parameters: specify joint dependency model if simulating 2 responses

if Mdl.nDmn>1
    % choice of joint dependency model
    Mdl.Jnt.Mth ='ASL';  % Options:
                         % - 'MVN': multivairate normal with dependency
                         % parameter rho
                         % - 'LGS': logistic with single parameter alp
                         % - 'ASL': Assymetric Logistic
    
    % Set parameters of joing dependency model                     
    switch Mdl.Jnt.Mth
        case 'MVN'  %multivariate normal N(0,Sig)
            %Sig = [1,rho; rho,1];
            Mdl.Jnt.Rho=0.6;  %correlation on [0,1];
        case 'LGS' %logistic with single parameter alp on [0,1]
            Mdl.Jnt.Alp=[0.9];
        case 'ASL' % Assymetric Logistic (misture of logistic and random)
            Mdl.Jnt.Alp=0.1;  %logistic single parameter alp on [0,1]
            Mdl.Jnt.Theta=[0,0.5]; %mixture probabilities for each margin both on [0,1]
                                 %...(proportion of 'random' data off the logistic relationship)      
    end
end



%%  Simulate
Dat=SimulateData(Mdl);
%% Save 
if ~exist('Output','dir')
   mkdir('Output') 
end
save('Output\Data','Dat','Mdl');


%%
XI=lhsdesign(100000,2).*max(Dat.Y);

F = ksdensity(Dat.Y,XI);

clf;
subplot(2,1,1)
plot(Dat.Y(:,1),Dat.Y(:,2),'k.')
hold on

subplot(2,1,2)
scatter(XI(:,1),XI(:,2),10,F,'filled')
colormap('jet')

