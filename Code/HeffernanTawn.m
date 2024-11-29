% Copyright © [2023] Shell Global Solutions International B.V. All Rights Reserved.
% SPDX-License-Identifier: Apache-2.0

classdef HeffernanTawn
    
    %functions for fitting the conditional extremes model of Heffernan and Tawn (2004) 
    % "A conditional approach to modelling multivariate extreme values"
    
    properties
        % parameters 
        Alp   %nBin x nDmn x nBoot or 1 x nBoot for slope parameter of Hefferenan
        Bet   %nBin x nDmn x nBoot
        Mu    %nBin x nDmn x nBoot
        Sig   %nBin x nDmn x nBoot  %Currently the standard deviation
        
        Rsd   % nBoot x 1 (cell array), residuals from each bootstrap
        Thr   % nBoot x 1, H&T conditional threshold
        NEP   % 1 x 1, non exceedence probability
        nBoot % 1 x 1, number of bootstraps (from marginal)
        X     % n x nBoot, Conditioned variable on Standard Margins
        Y     % n x nBoot, conditioning variable on Standard Margins
        A     % n x nBoot    %bin allocation
        RVSml    % Conditional Return Value Y | X
        RVMth = 1 %Conditional Return Value method
        % 1) Random draw from T-year max distribution (default)  P(Y|X_TYrMax)
        % 2) fixed value (by bin) to condition on to generate P(Y|X=x)
        RVFix  %nBin x nRtrPrd Fixed return value definition (only used in RVMth=2).
        RVFixMth %1 x 1 string with method used to get composite extremes.
        n     % 1 x 1, number of observations
        nDmn  % 1 x 1, number of dimensions
        SmpLclRsdOn=true;   % 1 x 1, flag for sampling residuals locally from cvr bin (1) or globally (0)
        kNearest=15;% resampling in the k nearest neighbours
        ResampleDistScale=sqrt(500); % resample distance scaling factor (written in terms of precision)    
        Sml   %Importance sampled simulation
        Alg= 'fit_newtonraphson';  %algorithm used to fit obj model either fit_fminsearch or fit_newtonraphson
        MinBet=-0.5; % minimum value that beta can possibly take        
        PrmInd   %index of parameters, alpha=1, beta=2, mu=3, sig=4;
        % delta parameter of the generalised Gaussian if Delta=1(2) corresponds to the Laplace(Gaussian) distribution
        Delta=2;
        AlphaLimit = [-1.2,1.2]; %limits for the range of alpha in the Heffernan and Tawn model 
    end
    
    
    properties(Dependent)
        nAsc;
    end
    
    properties(Hidden=true)
        FigureFolder='Figures';  %folder to save figures into (default Figures)      
    end
    
    properties(SetAccess=protected)
        nPrm     %total number of parameter       
        Prm_A    %n x 4 link between parameter
        nBin     %number of model nBins
        NonStat=[false, false, false, false];   %4 x 1 flag for non stationary parameters
        
        RsdInd   % nBoot x 1 (cell array), index of residuals from each bootstrap
        
        %% Cross Validation information
        
        CVMth=0; %0 Only Cross Validate smoothness for original dataset (fast);
        %1 Cross Validate smoothness for every bootstrap resample (slow),
        nCV=10; %no. cross-validation groups
        nSmth=10;   %no. smoothnesses tried in CV
        SmthLB=-4;   %lower bound (log10)  for smmothness range
        SmthUB=4;   %upper bound (log10)  for smmothness range       
        SmthSet     %Set of Candidate smoothness params in Cross Validation
        OptSmth     %Optimal smoothness
        CVLackOfFit %Lack of Fit of candidate smoothness parameters in CV
        MarginType = 'Laplace'; %Laplace or Gumbel
        TwoParamFlag = false;  %2P (alpha, beta) Hefferenan and Tawn model or 4P (alpha, beta , mu, sigma)
    end
    
    methods
        function obj=HeffernanTawn(Mrg,HTNEP,NonStationary,CV,SmpLclRsdOn)
            %obj=HeffernanTawn(Mrg,HTNEP,NonStationary,CV,SmpLclRsdOn)
            %INPUT
            % - Mrg 2 x 1, marignal model structure *output from stage 3
            % - HTNEP (scalar), Non-exceedance probability for conditional
            % - NonStationary flag for non stationary alpha
            %(Optional)
            %CV cross validation structure with control parameter
            %SmpLclRsdOn flag for whether we locally sample the residuals 
            
            %% Input checks
            if nargin == 0
                return
            end
            
            if ~isa(Mrg,'MarginalModel')
                error('Mrg should be a nDmn x 1 Marginal Model')
            end
            if numel(Mrg)<=1
                error('Mrg should be a nDmn x 1 Marginal Model')
            end
            obj.nDmn = numel(Mrg);
            
            %check have same number of observations for each margin
            for i = 2:obj.nDmn % loop over margins
                if size(Mrg(1).Y,1)~=size(Mrg(i).Y,1)
                    error('Marginal 1 and %d should have the same number of observations. Mrg1 nObs = %d; Mrg%d nObs = %d',i,size(Mrg(1).Y,1),i,size(Mrg(i).Y,1));
                end
            end
            %NEP must lead to a positive threshold
            switch Mrg(1).MarginType
                case 'Gumbel' %Gumbel margin
                    validateattributes(HTNEP, {'numeric'},{'<=',1,'>=',exp(-exp(0))},'HeffernanTawn','HTNEP',2);  %0<=Tau<=1
                case 'Laplace'%Laplace margin
                    validateattributes(HTNEP, {'numeric'},{'<=',1,'>=',0.5},'HeffernanTawn','HTNEP',2);  %0<=Tau<=1
            end
            
            validateattributes(NonStationary,{'logical','numeric'},{'numel',4},'HeffernanTawn','NonStationary',3)
            obj.NonStat=logical(NonStationary(:));  %flag for nonstationary obj model
            if nargin>=5
                validateattributes(SmpLclRsdOn,{'logical','numeric'},{'integer','scalar'},'HeffernanTawn','SmpLclRsdOn',5)
                obj.SmpLclRsdOn=logical(SmpLclRsdOn); %flag for sampling residuals locally from bin or globally
            end
            
            for i = 2:obj.nDmn % loop over margins
                %check have same number of bootstrap resamples for each margin
                if Mrg(1).nBoot~=Mrg(i).nBoot
                    error('Marginal 1 and %d should have the same number of bootstrap resamples. Mrg1 nBoot = %d; Mrg%d nBoot = %d',i,Mrg(1).nBoot,i,Mrg(i).nBoot);
                end
                %Check fitted each margin have the same bootstrap resamples
                if any(Mrg(1).BSInd(:)~=Mrg(i).BSInd(:))
                    error('Marginal 1 and %d should have the same bootstrap resamples',i);
                end
            end
            
            %% Smoothness parameters
            if nargin>=4  && ~isempty(CV) %% Cross Validation parameters (used for generalied Pareto fitting)
                if isfield(CV,'CVMth')
                    validateattributes(CV.CVMth, {'numeric'},{'binary'},'MarginalModel','CV.CVMth',7);
                    obj.CVMth=CV.CVMth;
                end
                if isfield(CV,'nCV')
                    validateattributes(CV.nCV, {'numeric'},{'scalar','positive'},'MarginalModel','CV.nCV',7);
                    obj.nCV=CV.nCV;      %number cross-validation groups
                end
                if isfield(CV,'nSmth')
                    validateattributes(CV.nSmth, {'numeric'},{'scalar','positive'},'MarginalModel','CV.nSmth',7);
                    obj.nSmth=CV.nSmth;    %number smoothnesses tried in CV
                end
                if isfield(CV,'SmthLB')
                    validateattributes(CV.SmthLB, {'numeric'},{'scalar'},'MarginalModel','CV.SmthLB',7);
                    obj.SmthLB=CV.SmthLB;   %lower bound (log10) for smoothness range
                end
                if isfield(CV,'SmthUB')
                    validateattributes(CV.SmthUB, {'numeric'},{'scalar'},'MarginalModel','CV.SmthUB',7);
                    obj.SmthUB=CV.SmthUB;   %upper bound (log10)  for smoothness range
                end
            end
            
            obj.FigureFolder=Mrg(1).FigureFolder;
            
            %% Pre-allocation
            obj.n=size(Mrg(1).BSInd,1);
            obj.nBoot=Mrg(1).nBoot;
            
            obj.NEP = rand(obj.nBoot,1).*(HTNEP(2)-HTNEP(1))+HTNEP(1);
            obj.nBin=Mrg(1).Bn.nBin;
            
            Prm={'Alp','Bet','Mu','Sig'};
            
            if obj.TwoParamFlag
                obj.NonStat(3:4)=false;
            end
            
            obj.PrmInd=[];
            
            if obj.nBin == 1 %all stationary since only 1 bin
                obj.NonStat=false(4,1);
                obj.nPrm=ones(4,1);
                obj.PrmInd=(1:4)';
            else
                for iP=1:4 %loop over parameters
                    if obj.NonStat(iP) %get number of model bins (number of alpha parameters).
                        obj.(Prm{iP})=NaN(obj.nBin,obj.nDmn-1,obj.nBoot);
                        obj.nPrm(iP)=obj.nBin;
                        obj.PrmInd=[obj.PrmInd;ones(obj.nBin,1).*iP];
                    else  %stationary model case
                        obj.(Prm{iP})=NaN(1,obj.nDmn-1,obj.nBoot);
                        obj.nPrm(iP)=1;
                        obj.PrmInd=[obj.PrmInd;iP];
                    end
                end
            end
            obj.Thr=NaN(obj.nBoot,obj.nDmn-1);
            obj.Rsd=cell(obj.nBoot,1); %different numbers of residuals for each bootstrap so store in cell
            obj.RsdInd=cell(obj.nBoot,1); %bin index of exceedences used in local sampling of residual;
            
            obj.Y=NaN(obj.n,obj.nDmn-1,obj.nBoot);
            obj.A=NaN(obj.n,obj.nBoot);
            obj.Prm_A=ones(obj.n,obj.nBoot,4);
            obj.OptSmth=NaN(obj.nBoot,1);
            
            if any(obj.NonStat)
                obj.SmthSet=logspace(obj.SmthLB,obj.SmthUB,obj.nSmth); %try range smoothness penalties for sigma varying by bin
            else
                obj.nSmth=1;
                obj.SmthSet=0;    %Switch off smoothness
            end
            obj.CVLackOfFit=NaN(obj.nSmth,obj.nBoot);
            
        end %HeffernanTawn constructor
        
        function obj=Fit(obj,Mrg)
            %Fit Heffernan and Tawn model and compute return values
            %
            %INPUT
            %obj object
            %Mrg Marginal model object
            %
            %OUPTUT
            %obj fitted obj model with return values
            
            obj.X=Margins(Mrg(1));
            
            %% Fit H&T Model
            for iBt=1:obj.nBoot           %loop over bootstrap resamples
                fprintf('Fitting for bootstrap sample %d of %d\n',iBt,obj.nBoot);
                %transform conditioned variable to Standard margins for iBt'th bootstrap
                obj.Thr(iBt)=Mrg(1).INV_Standard(obj.NEP(iBt));
                IExc= obj.X(:,iBt)>obj.Thr(iBt);   %threshold exceedences
                % dealing with excluding with pairs that are not
                % transformed as they do not exist
                %transform conditioning variable to Standard margins for iBt'th bootstrap
                
                INotNaN=Mrg(1).BSInd(:,iBt)>0;
                J=Mrg(1).BSInd(INotNaN,iBt);
                obj.A(INotNaN,iBt)=Mrg(1).Bn.A(J);
                for iP=1:4
                    if obj.NonStat(iP)
                        obj.Prm_A(INotNaN,iBt,iP)=Mrg(1).Bn.A(J);
                    end
                end             
                
                obj.RsdInd{iBt}=obj.A(IExc,iBt);
                for iDmn=2:obj.nDmn
                    obj.Y(:,iDmn-1,iBt)=Margins(Mrg(iDmn),iBt);
                end
                %% Fit Model
                obj=Fit_private(obj,IExc,iBt);
            end
            
            %% Return Value Simulation
            obj=RVSimulations(obj,Mrg);
        end %Fit
        
        function Sml=SimulateIS(obj,Mrg,nRls)
            % simulate importance sampling in X MC in Y|X
            % Simulation modifed to force same reponses for all bins nRls
            %INPUTS:
            % - nDmn x 1 Mrg, MarginalModel class Mrg (used to transform back to orginal scale)
            % - 1 x 1 nRls, number of realisations under the model
            %OUTPUT:
            % - Data structure Sml, data simulated from H&T model on Original,
            % Uniform and Standard margins
            
            Sml.A=reshape((1:obj.nBin).*ones(nRls,1),[],1);  %bin allocation
            Sml.I=reshape(randi(obj.nBoot,nRls,1).*ones(1,obj.nBin),[],1);%decide which bootstrap sample to use over all bootstraps obj.nBoot
            %% Simulate covariate
            RatNrm=Mrg(1).Rat./sum(Mrg(1).Rat,1); %nBin x nBoot
            %Simulate covariate with the right rate
            Sml.X=SampleCovariateFromBin(Mrg(1).Bn,Sml.A);  %Sample covariate value uniformly from  within bin
            
            if obj.nBin==1
                J=Sml.I;
                Sml.W=RatNrm(J)'; %nRls x 1 %probability weight for chosen bin;
            else
                J=sub2ind([obj.nBin,obj.nBoot],Sml.A,Sml.I);  %common index across bin allocation and bootstrap
                Sml.W=RatNrm(J); %nRls x 1 %probability weight for chosen bin;
            end
            %% Generate X uniformly on standard margins
            Sml.StnMrg= NaN(nRls.*obj.nBin,obj.nDmn); %simulation resluts on standard margins
            
            BndStnMrg=[Mrg(1).INV_Standard(1e-7),Mrg(1).INV_Standard(1-1e-7)];
            G=rand(nRls.*obj.nBin,1).*(BndStnMrg(2)-BndStnMrg(1))+BndStnMrg(1);
            Sml.logfog=Mrg(1).LogPDF_Standard(G)+log(BndStnMrg(2)-BndStnMrg(1));
            Sml.logfog(:,2:obj.nDmn)=0;
            Sml.logfog_Stn=Sml.logfog;
            %% Simulate on standard margins        
            Sml.StnMrg=obj.SampleYgX_Stn(G,Sml.I,Sml.A);
            
            %% Transform to uniform
            Sml.Unf=CDF_Standard(Mrg(1),Sml.StnMrg);
            
            %% Transform back to original margin
            Sml.Org=NaN(size(Sml.Unf));
            for iDmn=1:obj.nDmn %loop over dimension
                Sml.Org(:,iDmn)=Mrg(iDmn).INV(Sml.Unf(:,iDmn),Sml.I,Sml.A);
            end
            
        end %SimulateIS
        
        function Sml=SimulateMC(obj,Mrg,smlMethod,nRls,A,nA,RVX)
            %Compute conditional quantiles of Y|X for given X values
            % function to simulate Monte carlo samples from the Heffernan
            % and Tawn model - with the option to specify X in 3 different
            % ways
            %
            %Inputs
            % obj Heffernan and Tawn model - fit upfront
            % Mrg marginal model object, which includes the precalculated values of X
            % smlMethod - method choice:
            % randX - X is generated at random from the marginal model
            % retValX - X is generated from the return value model
            % userX - user provided X, which is provided as a user input
            % in that case RVX [nA X nBt x nRtr]
            % A vector of bins to compute; for omni should be A
            % nA number of bins (default max(A));
            % in the omnidirectional case we would have nA as
            % (optional only used in userX case)
            % RVX - [nA x nBt x nRtr] return value calculations of X
            %
            % Outputs
            % Sml = Populated structure of Y|X
            % Orig: [nA x nRls x nDmn x nRtr] cond quantiles of Y given X on original margins
            % StnMrg: [nA x nRls x nDmn x nRtr] cond quantiles of Y given X on std margins
            % Unif: [nA x nRls x nDmn x nRtr] cond quantiles of Y given X on uniform margins
            % A: [nA x nRls] vector of bins versus the number of realisations
            % I: [nRls x 1] bootstrap index
            
            %% Input checks
            % validate
            if nargin<=2
                smlMethod="randX";
            else
                smlMethod=validatestring(smlMethod,["randX","retValX","userX"]);
            end
            
            % optional arguments for nRls, A and nA
            if nargin<=3
                nRls=1000; %default number of MC draws
            end
            Sml.nRls=nRls;
            
            %Input A and nA super bins definitions
            if nargin<=4
                A=ones(obj.nBin,1); %omni situation
                %A=(1:obj.nBin)';% per original bin
            end
            
            %Input 5 number of super bins
            if nargin<=5
                nA=max(A);
            end
            if nA<max(A)
                error('Inconsistency in number of bins definition')
            end
            
            %Calculation of the omnidirectional version
            if nA>1
                OmniOn=true;
            else
                OmniOn=false;
            end
            
            % populating conditional values of X
            if nargin<=6 && strcmp(smlMethod,"userX")
                error('User defined return values need to be provided in the case where smlMethod is UserX.')
            end
            
            if strcmp(smlMethod,"userX")
                OmniOn=false; %this needs to be handled outside in this case
            end
            
            %% switching between the different simulation methods
            switch smlMethod
                case "randX" %sample uniform values of X
                    nRtr=1;
                    Sml=Mrg(1).sample_MC(nRls,A,nA);
                case "retValX" %X defined in terms of the return values
                    nRtr=Mrg(1).nRtr;
                    Sml=Mrg(1).RVSml;
                    if isempty(Sml) %no marginal RV
                        return
                    end
                case "userX" %user defined X
                    nRtr=size(RVX,3);
                    Sml=Mrg(1).sample_MC_CondX(nRls,A,nA,RVX);
            end          
            
            % arrays for original, uniform and standard margins
            % empty array for the sector identification
            tOrg=NaN(nA*nRls*nRtr,obj.nDmn);
            tUnf=NaN(nA*nRls*nRtr,obj.nDmn);
            %unwrap marginal output
            tA=Sml.A(:);
            tI=Sml.I(:);
            % find nans in A
            I_nanA=~isnan(tA);
            tUnf(:,1)=Sml.Unf(:);
            tOrg(:,1)=Sml.Org(:);
            
            %% sample covariates
            Sml.X=NaN(nA*nRls*nRtr,Mrg(1).nCvr);
            Sml.X(I_nanA,:)=SampleCovariateFromBin(Mrg(1).Bn,tA(I_nanA));
            tX=Sml.X; %[(nA*nRls*nRtr) x nCvr] copied to allow omni case
            %% simulate on standard scale
            % [nA nRls 1 nRtr] into something that [(nA*nRls*1*nRtr) x 1]
            tStnMrg=obj.SampleYgX_Stn(INV_Standard(Mrg(1),tUnf(:,1)),tI,tA); % [(nA*nRls*nRtr) x nAsc]
            %% Transform to uniform
            tUnf=Mrg(1).CDF_Standard(tStnMrg);
            %% Transform back to original margin
            for iDmn=1:obj.nDmn %loop over dimension - populating only the associated margins
                tOrg(I_nanA,iDmn)=Mrg(iDmn).INV(tUnf(I_nanA,iDmn),tI(I_nanA),tA((I_nanA)));
            end
            %% Identification of maximum value and corresponding margin
            if OmniOn
                [~,J]=max(reshape(tOrg(:,1),[nA,nRls*nRtr]),[],1);
                Sml.J=sub2ind([nA,nRls*nRtr],J,1:(nRls*nRtr));
            end
            
            % reshape data back into [nA x nRls x nDmn x nRtr]
            Sml.StnMrg=permute(reshape(tStnMrg,[nA,nRls,nRtr,obj.nDmn]),[1,2,4,3]); 
            Sml.Org=permute(reshape(tOrg,[nA,nRls,nRtr,obj.nDmn]),[1,2,4,3]); 
            Sml.Unf=permute(reshape(tUnf,[nA,nRls,nRtr,obj.nDmn]),[1,2,4,3]); 
            Sml.X=permute(reshape(tX,[nA,nRls,nRtr,Mrg(1).Bn.nCvr]),[1,2,4,3]);
            %%Omnidirectional calculation addition
            % nA in this case is nA+1 due to the addition of the omni
            if OmniOn
                Sml.StnMrg=cat(1,Sml.StnMrg,permute(reshape(tStnMrg(Sml.J,:),[1,nRls,nRtr,obj.nDmn]),[1,2,4,3])); % [nA x NRls x nAsc x nRtr]
                Sml.Org=cat(1,Sml.Org,permute(reshape(tOrg(Sml.J,:),[1,nRls,nRtr,obj.nDmn]),[1,2,4,3])); % [nA x NRls x nAsc x nRtr]
                Sml.Unf=cat(1,Sml.Unf,permute(reshape(tUnf(Sml.J,:),[1,nRls,nRtr,obj.nDmn]),[1,2,4,3])); % [nA x NRls x nAsc x nRtr]
                Sml.X=cat(1,Sml.X,permute(reshape(tX(Sml.J,:),[1,nRls,nRtr,Mrg(1).Bn.nCvr]),[1,2,4,3])); % [nA x NRls x nCvr x nRtr]
            end
            % Populated Sml object
        end %SimulateMC
        
        function Plot(obj,Mrg)         
            %% Simulate under model
            nRls=10000;
            MCSml=SimulateMC(obj,Mrg,"randX",nRls);
            
            %% Plot data and simulation
            figure;
            clf;
            %summary of output for a range of quantiles
            [sOrig,sStnMrg]=obj.GetSummary('quantile',[0.05,0.37,0.5,0.975]);
            
            for iDmn = 2:obj.nDmn
                
                subplot(2,obj.nDmn-1,iDmn-1)
                plot(MCSml.StnMrg(1,:,1),MCSml.StnMrg(1,:,iDmn),'.','color',[1,1,1]*0.8)
                hold on
                plot(obj.X(:,1),obj.Y(:,iDmn-1,1),'k.','markersize',10)
                hold on
                plot(squeeze(sStnMrg(end,1,:,2)),squeeze(sStnMrg(end,iDmn,:,[1,3,4])),'r--','linewidth',2);
                xlabel(sprintf('%s: Conditioning variable',Mrg(1).RspLbl))
                ylabel(sprintf('%s: Conditioned variable',Mrg(iDmn).RspLbl));
                title(sprintf('%s|%s: Conditional extremes simulation on %s margins',Mrg(iDmn).RspLbl,Mrg(1).RspLbl,Mrg(1).MarginType))
                legend('Simulation','Data','location','NorthWest')
                axis tight
                box on
                grid on

                subplot(2,obj.nDmn-1,obj.nDmn-1+iDmn-1)
                plot(MCSml.Org(1,:,1),MCSml.Org(1,:,iDmn),'.','color',[1,1,1]*0.8)
                hold on
                plot(Mrg(1).Y,Mrg(iDmn).Y,'k.','markersize',10)
                hold on
                plot(squeeze(sOrig(end,1,:,2)),squeeze(sOrig(end,iDmn,:,[1,3,4])),'r--','linewidth',2);
                xlabel(sprintf('%s: Conditioning variable',Mrg(1).RspLbl))
                ylabel(sprintf('%s: Conditioned variable',Mrg(iDmn).RspLbl));
                title(sprintf('%s|%s: Conditional extremes simulation on original scale',Mrg(iDmn).RspLbl,Mrg(1).RspLbl))
                axis tight
                grid on
            end
            savePics(fullfile(obj.FigureFolder,'Stg4_HT_1_SmlvsData'))
            
            %% Residual Diagnostics orignal sample
            figure;
            clf;
            tRsd=cell2mat(obj.Rsd);
            
            %histogram of residuals in ith variable
            for iDmn = 1:(obj.nDmn-1)
                subplot((obj.nDmn-1),2,2*iDmn-1)
                if verLessThan('Matlab','8.5')
                    hist(tRsd(:,iDmn)); %#ok
                    h = findobj(gca,'Type','patch');
                    set(h,'FaceColor',[1 1 1]*0.5);
                    set(h,'linestyle','none')
                else
                    histogram(tRsd(:,iDmn),'edgecolor','none','facecolor','k','normalization','pdf')
                end
                axis tight
                title(sprintf('%s|%s: Empirical density of residuals',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl))
                xlabel('Residual')
                ylabel('Empirical density');
                grid on;
            end
            %qq plot of residuals in ith variable
            for iDmn = 1:(obj.nDmn-1)
                
                subplot(obj.nDmn-1,2,2*iDmn);
                h=qqplot(tRsd(:,iDmn));
                if verLessThan('Matlab','8.5')
                    h = findobj(gca,'Type','patch');
                    set(h,'Marker','.');
                    set(h,'MarkerEdgeColor',[0 0 0]);
                    set(h,'linestyle','--')
                    set(h,'LineWidth',1)
                else
                    h(1).Marker='.';
                    h(1).MarkerEdgeColor=[0,0,0];
                end
                ylabel('Empirical quantiles of residuals')
                xlabel('Standard normal quantiles')
                axis tight
                box on
                title(sprintf('%s|%s: Residuals QQ plot',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl))
                grid on;
            end
            savePics(fullfile(obj.FigureFolder,'Stg4_HT_2_ResidualDiagnostics1'))
            
            
            %% Residual diagnostics: residuals against covariates
            figure;clf;
            iP = 0; %counter on subplot number
            for iDmn = 1:(obj.nDmn-1)
                for iC=1:Mrg(1).nCvr
                    iP = iP +1;
                    subplot(obj.nDmn-1,Mrg(1).nCvr,iP)
                    IExc=obj.X(:,1)>obj.Thr(1);   %threshold exceedences
                    IExc=IExc(Mrg(1).BSInd(:,1)>0);
                    plot(Mrg(1).X(IExc,iC),obj.Rsd{1}(:,iDmn),'k.')
                    axis tight
                    box on
                    grid on
                    xlabel(Mrg(1).CvrLbl(iC))
                    title(sprintf('%s|%s: Residuals on %s',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl,Mrg(1).CvrLbl{iC}))
                    grid on;
                end
            end
            savePics(fullfile(obj.FigureFolder,'Stg4_HT_2_ResidualDiagnostics2'))
            
            %% Parameter Plot (bootstrap uncertainty)
            LgnLbl=cell(obj.nDmn-1,1);
            for iDmn = 1:(obj.nDmn-1)
                LgnLbl{iDmn} = sprintf('%s|%s',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl);
            end
            Lbl={'\alpha','\beta','\mu','\sigma'};
            
            for iDmn = 1:(obj.nDmn-1) %different parameters figure for each dimension
                figure;
                clf;
                
                PrmLbl={'Alp','Bet','Mu','Sig'};
                %loop over parameters
                for iP=1:4          
                    if obj.NonStat(iP)
                        for iC=1:Mrg(1).nCvr
                            if Mrg(1).nCvr>2
                                subplot(4,Mrg(1).nCvr,(iP-1)*Mrg(1).nCvr+iC)
                            else
                                subplot(2*Mrg(1).nCvr,2,2*(iC-1)+1);
                            end
                            hold on
                            PlotParameter(Mrg(1).Bn,squeeze(obj.(PrmLbl{iP})(:,iDmn,:)),iC,'color',[0,0,0],'linewidth',2)
                            xlabel(Mrg(1).CvrLbl(iC))
                            ylabel(Lbl{iP})
                            box on
                            grid on
                        end
                    else
                        subplot(2,2,iP)
                        Prm=squeeze(obj.(PrmLbl{iP})(1,iDmn,:));
                        
                        histogram(Prm,'edgecolor','none','facecolor','k','normalization','pdf')
                        
                        hold on
                        xlabel(Lbl{iP})
                        ylabel('Empirical density');
                        grid on;
                    end
                end
                axes('position',[0.1300    0.1100    0.7750    0.8150]);
                axis off
                title(sprintf('%s|%s: Conditional extremes model parameters',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl))
                
                savePics(fullfile(obj.FigureFolder,sprintf('Stg4_HT_3_Parameters_%s',Mrg(iDmn+1).RspLbl)))
            end
            
            %% Threshold plot (alpha as function of NEP)
            %Alpha varies by covariate bin, so to avoid over-crowded plots; have a different version of this plot for each dimension.
            nPlt1=ceil(sqrt(obj.nPrm(1))); %max size nPlt x nPlt
            nPlt2=ceil(obj.nPrm(1)./nPlt1);
            
            for iDmn = 1:(obj.nDmn-1)
                figure;
                clf;
                for iAlp = 1:obj.nPrm(1)
                    subplot(nPlt2,nPlt1,iAlp)
                    tAlp=squeeze(obj.Alp(iAlp,iDmn,:));
                    plot(obj.NEP,tAlp,'k.','markersize',20)
                    
                    if Mrg(1).nBoot>=20
                        hold on
                        nNEP=10; %number of bins for threshold diagnostic plot
                        NEPEdg=linspace(min(obj.NEP),max(obj.NEP),nNEP+1);
                        NEPBin=(NEPEdg(1:end-1)+NEPEdg(2:end))/2;
                        [~,tA] = histc(obj.NEP,NEPEdg);
                        tA(tA>nNEP)=nNEP;
                        
                        
                        S=accumarray(tA,tAlp,[nNEP,1],@median,NaN);
                        Slb=accumarray(tA,tAlp,[nNEP,1],@(x)quantile(x,0.025),NaN);
                        Sub=accumarray(tA,tAlp,[nNEP,1],@(x)quantile(x,0.975),NaN);
                        
                        NanI=isnan(S);
                        plot(NEPBin(~NanI),S(~NanI),'r-','linewidth',2)
                        plot(NEPBin(~NanI),Slb(~NanI),'r--','linewidth',2)
                        plot(NEPBin(~NanI),Sub(~NanI),'r--','linewidth',2)
                    end
                    grid on
                    
                    if iAlp == 1
                        xlabel('obj NEP')
                        ylabel('\alpha')
                        if obj.nPrm(1)==1 %stationary case
                            title(sprintf('%s|%s: Conditional extremes parameter stability by threshold',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl))
                        else
                            title(sprintf('%s|%s: Conditional extremes parameter stability by threshold \n Bin %s',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl,Mrg(1).Bn.BinLbl{iAlp}))
                        end
                    else
                        title(sprintf('Bin %s',Mrg(1).Bn.BinLbl{iAlp}))
                    end
                end
                hold off;
                savePics(fullfile(obj.FigureFolder,sprintf('Stg4_HT_4_AlphaThresholdStability_%s',Mrg(iDmn+1).RspLbl)))
            end
            
            
            %% Lack of fit plot
            if any(obj.NonStat)
                figure;
                clf;
                if obj.nBoot>1 && obj.CVMth==1
                    plot(log10(obj.SmthSet),nanmedian(obj.CVLackOfFit,2),'k-','linewidth',2)
                    hold on
                    plot(log10(obj.SmthSet),quantile(obj.CVLackOfFit,0.025,2),'k--','linewidth',2)
                    plot(log10(obj.SmthSet),quantile(obj.CVLackOfFit,0.975,2),'k--','linewidth',2)
                else
                    plot(log10(obj.SmthSet),obj.CVLackOfFit(:,1),'k-','linewidth',2)
                end
                axis tight
                hold on
                grid on
                plot(log10(median(obj.OptSmth))*[1,1],ylim,'r--','linewidth',2)
                ylabel('Lack of fit')
                xlabel('$\log_{10}(\tilde{\lambda)}$','Interpreter','latex')
                title(sprintf('Conditional extremes cross-validation lack of fit'))
                grid on;
                savePics(fullfile(obj.FigureFolder,'Stg4_HT_5_SectorGoodnessOfFit'))
            end
            
            %% Cond. RV CDFs
            ColMat=hsv(obj.nBin);
            
            figure;
            clf;
            c=0;
            for iDmn=2:obj.nDmn %loop over associated variables
                for iRtr=1:Mrg(1).nRtr
                    c=c+1;
                    
                    subplot(obj.nDmn-1,Mrg(1).nRtr,c)
                    grid on
                    % Per bin: -----------
                    if obj.nBin<16  && obj.nBin > 1  %plot directional CDFs when there are no more than 15 bins
                        
                        tX=sort(obj.RVSml.Org(1:obj.nBin,:,iDmn,iRtr),2,'MissingPlacement','first');
                        
                        for iBin = 1:obj.nBin
                            plot(tX(iBin,:),linspace(0,1,obj.RVSml.nRls),'linewidth',2,'color',ColMat(iBin,:))
                            hold on
                        end
                    end
                    
                    %Omni: --------------
                    tX=sort(obj.RVSml.Org(end,:,iDmn,iRtr),2,'MissingPlacement','first');
                    plot(tX,linspace(0,1,obj.RVSml.nRls),'k-','linewidth',2)
                    hold on
                    talPrb=1e-3;
                    lb=min(quantile(obj.RVSml.Org(:,:,iDmn,iRtr),talPrb,2));
                    ub=max(quantile(obj.RVSml.Org(:,:,iDmn,iRtr),1-talPrb,2));
                    xlim([lb ub]);
                    
                    plot(xlim,[0.5,0.5],'--k')
                    plot(xlim,[exp(-1),exp(-1)],'--k')
                    
                    ylabel('Cumulative probability')
                    xlabel(sprintf('%s | %s',Mrg(iDmn).RspLbl, Mrg(1).RspLbl))
                    grid on
                    
                    if (obj.nBin<16  && obj.nBin >1 ) && (c == 1)
                        legend([Mrg(1).Bn.BinLbl(:);'Omni'],'location','best','fontsize',6);
                    end
                    
                    if iDmn==2
                        title(sprintf('Conditional maximum over %g years',Mrg(iDmn).RtrPrd(iRtr)))
                    end
                end
            end
            savePics(fullfile(obj.FigureFolder,'Stg4_HT_6_ConditionalReturnValueCDF'))
        end %plot
        
        function [sOrig,sStnMrg]=GetSummary(obj,summarystat,quantilelvl)
            % function to summarise monte Carlo simulation from a Heffernan and
            % Tawn model
            % summarystat default='median' other options 'mean', 'quantile'
            % quantilelvl default = 0.5
            %
            %
            % OUTPUT
            % sOrig  nA x nDmn x nRtr  (x nQnt)
            % sStnMrg  nA x nDmn x nRtr  (x nQnt)
            
            if isempty(obj.RVSml)
                obj.RVSml=obj.Sml;
            end
            
            % default to the median calculation
            if nargin<=2
                summarystat='median';
            end
            
            if nargin>=3
                qLvl=quantilelvl;
            else
                qLvl=0.5;
            end
                        
            %original margin
            switch summarystat
                case 'mean'
                    sOrig=permute(nanmean(obj.RVSml.Org,2),[1,3,4,2]);
                case {'median','quantile'}
                    obj.RVSml.Org(isnan(obj.RVSml.Org))=-Inf;
                    sOrig=permute(quantile(obj.RVSml.Org,qLvl,2),[1,3,4,2]);
                    sOrig(isinf(sOrig))=NaN;
            end
            
            %if also want standard margin out
            if nargout>1
                switch summarystat
                    case 'mean'
                        sStnMrg=permute(nanmean(obj.RVSml.StnMrg,2),[1,3,4,2]);
                    case {'median','quantile'}
                        obj.RVSml.StnMrg(isnan(obj.RVSml.StnMrg))=-Inf;
                        sStnMrg=permute(quantile(obj.RVSml.StnMrg,qLvl,2),[1,3,4,2]);
                        sStnMrg(isinf(sStnMrg))=NaN;
                end
            end
        end % GetSummary
        
        function nAsc = get.nAsc(obj) %get no associated variables
            nAsc = obj.nDmn-1;
        end %nAsc
        
    end %methods
    methods(Access = protected)
              
        function StnMrg=SampleYgX_Stn(obj,X,I,A)
            %Input
            %X  nRls x 1 conditioning variable on standard margins
            %I  nRls x 1 bootstrap index vector
            %A  nRls x 1 bin index
            
            %Output
            %StnMrg nRls x nDmn   
            
            
            APrm=A.*(obj.NonStat');
            APrm(APrm==0)=1; %nRls x 4
            
            nRls=numel(X);            
            StnMrg=NaN(nRls,obj.nDmn);
            StnMrg(:,1)=X; %use random gumbel/laplace
            
            for iBt=1:obj.nBoot
                % indicate exceedances 
                IBlw=X<=obj.Thr(iBt) & I==iBt;%threshold on the standard scale
                IAbv=X>obj.Thr(iBt) & I==iBt;%threshold on the standard scale
                % threshold exceedances      
                if any(IAbv)
                    Z=obj.BinResampling(obj.Rsd{iBt},obj.RsdInd{iBt},A(IAbv),[],obj.kNearest);
                    %nRsd x nDmn-1 joint residuals in each dimension
                    idNaN=isnan(Z(:,1));
                    % if we don't have any samples to take in that bin
                    if any(idNaN)
                        iSamp=randi(size(obj.Rsd{iBt},1),sum(idNaN),1);
                        Z(idNaN,:)=obj.Rsd{iBt}(iSamp,:);
                    end                   
                    %get conditioning value from H&T model
                    StnMrg(IAbv,2:end)=obj.Alp(APrm(IAbv,1),:,iBt).* X(IAbv) +...
                        X(IAbv).^obj.Bet(APrm(IAbv,2),:,iBt).* ( obj.Mu(APrm(IAbv,3),:,iBt)+ obj.Sig(APrm(IAbv,4),:,iBt).*Z);
                end
                    % those below the exceedances
                if any(IBlw)
                    IBlwObjX=obj.X(:,iBt)<=obj.Thr(iBt);
                    Z=[obj.X(IBlwObjX,iBt),obj.Y(IBlwObjX,:,iBt)];
                    if ~isempty(Z)
                        ZSamp=obj.BinResampling(Z,obj.A(IBlwObjX,iBt),A(IBlw),X(IBlw),obj.kNearest);
                        ZSamp(isnan(ZSamp))=-Inf;
                        StnMrg(IBlw,2:end)=ZSamp(:,2:end);
                    end
                end
                StnMrg(isnan(StnMrg))=-Inf;
                % TODO check X for NaNS and treat it differently     
            end %iBoot
        end %Sample
        
        function ZSamp=BinResampling(obj,BinPopZ,BinPopA,A,Z,~)
            % obj=BinResampling(obj,Z,A,kNearest)
            % Function to deal deal with sampling those instances where
            %Input
            %BinPopZ  nRls x nAsc indices of the population
            %BinPopA  nRls x 1 indices of the population
            %A  nSmp x 1 specific index of the bin
            %(Optional)
            % Z nSmp x 1 also match based on distance in Z (only matches
            % first column of BinPOpZ default is [] i.e. not used
            %Output
            %ZSamp nSmp x nDmn location of the chosen index
            if nargin<5
                Z=[];
            end
            nSmp = numel(A);
            
            % check if the number of values to sample from
            if size(BinPopA,1)>0
                ZSamp=NaN(nSmp,size(BinPopZ,2));
                for iSmp=1:nSmp
                    % if samples exist in the bin
                    if sum(BinPopA==A(iSmp))>0
                        tZ=BinPopZ(BinPopA==A(iSmp),:);
                        if isempty(Z)
                            ISamp=randi(size(tZ,1),1);
                        else
                            ISamp=CovResampling(obj,tZ,Z(iSmp),1);
                        end
                        ZSamp(iSmp,:)=tZ(ISamp,:);
                    else
                        ZSamp(iSmp,:)=NaN;
                    end
                end %iSmp
            else
                ZSamp=[];
            end
        end %BinResampling
                
        function ISamp=CovResampling(obj,BinPopZ,Z,nSmp)
            %index covariate resampling 
            %use squared exponential distance weighting.
            %Input
            %BinPopZ  nRls x nAsc indices of the population
            % Z nSmp x 1 also match based on distance in Z (only matches
            % nSmp number of samples 
            % first column of BinPopZ default is [] i.e. not used
            %Output
            %ISamp nSmp x nAsc indices of for resampling
            if isempty(Z)
                ISamp=randi(size(BinPopZ,1),nSmp,1);
            else
                ISamp=NaN(nSmp,1);
                nu = obj.ResampleDistScale;
                for iSmp=1:nSmp
                    Wunscaled = exp(-0.5.*nu.*abs(BinPopZ(:,1)-Z(iSmp))); %k x 1 distance to sampled Z
                    if sum(Wunscaled)==0
                        ISamp(iSmp)=randi(size(BinPopZ,1),1,1);
                    else
                        W = Wunscaled./sum(Wunscaled);
                        ISamp(iSmp)=randsample(size(BinPopZ,1),1,true,W);
                    end
                end %iSmp
            end
        end %CovResampling
      
        function obj=Fit_private(obj,IExc,iBt)
            % obj=Fit_private(obj,IExc,iBt)
            % Fit single bootstrap of obj model
            % INPUTS:
            % - n x 1 IExc logical array of exceedances
            % - 1 x 1 iBt bootstrap number
            % OUTPUTS:
            % - obj.Alp(:,:,iBt) fitted H&T parameters for iBt'th bootstrap
            % - obj.Bet(:,:,iBt) fitted H&T parameters for iBt'th bootstrap
            % - obj.Mu(:,:,iBt) fitted H&T parameters for iBt'th bootstrap
            % - obj.Sig(:,:,iBt) fitted H&T parameters for iBt'th bootstrap
            % - obj.Rsd{:,iBt} residuals for iBt'th bootstrap
            
            %% Get exceedances
            nExc=sum(IExc);
            tX=obj.X(IExc,iBt); %exceedances for iBt'th bootstrap
            tY=obj.Y(IExc,:,iBt);  %conditioning value given exceedance in conditioned value
            tA=permute(obj.Prm_A(IExc,iBt,:),[1,3,2]);  %nExc x 4
            if sum(isnan(tA))>0
                warning('NaNs in allocation')
            end
            
            %% Flag for do Cross Validation
            if (obj.CVMth==0 &&  iBt>1) || obj.nSmth==1
                %% Cross validation off
                if obj.nSmth==1  %only one smoothness
                    obj.OptSmth(iBt)=obj.SmthSet;
                else  %use first bootstrap instead of cross validating every time
                    obj.OptSmth(iBt)=obj.OptSmth(1);
                end
            else
                fprintf('Starting Cross Validation to estimate obj alpha smoothness:\n');
                LackOfFit=NaN(obj.nCV,obj.nSmth); %lack of fit
                ICV=randi(obj.nCV,nExc,1); %split nExc observations into nCV cross validation groups
                for iCV=1:obj.nCV  %cross validation loop
                    fprintf('.')
                    %Pull out fit set
                    tXFit=tX(ICV~=iCV); %main response
                    tYFit=tY(ICV~=iCV,:); %associated variable
                    tAFit=tA(ICV~=iCV,:);
                    %Pull out prediction set
                    tXPrd=tX(ICV==iCV);
                    tYPrd=tY(ICV==iCV,:);
                    tAPrd=tA(ICV==iCV,:);
                    for iLam=obj.nSmth:-1:1 %Smoothness penalties
                        PrmHat=obj.(obj.Alg)(tXFit,tYFit,tAFit,obj.SmthSet(iLam));
                        LackOfFit(iCV,iLam)=obj.likelihood(tXPrd,tYPrd,tAPrd,PrmHat);
                    end %SigPen
                end %CV
                fprintf('\n')
                obj.CVLackOfFit(:,iBt)=sum(LackOfFit,1)';
                [~,tI]=min(obj.CVLackOfFit(:,iBt));
                obj.OptSmth(iBt)=obj.SmthSet(tI)';
            end
            
            %% fit model: find optimal parameter values for heffernan and tawn model
            tPrm=obj.(obj.Alg)(tX,tY,tA,obj.OptSmth(iBt));
            [obj.Alp(:,:,iBt), obj.Bet(:,:,iBt),obj.Mu(:,:,iBt), obj.Sig(:,:,iBt)]=p2Prm(obj,tPrm);
            %% keep residuals
            tXpB=tX.^obj.Bet(tA(:,2),:,iBt);
            mn=tY-obj.Alp(tA(:,1),:,iBt).*tX-obj.Mu(tA(:,3),:,iBt).*tXpB;
            std= obj.Sig(tA(:,4),:,iBt).*tXpB;
            obj.Rsd{iBt}=mn./std;
            if any(~isreal(obj.Rsd{iBt}(:)))
                error('Complex Residuals found!!\n')
            end
            
        end %fit
        
        function PrmHat=fit_fminsearch(obj,X,Y,A,Lam)
            %Fit Heffernan and Tawn model using simple fminsearch algorithm, works well for small
            %number of bins
            %INPUT
            %X    %conditioning data
            %Y    %conditioned data
            %A    %bin allocation
            %Lam  %smoothness parameter
            
            %% Starting Solution
            p0=obj.startingSolution(X,Y,A,Lam);
            
            opts=optimset('display','off');
            
            PrmHat=fminsearch(@(p)obj.likelihood(X,Y,A,p,Lam),p0,opts);
            
            if obj.TwoParamFlag
                PrmHat=[PrmHat;zeros(obj.nDmn-1,1)';ones(obj.nDmn-1,1)'];
            end
        end
        
        function p=fit_newtonraphson(obj,X,Y,A,Lam)
            %Fit Hefferenan and Tawn model using newton raphSon algorithm, works well for
            %bigger cases
            %INPUT
            %X    %conditioning data
            %Y    %conditioned data
            %A    %bin allocation
            %Lam  %smoothness parameter
            %OUTPUT
            %p    %parameters of HT expect [nA x nDmn-1]
            b0=0;
            % determine and suitable starting solutiom
            % Hessian checks are also perforemd
            p=obj.startingSolution(X,Y,A,b0,Lam);
            
            FitInd=[1,2,1,2];
            if obj.TwoParamFlag
                FitInd=FitInd(obj.PrmInd(obj.PrmInd<=2));
            else
                FitInd=FitInd(obj.PrmInd);
            end
            
            maxIter=1000; %maximum number of iterartions of solver
            tol=1e-2; %tolerance for convergence
            Dlt0=[0.2,0.01];%[0.2,0.2];%step size on [0,1])
            Dlt=Dlt0; %step size on [0,1])
            
            for iIt=1:maxIter %iteration loop
                Chg=0;
                % iterate between pairs (alpha,mu) & (beta,sigma)
                for iP=1:2
                    %G p x nDmn, H p x p nDmn
                    [G,H]=Gradient(obj,iP,X,Y,A,p,Lam);
                    minEig=eigenCheck(obj,X,Y,A,p,Lam);
                    if minEig<0 %trying to deal with a non positive semi-definite covariance matrix
                        continue
                    end
                    newp=p;
                    for iDmn=1:obj.nDmn-1
                        newp(FitInd==iP,iDmn)= p(FitInd==iP,iDmn)-Dlt(iP).*(H(:,:,iDmn)\G(:,iDmn));
                    end
                    
                    if obj.checkParam(newp)
                        %parameters are bad so reset and reduce step size
                        Dlt(iP)=Dlt(iP).*0.1; %
                        Chg=Inf; %need another loop
                    else %good case
                        Dlt(iP)=Dlt0(iP); %reset step size on good step
                        Chg=Chg+sum(abs((newp(:)-p(:))));
                        p=newp; %set new parameters
                    end
                    
                end
                
                %fprintf('%g\n',Chg);
                %% Check parameters
                if Chg<tol %stopping criteria
                    %   fprintf('Converged in %g iterations\n',iIt)
                    break
                end
            end %iterations
            
            if obj.TwoParamFlag
                p=[p;zeros(1,obj.nDmn-1);ones(1,obj.nDmn-1)];
            end
        end   %fit_newtonraphson
        
        function p0=startingSolution(obj,X,Y,A,b0,Lam)
            % function to calculate starting solution
            % can do this for all dimensions at the same time
            %% Assume b=0 for starting solution
            if nargin<=4
                b0=zeros(1,obj.nAsc);
            end
            % update starting solution for the parameters
            p0=obj.updateStartingSolution(X,Y,A,b0);
            
            if obj.TwoParamFlag
                p0=p0(obj.PrmInd<=2);
            end
            
            % while loop for checking the hessians
            % if the eigenvalues aren't
            cnt=0;
            minEig=obj.eigenCheck(X,Y,A,p0,Lam);
            if any(minEig<=0)
                while any(minEig<=0)
                    cnt=cnt+1;
                    b0=zeros(obj.nPrm(2),obj.nAsc);
                    for iAsc=1:obj.nAsc
                        b0(:,iAsc)=b0(:,iAsc)-0.1; %TODO is this adjustment
                        % determine the starting solution for the parameters
                        p0=obj.updateStartingSolution(X,Y,A,b0);
                        % check the updated eigen values solution
                        minEig=obj.eigenCheck(X,Y,A,p0,Lam);
                    end %iD
                    if cnt==100 || min(b0(:)) < obj.MinBet
                        %if we struggle to find a good starting solution
                        %we stop at the original soln assuming b=0
                        b0(b0<obj.MinBet)=obj.MinBet;
                        p0=obj.updateStartingSolution(X,Y,A,b0);
                        return
                    end
                end
            end
            
        end %startingSolution
        
        function [Alp,Bet,M,S]=p2Prm(obj,p)
            %convert optimization p to separate parameters
            if obj.TwoParamFlag
                Alp=p(obj.PrmInd==1,:);
                Bet=p(obj.PrmInd==2,:);
                M=zeros(1,obj.nAsc-1);
                S=ones(1,obj.nAsc-1);
                
            else
                Alp=p(obj.PrmInd==1,:); %nBin x nD-1
                Bet=p(obj.PrmInd==2,:); %1 xnD-1
                M=p(obj.PrmInd==3,:);  %1 x nD -1
                S=p(obj.PrmInd==4,:); % 1 x nD-1 - standard deviation
            end
        end
        
        function PLOGL=likelihood(obj,X,Y,A,p,L)
            %function PLOGL=likelihood(obj,X,Y,A,p,L)
            %compute penalised Heffernan and Tawn likelihood for data on Gumbel scale.
            %INPUT
            %-X         n x 1 conditioned value
            %-Y         n x 1 x (nDmn-1) conditioning value
            %-A         n x 1 x (nDmn-1) bin allocation
            %-p         nBin+3 x (nDmn-1) parameter values
            %-L         smoothness parameter (alpha)
            %-MrgTp     Gumbel or Laplace
            %-TwoPrmFlag TwoParameters obj Model if true, else 4 parameter
            if nargin<=5
                L=0; %no penalty
            end
            
            if obj.checkParam(p)
                PLOGL=Inf;
                return
            end
            
            [tAlp,tBet,M,S]=obj.p2Prm(p);
            
            Alp_Dat = tAlp(A(:,1),:);
            Bet_Dat = tBet(A(:,2),:);
            M_Dat = M(A(:,3),:);
            S_Dat = S(A(:,4),:);  %B*Alp
            
            Xb=bsxfun(@power,X,Bet_Dat);  %X^b
            Std=bsxfun(@times,S_Dat,Xb); %standard deviation
            
            Y = reshape(Y,[],obj.nDmn-1);
            Kappa=sqrt(gamma(1./obj.Delta)./gamma(3./obj.Delta));
            StdTerm=bsxfun(@times,Std,Xb);
            NLOGL=sum(obj.n.*(-log(obj.Delta)+log(2)+gammaln(1./obj.Delta))+ sum((abs(Y-bsxfun(@times,Alp_Dat,X)-bsxfun(@times,M_Dat,Xb))./(Kappa.*StdTerm)).^obj.Delta+log(Kappa.*StdTerm)));
            if obj.nBin>1
                PLOGL=NLOGL; %Penalised Likelhiood
                for iP=1:4
                    if obj.nPrm(iP)>1
                        switch iP
                            case 1
                                PLOGL=PLOGL+L.*mean((tAlp(:)-mean(tAlp(:))).^2);
                            case 2
                                PLOGL=PLOGL+L.*mean((tBet(:)-mean(tBet(:))).^2);
                            case 3
                                PLOGL=PLOGL+L.*mean((M(:)-mean(M(:))).^2);
                            case 4
                                PLOGL=PLOGL+L.*mean((S(:)-mean(S(:))).^2);
                        end
                    end
                end
            else
                PLOGL=NLOGL;
            end
            if isnan(PLOGL)
                warning('NaNs in the likelihood');
            end
        end %likelihood
        
        function [G,H]=Gradient(obj,iP,X,Y,A,p,L)
            %[G,H]=Gradient(iP,X,Y,A,p,L,TwoPrmFlg)
            %compute gradients for the Heffernan and Tawn model
            %INPUT
            % iP        either 1 for [alp and mu] or 2 for [beta and sigma]
            %-X         n x 1 conditioned value
            %-Y         n x 1 x (nDmn-1) conditioning value
            %-A         n x 1 x (nDmn-1) bin allocation
            %-p         nBin+3 x (nDmn-1) parameter values
            %-L         smoothness parameter (alpha)
            %
            %OUTPUT
            % G       p x (nDmn-1)
            % H       p x p x (nDmn-1)
            
            if nargin<=5
                L=0; %no penalty
            end
            
            [Alp,Bet,M,S]=obj.p2Prm(p);
            
            if obj.nBin==1
                L=0; %switch off penalty
            end
            
            % Get paramters at data
            Alp_Dat = Alp(A(:,1),:);
            Bet_Dat = Bet(A(:,2),:);
            if obj.TwoParamFlag
                M_Dat=0;
                S_Dat=1;
            else
                M_Dat = M(A(:,3),:);
                S_Dat = S(A(:,4),:);  %B*Alp
            end
            Kappa=sqrt(gamma(1./obj.Delta)./gamma(3./obj.Delta));
            Xb=X.^Bet_Dat;  %X^b
            MeanTerm=Alp_Dat.*X+(Xb).*M_Dat; %mean term of the distribution
            StdTerm=Kappa.*(Xb).*S_Dat;%standard deviation of the distribtuion
            R=Y-MeanTerm;            
            %% First Derivatives
            switch iP
                case 1 %alpha and mu
                    G = NaN(sum(obj.nPrm([1,3])),obj.nDmn-1);
                    %dL/dalp
                    t1=-X.*obj.Delta.*R.*(abs(R).^(obj.Delta-2))./((StdTerm).^obj.Delta);
             
                    for iDmn=1:obj.nDmn-1
                        %Penalty term
                        if obj.NonStat(1)
                            P=Alp(:,iDmn);
                        else
                            P=0;
                        end    
                        G(1:obj.nPrm(1),iDmn) = accumarray(A(:,1),t1(:,iDmn),[obj.nPrm(1),1],@sum,0)+L.*P;
                    end
                    
                    %dL/dmu
                    if obj.TwoParamFlag
                        G(obj.nPrm(1)+1,:)=[];
                    else
                        t1=-Xb.*obj.Delta.*R.*abs(R).^(obj.Delta-2)./((StdTerm).^obj.Delta);
                        for iDmn=1:obj.nDmn-1
                            %Penalty term
                            if obj.NonStat(3)
                                P=M(:,iDmn);
                            else
                                P=0;
                            end   
                            G(obj.nPrm(1)+1:end,iDmn) = accumarray(A(:,3),t1(:,iDmn),[obj.nPrm(3),1],@sum,0)+L.*P;
                        end
                        
                    end
                case 2 %beta and sigma
                    G = NaN(sum(obj.nPrm([2,4])),obj.nDmn-1);
                    %dL/dbeta
                    t1=log(X)-obj.Delta.*R.*abs(R).^(obj.Delta-2).*M_Dat.*Xb.*log(X)./((StdTerm).^obj.Delta)...
                        -obj.Delta.*log(X).*abs(R).^(obj.Delta)./((StdTerm).^obj.Delta);
                    
                    for iDmn=1:obj.nDmn-1
                            %Penalty term
                            if obj.NonStat(2)
                                P=Bet(:,iDmn);
                            else
                                P=0;
                            end   
                        G(1:obj.nPrm(2),iDmn) = accumarray(A(:,2),t1(:,iDmn),[obj.nPrm(2),1],@sum,0)+L.*P;
                    end
                    
                    %dL/dsigma
                    if obj.TwoParamFlag
                        G(obj.nPrm(2)+1,:)=[];
                    else
                        t1=(1./S_Dat)-(obj.Delta.*abs(R).^(obj.Delta)./(S_Dat.*(StdTerm).^obj.Delta));
                        for iDmn=1:obj.nDmn-1
                            %Penalty term
                            if obj.NonStat(4)
                                P=S(:,iDmn);
                            else
                                P=0;
                            end   
                            G(obj.nPrm(2)+1:end,iDmn) = accumarray(A(:,4),t1(:,iDmn),[obj.nPrm(4),1],@sum,0)+L.*P;
                        end
                    end
            end
            
%             tDel=obj.Delta;
%             epsilon=1;%1e-1; %deal with delta==1 issues
%             %
%             if abs(obj.Delta-1)<epsilon
%                 tDel=1+epsilon;
%             end
            
            
            %% Expected Fisher information matrix
            if nargout==2
                H=NaN(size(G,1),size(G,1),size(G,2));
                switch iP
                    case 1 %alpha and mu
                        
                        for iDmn=1:obj.nDmn-1 %loop iDmn
                            t1=(1./(S_Dat.^2)).*X.^(2-2*Bet_Dat);
                            %Penalty term
                            if obj.NonStat(1)
                                P=1;
                            else
                                P=0;
                            end   
                            E11 = diag(accumarray(A(:,1),t1(:,iDmn),[obj.nPrm(1),1],@sum,0)+L.*P);
                            if  obj.TwoParamFlag
                                H(:,:,iDmn)=E11;
                            else
                                %t1=(tDel-1).*tDel.*X.^(2*Bet_Dat).*abs(R).^(tDel-2)./((StdTerm).^tDel);
                                t1=1./(S_Dat.^2);
                                %Penalty term
                                if obj.NonStat(3)
                                    P=1;
                                else
                                    P=0;
                                end   
                                E33=diag(accumarray(A(:,3),t1(:,iDmn),[obj.nPrm(3),1],@sum,0)+L.*P);
                                %t1=(tDel-1).*tDel.*X.^(Bet_Dat+1).*abs(R).^(tDel-2)./((StdTerm).^tDel);
                                t1=(X.^(1-Bet_Dat))./(S_Dat.^2);
                                if obj.NonStat(1) &&  obj.NonStat(3) %both non stationary
                                    E13 = diag(accumarray(A(:,1),t1(:,iDmn),[obj.nPrm(1),1],@sum,0)); % nBin  x nBin
                                elseif obj.NonStat(1)
                                    E13 = accumarray(A(:,1),t1(:,iDmn),[obj.nPrm(1),1],@sum,0); % nBin  x 1
                                elseif obj.NonStat(3)
                                    E13 = accumarray(A(:,3),t1(:,iDmn),[obj.nPrm(3),1],@sum,0)'; % 1 x nBin
                                else %both stationary
                                    E13 = sum(t1(:,iDmn)); % 1 x 1
                                end
                                H(:,:,iDmn)=[E11,E13;E13',E33]; %[nBin +1 x nBin +1]
                            end
                        end %iDmn
                    case 2 %beta and sigma
                        %d2L/DLbeta2
                        %t1=(log(X).^2).*(2 + M_Dat.^2./V_Dat);
                        %Y_m_aX=(Y-Alp_Dat.*X);
                        %t1=((obj.Delta.*(log(X).^2).*Y_m_aX)./((StdTerm).^obj.Delta)).*...
                        %		(obj.Delta.*(Y_m_aX)-Xb.*M_Dat).*abs(R).^(obj.Delta-2);
                        for iDmn=1:obj.nDmn-1 %loop iDmn
                            t1=(log(X).^2).*(2 + M_Dat.^2./S_Dat.^2);
                            %Penalty term
                            if obj.NonStat(2)
                                P=1;
                            else
                                P=0;
                            end  
                            E22 = diag((accumarray(A(:,2),(t1(:,iDmn)),[obj.nPrm(2),1],@sum,0))+L.*P);
                            
                            if  obj.TwoParamFlag
                                H(:,:,iDmn)=E22;
                            else
                                
                                %d2L/DLsigma2
                                %t1=(1./(2.*V_Dat(:,iDmn).^2));
                                %t1=(obj.Delta+1).*obj.Delta.*abs(R).^(obj.Delta)./(((StdTerm).^obj.Delta).*(S_Dat.^2))-1./(S_Dat.^2);
                                t1=(2./(S_Dat.^2));
                                %Penalty term
                                if obj.NonStat(4)
                                    P=1;
                                else
                                    P=0;
                                end  
                                E44=diag((accumarray(A(:,4),(t1(:,iDmn)),[obj.nPrm(4),1],@sum,0))+L.*P);
                                
                                %d2L/DLbetasigma
                                %t1=(1./V_Dat).*log(X);
                                %t1=zeros(size(X));%obj.Delta.^2.*log(X).*Y_m_aX.*R.*abs(R).^(obj.Delta-2)./((StdTerm).^(obj.Delta).*S_Dat);
                                t1=2.*log(X)./S_Dat;
                                if obj.NonStat(2) &&  obj.NonStat(4) %both non stationary
                                    E24 = diag(accumarray(A(:,2),(t1(:,iDmn)),[obj.nPrm(2),1],@sum,0)); % nBin  x nBin
                                elseif obj.NonStat(2)
                                    E24 = (accumarray(A(:,2),(t1(:,iDmn)),[obj.nPrm(2),1],@sum,0)); % nBin  x 1
                                elseif obj.NonStat(4)
                                    E24 = (accumarray(A(:,4),(t1(:,iDmn)),[obj.nPrm(4),1],@sum,0))'; % 1 x nBin
                                else %both stationary
                                    E24 = (sum(t1(:,iDmn))); % 1 x 1
                                end
                                
                                H(:,:,iDmn)=[E22,E24;E24',E44]; %[2 x 2]
                            end
                        end %iDmn
                end
            end
        end%Gradient

        function Bad=checkParam(obj,varargin)
            %Check parameter values are in domain
            %TODO add Keef constraints??
            Bad=false;
            
            [tAlp,tBet,~,S]=obj.p2Prm(varargin{:});
            
            %nD checks need unwrap Alp(:)
            if any(tAlp(:)>obj.AlphaLimit(2)) || any(tAlp(:)<obj.AlphaLimit(1)) || ...
                    any(S(:)<0) || any(tBet(:)>1)  || any(tBet(:)<(obj.MinBet-eps))
                Bad=true;
                return
            end
           
            if strcmp(obj.MarginType,'Gumbel')
                if any(tAlp(:)<0)
                    Bad=true;
                    return
                end
            end
        end %checkParam
        
        function minEig=eigenCheck(obj,X,Y,A,p0,Lam)
            % store min eigen values
            min_eigen=zeros(obj.nAsc,2);
            % test alpha and mu
            [~,H1]=obj.Gradient(1,X,Y,A,p0,Lam);
            % test beta and sigma
            [~,H2]=obj.Gradient(2,X,Y,A,p0,Lam);
            if any(isinf(H1(:))) || any(isinf(H2(:))) || any(isnan(H1(:))) || any(isnan(H2(:)))
                minEig=-inf;
            else
                for iD=1:obj.nDmn-1
                    if obj.nDmn-1 == 1
                        min_eigen(iD,1)=min(eig(H1(:,:)));
                        min_eigen(iD,2)=min(eig(H2(:,:)));
                    else
                        min_eigen(iD,1)=min(eig(H1(:,:,iD)));
                        min_eigen(iD,2)=min(eig(H2(:,:,iD)));
                    end
                end %iD
                minEig=min(min_eigen,[],2);
            end
        end %eigenCheck
        
        function p=updateStartingSolution(obj,X,Y,A,b0)
            % X input variables [nObs x 1]
            % Y associated variables [nObs x nAsc]
            % A bin allocation [nObs x 4]
            % b0 beta parameter [1 x nAsc]
            Tol=1e-6;
            %% Alp and Mu | b
            if any(obj.NonStat)
                B=A(:,find(obj.NonStat,1))==1:obj.nBin; % B is [nObs x obj.nBin]
            end
            % if alpha is non stationary
            if obj.NonStat(1)
                BX=B.*X; % BX is [nObs x obj.nBin]
            else
                BX=X;%BX is always [nObs x 1]
            end
            % pre allocation of parameters
            a0=NaN(obj.nPrm(1),obj.nAsc);
            b0=ones(obj.nPrm(2),obj.nAsc).*b0;
            m0=NaN(obj.nPrm(3),obj.nAsc);
            v0=ones(obj.nPrm(4),obj.nAsc);
            % empty residual vector
            R=zeros(size(Y));
            
            for iAsc=1:obj.nAsc
                % if mu is non stationary
                if obj.NonStat(3)
                    BXpb=B.*X.^b0(iAsc);% BXpb will be [nObs x obj.nBin]
                else
                    BXpb=X.^b0(iAsc); % BXpb will be [nObs x 1]
                end
                Q=[BX,BXpb];
                %RegTerm=Lam*eye(obj.nPrm(1)+obj.nPrm(3));
                RegTerm=Tol*eye(obj.nPrm(1)+obj.nPrm(3));
                pAlpMu=(Q'*Q + RegTerm)\Q'*Y(:,iAsc); % [(obj.nPrm(1) + obj.nPrm(3)) x 1] parameter estimates for alpha and mu
                
                a0(:,iAsc)=max(pAlpMu(1:obj.nPrm(1)),obj.AlphaLimit(1)+1e-6);
                a0(:,iAsc)=min(pAlpMu(1:obj.nPrm(1)),obj.AlphaLimit(2)-1e-6);                
                
                m0(:,iAsc)=pAlpMu(obj.nPrm(1)+1:end);
                %non stationarity in b
                if obj.NonStat(2)
                    %indexing over the different directional bins
                    diffY=(Y(:,iAsc)-Q*pAlpMu);
                    for iB = 1:obj.nBin
                        R(B(:,iB)==1,iAsc)=diffY(B(:,iB)==1)./X(B(:,iB)==1).^b0(iAsc);
                    end
                else
                    %stationary in b
                    R(:,iAsc)=(Y(:,iAsc)-Q*pAlpMu)./(X.^b0(iAsc));
                end
                %% variance | alp,mu,b
                %RT 08/01/2020 issues with the starting values if we use a
                %stationary variance across bins
                %scaling the variance by a factor also helps with the
                %convergence of the model
                if obj.NonStat(4)%TODO for either the variance and std
                    for i = 1:obj.nBin
                        %handle situations where we have bad data or there
                        %is no data in the bin
                        v0(i,iAsc)=std(R(B(:,i)==1,iAsc),1).*0.1;
                        if v0(i,iAsc)==0 || isempty(R(B(:,i)==1,iAsc))
                            v0(i,iAsc)=std(R(:,iAsc),1).*0.1;
                        end
                        
                    end
                end
            end
            % stationary case for v0
            if obj.NonStat(4)==0
                v0=ones(obj.nPrm(4),obj.nAsc).*std(R,1);
            end
            % concatenate all the parameters together
            p=[a0;b0;m0;v0];
            
        end %updateStartingSolution
                          
        function obj=RVSimulations(obj,Mrg)
            % Run conditional return value simulation
            fprintf('Computing Conditional return value\n')
            nRls=10000;
            % Monte carlo simulation from the fitted model
            %Compute conditional quantiles of Y|X for given X values
            obj.RVSml=obj.SimulateMC(Mrg,"retValX",nRls,(1:obj.nBin)',obj.nBin);
            
        end %RVSimulations
           
    end %Methods private
    
end %classdef

