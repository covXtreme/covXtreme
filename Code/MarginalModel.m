% Copyright © [2023] Shell Global Solutions International B.V. All Rights Reserved.
% SPDX-License-Identifier: Apache-2.0

classdef MarginalModel
    %Fit Piecewise constant EVA Model
    %
    % o Estimates piecewise constant threshold using local bin quantile given non-exceedance
    % probability (NEP).
    % o Annual rate calculated using local Possion MLE for each bin = (no.% obs)/(no. years of data)
    % o Generalised Pareto distribution fitted above the threshold:
    %     - Shape assumed constant across all bins
    %     - Scale can varies across covariate bins (constant within bins)
    %     - Variation in scale parameter across bins controlled by smoothness penalty
    %     estimated using cross validation
    
    properties
        %% Data
        X        %nDat x nCvr, covariate direction data (X) [0,360]
        Y        %nDat x 1, response data (Y)
        Yrs=1;   %1 x 1, number of years of data
        RspLbl      %1 x 1, (string), response label
        RspSavLbl   %1 x 1 (string), rseponse label for saving plots (Avoid special characters)
        CvrLbl      %nCvr x 1, (string), covariate label vecotr
        nBoot=100;  %1 x 1, number of bootstrap resamples
        RtrPrd=100; %nRtr x 1 Return Period  in years
        
        %% Parameters
        Bn   % Bin class defining binning properties
        Scl      %nBin x 1, Fitted Generalised Paraeto Scale nBin x nB
        Shp      %1 x nBoot, Fitted Generalised Paraeto Shape (stationary shape from PPC marginal fitting)
        %nBin x nBoot (nonstationary shape)
        Omg      %nBin x nBoot   gamma shape parameter
        Kpp      %nBin x nBoot   gamma scale parameter
        GmmLct   %nBin x nBoot   gamma location parameter
        
        NEP      %nBoot x 1, Non Exceedence Probability Quantile threshold
        Thr      %nBin x nBoot, Extreme value threshold nBin x nB
        Rat      %nBin x nBoot, count no. observations in each bin nBin x nB
        BSInd    %nDat x nBoot, bootstrap index;
        
        nCvr      % 1 x 1, number of covariates
        nDat      % 1 x 1, number of data obs         
        %% Return Value
        RVSml    %Return value simulation
        
        %% Margin choice 
        MarginType = 'Laplace' %Laplace or Gumbel
    end
    
    properties(Hidden=true)
        FigureFolder='Figures';  %folder to save figures into (default Figures)
        BnMax         %bin max (used in check for upper endpoint)
    end
    
    properties(Dependent)
        nRtr    % 1 x 1 number of return periods
    end
    
    properties (Access=protected)
        IsPrd   %nCvr x 1  flag for if it is a periodic covariate
        %% Cross Validation (CV) options
        CVMth=0;      %1 x 1 (boolean), CV method
        %=0 Cross Validate smoothness parameter on original dataset only (fast);
        %=1 Cross Validate smoothness for every bootstrap resample (slow),
        nCV=10;       %1 x 1, no. cross-validation groups
        nSmth=10;     %1 x 1, no. smoothness tried in CV
        SmthLB=-4;    %1 x 1, lower bound (log10) for smoothness range
        SmthUB=4;     %1 x 1, upper bound (log10) for smoothness range
        SmthSet       %nSmth x 1, Set of candidate smoothness parameters in CV
        OptSmth       %1 x 1, Optimal smoothness from SmthSet
        CVLackOfFit   %nSmth x nBoot, Lack of Fit of candidate smoothness parameterss in CV
        
    end
    
    methods
        function obj=MarginalModel(Dat,iDmn,NEP,Bn,nB,Yrs,RtrPrd,CV,MarginType)
            %% INPUTS:
            % - Dat structure  (from stage 1)
            %     - Dat.X     nDat x nCvr  covariate values
            %     - Dat.Y     nDat x nDmn  response data (Y)
            %     - Dat.IsPrd   nCvr x 1 vector if covairte is periodic or not
            %     - Dat.CvrLbl    char string for response label
            %     - Dat.RspLbl    char string for response label
            % - iDmn  Dimension of Y to use as reponse
            % - NEP   1x1 Non Exceedence Probability Quantile threshold
            % - Bn CovariateBinning class containing binning information
            %% OPTIONAL INPUTS:
            % - nB number of bootstrap resamples (assumed 1 if not specified)
            % - Yrs  number of years of data (assumed 1 if not specified)
            % - RtrPrd  return period (years) (assumed to be 100 if non speicifed)
            % - CV cross validation structure with control parameters
            % - MarginType, 'Gumbel', 'Laplace' (default is Laplace)
            %% OUTPUTS:
            % - obj, Marginal Model class containing details of data and
            % fitted piecewise constant marginal model
            
            if nargin == 0 %empty constructor
                return;
            end
            
            %% Check inputs
            [obj.nDat,obj.nCvr]=size(Dat.X);
            validateattributes(Dat.IsPrd,{'numeric','logical'},{'numel',obj.nCvr,'integer'},'BinAllocation','IsPrd',1);
            obj.IsPrd=Dat.IsPrd;
            
            for iC=1:obj.nCvr
                if obj.IsPrd(iC)  %periodic check
                    validateattributes(Dat.X(:,iC), {'numeric'},{'vector','<=',360,'>=',0},'MarginalModel','X',1);
                end
            end
            obj.X=Dat.X;
            % dimension check
            validateattributes(iDmn,{'numeric'},{'numel',1,'integer','positive'},'BinAllocation','Edg',2);
            % response data check
            validateattributes(Dat.Y(:,iDmn), {'numeric'},{'vector','numel', obj.nDat},'MarginalModel','Y',1);
            obj.Y=Dat.Y(:,iDmn);
            obj.RspLbl=Dat.RspLbl{iDmn};
            obj.RspSavLbl=Dat.RspLbl{iDmn};
            obj.CvrLbl=Dat.CvrLbl;          
            % validate non exceedance probability
            validateattributes(NEP, {'numeric'},{'<=',1,'>=',0},'MarginalModel','NEP',3);  %0<=Tau<=1 NEP range
            % checking that the previous stage has been run
            if ~isa(Bn,'CovariateBinning')
                error('Input Bn, should be of class: CovariateBinning');
            end
            obj.Bn=Bn;
            
            %% Optional inputs
            % validate number of bootstraps
            if nargin>=5
                validateattributes(nB, {'numeric'},{'scalar','nonnegative'},'MarginalModel','nB',5);
                obj.nBoot=nB;  %number of bootstraps
                if obj.nBoot==0  %deal with the case where no bootstraps are provided
                    obj.nBoot=1;
                end
            end
            obj.NEP = [range(NEP)/2+min(NEP) ;sort(rand(obj.nBoot-1,1)).*range(NEP)+min(NEP)];  
            %sample NEPs over range with middle of range first
            % validate number of years of data 
            if nargin>=6
                validateattributes(Yrs, {'numeric'},{'scalar','positive'},'MarginalModel','Yrs',6);
                obj.Yrs=Yrs;
            end
            % validate return period choice 
            if nargin>=7
                validateattributes(RtrPrd, {'numeric'},{'vector','positive'},'MarginalModel','RV',7);
                obj.RtrPrd=sort(RtrPrd);
                
            end
            % validate margin type 
            if nargin>=9
                validateattributes(MarginType,{'string','char'},{},'MarginalModel','MarginType',9)
                obj.MarginType = validatestring(MarginType,{'Gumbel','Laplace'});
            end
            
            %% Check cross-validation inputs
            if nargin>=8 %% Cross Validation parameters (used for generalied Pareto fitting
                if isfield(CV,'CVMth')
                    validateattributes(CV.CVMth, {'numeric'},{'binary'},'MarginalModel','CV.CVMth',8);
                    obj.CVMth=CV.CVMth;
                end
                if isfield(CV,'nCV')
                    validateattributes(CV.nCV, {'numeric'},{'scalar','positive'},'MarginalModel','CV.nCV',8);
                    obj.nCV=CV.nCV;      %number cross-validation groups
                end
                if isfield(CV,'nSmth')
                    validateattributes(CV.nSmth, {'numeric'},{'scalar','positive'},'MarginalModel','CV.nSmth',8);
                    obj.nSmth=CV.nSmth;    %number smoothnesses tried in CV
                end
                if isfield(CV,'SmthLB')
                    validateattributes(CV.SmthLB, {'numeric'},{'scalar'},'MarginalModel','CV.SmthLB',8);
                    obj.SmthLB=CV.SmthLB;   %lower bound (log10)  for smmothness range
                end
                if isfield(CV,'SmthUB')
                    validateattributes(CV.SmthUB, {'numeric'},{'scalar'},'MarginalModel','CV.SmthUB',8);
                    obj.SmthUB=CV.SmthUB;   %upper bound (log10)  for smmothness range
                end
            end
            
            
            %% Preallocate Output Arrays
            obj.Scl=NaN(obj.Bn.nBin,obj.nBoot);  %Fitted Generalised Pareto Scale [nBin x nB]
            obj.Shp=NaN(obj.nBoot,1);    %Fitted Generalised Paraeto Shape  [nB x 1]
            obj.Omg=NaN(obj.Bn.nBin,obj.nBoot); %Gamma Shape Parameter [nBin x nB]
            obj.Kpp=NaN(obj.Bn.nBin,obj.nBoot); %Gamma Scale Parameter [nBin x nB]
            obj.GmmLct=NaN(obj.Bn.nBin,1); %Gamma Location Parameter [nBin x 1]
            obj.Thr=NaN(obj.Bn.nBin,obj.nBoot); %Extreme value threshold [nBin x nB]
            obj.Rat=NaN(obj.Bn.nBin,obj.nBoot); %Annual Rate no. exceedence observations in each bin [nBin x nB]
            obj.OptSmth=NaN(obj.nBoot,1);  %Fitted smoothness parameter (smoothness in scale param across bins) [nB x 1]
            obj.CVLackOfFit=NaN(obj.nSmth,obj.nBoot); %Lack of fit associated with different smoothness parameters [nSmth x nB]
            obj.BSInd=NaN(numel(obj.Y),obj.nBoot);  %Bootstrap sample indices nDat x nBoot
            obj.SmthSet=logspace(obj.SmthLB,obj.SmthUB,obj.nSmth); %try range smoothness penalties for sigma varying by bin [1 x nSmth]
            
            %% Fit model
            obj = Fit(obj);
            
            %% Return Value
            RVMCRls=10000;
            obj.RVSml=sample_RV_MC(obj,RVMCRls,(1:Bn.nBin)',obj.Bn.nBin);
            %% Plots
            if ~exist('Figures','dir')
                mkdir('Figures');
            end
            Plot(obj)
            
        end %MarginalModel constructor
        
        function obj=Fit(obj)
            %obj=Fit(obj)
            %Bootstrap fit the threshold, rate and generalised Pareto model.
            %Generalised Pareto estimates smoothness of scale using k-fold cross
            %validation
            
            rng(1); %reset random seed to get same bootstrap sample for all marginals
            obj.BSInd=[(1:length(obj.Y))',randi(length(obj.Y),length(obj.Y),obj.nBoot-1)];
            rng('shuffle');
            obj.BnMax=accumarray(obj.Bn.A,obj.Y,[obj.Bn.nBin,1],@max,NaN);
            
            obj.GmmLct=accumarray(obj.Bn.A,obj.Y,[obj.Bn.nBin,1],@(x)min(x)-range(x)*0.001,0); %threshold (simple quantile in each bin)
            
            %% Bootstrap loop
            for iBt = 1:obj.nBoot %will end up with obj.nBoot+1 sample parameters (1 for orig data, obj.nBoot boostrap resampled data)
                if obj.nBoot>1
                    fprintf('Fitting for bootstrap sample %d of %d\n',iBt,obj.nBoot);
                else
                    fprintf('Fitting sample\n');
                end
                %% Bootstrap Allocation
                [tY,A]=GetBootstrapSample(obj,obj.BSInd(:,iBt));
                
                %% Fit gamma distribution
                for iBn=1:obj.Bn.nBin
                    I=A==iBn;
                    if any(I)
                        p=gamfit(tY(I)-obj.GmmLct(iBn));
                        % store parameters
                        obj.Omg(iBn,iBt)=p(1);
                        obj.Kpp(iBn,iBt)=p(1).*p(2); %orthogonal parameterisation.
                    end
                end
                
                %% Threshold
                obj.Thr(:,iBt)=obj.gaminv(obj.NEP(iBt),obj.Omg(:,iBt),obj.Kpp(:,iBt),obj.GmmLct); %threshold (simple quantile in each bin)
                
                %% Get Exceedences
                IExc=(tY>obj.Thr(A,iBt));   %Index of exceedances
                AExc=A(IExc);  %bins allocation of exceedances
                
                Rsd=tY(IExc)-obj.Thr(AExc,iBt);   %Y-u above threshold used for gpfit
                ObsMax=obj.BnMax(AExc)-obj.Thr(AExc,iBt); %observation max used in upper end point of gp
                
                %% Rate of occurence in each bin
                obj.Rat(:,iBt)=accumarray(A,A,[obj.Bn.nBin,1],@numel)./obj.Yrs;
                
                %% Generalised Pareto Fit
                obj=GPCrossValidation(obj,Rsd,ObsMax,AExc,iBt);
                
            end
        end %BootMargModel
        
        function [YMrg,YUnif]=Margins(obj,iBt)
            %% Transform response Y to Uniform margins (GP -> Unif -> Standard margins)
            %% INPUT
            % - (Optional) [1 x 1] iBt, index on bootstrap resample (default 1)
            %% OUTPUT
            % - [nDat x 1] YUnif,  response data (Y) on Uniform margins
            % - [nDat x 1] YMrg,  response data (Y) on Standard margins
            % (Gumbel/Laplace)
            
            if nargin==1
                iBt=1:obj.nBoot;  %default to orginal data
            end
            
            YUnif=NaN(numel(obj.Y),numel(iBt));
            YMrg=NaN(numel(obj.Y),numel(iBt));
            
            for i=1:numel(iBt)
                jBt=iBt(i);
                tY=obj.Y(obj.BSInd(:,jBt));
                tA=obj.Bn.A(obj.BSInd(:,jBt));
                
                if obj.Shp(jBt)<0
                    if max(tY-obj.Thr(tA,jBt)+obj.Scl(tA,jBt)./obj.Shp(jBt))>0
                        error('Upper end point invalid');
                    end
                end
                YUnif(:,i)=CDF(obj,tY,tA,jBt);
                
                %% transform from uniform to standard margins
                YMrg(:,i)=INV_Standard(obj,YUnif(:,i));
                
                J = YUnif(:,i) == 1;
                if any(J)
                    Q=SURVIVOR(obj,tY(J),tA(J),jBt);
                    YMrg(J,i)=INV_Standard_survivor(obj,Q);
                end
            end
            
            if any(isinf(YMrg(:)))
                error('Bad transformation to %s',obj.MarginType)
            end
        end %Margins
        
        function Sml=sample_MC(obj,nRls,A,nA)
            %% INPUTS
            % Simulate Monte Carlo draws from marginal model
            % nRls number of realisations 
            % A [index of subper bins] bins to simulate
            % nA number of super bins
            %
            %% OUTPUTS
            % Sml.A    [nA x nRls]
            % Sml.I    [nA x nRls]
            % Sml.Unf  [nA x nRls]
            % Sml.Org  [nA x nRls]
            
            %Input A and nA super bins definitions
            if nargin<=2
                A=ones(obj.Bn.nBin,1); %omni situation                
            end    
            uA=unique(A(A>0)); %find unique bins (handles rescricted domain case)         
            %Input 4 number of super bins
            if nargin<=3
                nA=max(A);
            end
            %% Choose bootstraps
            Sml.I=randi(obj.nBoot,1,nRls);%decide which bootstrap sample to use over all bootstraps HT.nBoot
            Sml.I=reshape(repmat(Sml.I,nA,1),[],1);
            %% Choose bins from possion rate model
            Sml.A=NaN(nA,nRls);
            Sml.nRls=nRls;
            % simulate bin and sector allocation
            if obj.Bn.nBin>1
                for iA=1:numel(uA)
                    tA=find(A==uA(iA));% local index within the loop
                    if numel(tA)==1
                        Sml.A(uA(iA),:)=tA;
                    else
                        %Simulate covariate with right rate
                        tRat=obj.Rat(tA,Sml.I);  %get rate of all observations
                        RatCdf=cumsum(tRat)/sum(tRat); %get rate cdf
                        tJ=sum(rand(1,nRls)>RatCdf,1)+1; %local index within the loop
                        Sml.A(uA(iA),:)=tA(tJ);%bin allocation
                    end
                end %iA
            else
                Sml.A=ones(1,nRls);
            end
            %Simulate data on uniform margins
            Sml.Unf=rand(nA,nRls);
            %Probability integral to get the data onto original margins 
            Sml.Org=reshape(obj.INV(Sml.Unf(:),reshape(Sml.I,[],1),reshape(Sml.A,[],1)),nA,nRls);
        end %sample_MC
        
        function Sml=sample_RV_MC(obj,nRls,A,nA)
            % Simulate Monte Carlo draws from return value distribution
            %% INPUTS 
            % nRls number of realisations
            % A [index of super bins] bins to simulate
            % nA number of super bins
            %
            %% OUTPUTS
            % Sml.A    [nA x nRls x nRtr]
            % Sml.I    [nA x nRls x nRtr]
            % Sml.Unf  [nA x nRls x nRtr]
            % Sml.Org  [nA x nRls x nRtr]
            
            %Input A and nA super bins definitions
            if nargin<=2
                A=ones(obj.Bn.nBin,1); %omni situation                
            end      
            uA=unique(A(A>0)); %find unique bins (handles restricted domain case)
            %Input 4 number of super bins
            if nargin<=3
                nA=max(A);
            end
            Sml.nRls=nRls;
            %% Choose bootstraps
            Sml.I=randi(obj.nBoot,1,nRls);%decide which bootstrap sample to use over all bootstraps HT.nBoot
            %preallocate
            Sml.A=NaN(nA,nRls,obj.nRtr);
            [Sml.Unf,Sml.Org]=deal(NaN(nA,nRls,obj.nRtr));
            
            for iRtr=1:obj.nRtr
                %adjust for return period  (inv of exp(-L*T*(1-C));
                rho=obj.Rat(:,Sml.I);%annual rate of occurence
                LT=rho*obj.RtrPrd(iRtr); %poisson Rate
                UX=rand(obj.Bn.nBin,nRls);
                U=1+log(UX)./(LT); %[nBin x nRls]
                U(U<0)=NaN; %non-occurence
                tOrgAllBin=obj.INV(U,Sml.I); %[nBin x nRls]
                for iA=1:numel(uA) %loop over each bin and take max                    
                    tA=find(A==uA(iA));% local index within the loop
                    [Sml.Org(uA(iA),:,iRtr),ind_max]=max(tOrgAllBin(tA,:),[],1);                    
                    Sml.A(uA(iA),:,iRtr)=tA(ind_max);
                    ind_jnt=sub2ind([obj.Bn.nBin,nRls],Sml.A(uA(iA),:,iRtr),1:nRls); %need joint index
                    Sml.Unf(uA(iA),:,iRtr)=U(ind_jnt);
                end
                % UX should be in the range [P0, 1] where P0 is the non occurence rate.
            end %iRtr
            
            Sml.I=repmat(Sml.I,nA,1,obj.nRtr);
        end %sample_RV_MC
        
        function Sml=sample_MC_CondX(obj,nRls,A,nA,RspCond)
            % Simulate Monte Carlo draws conditional on chosen input XCond
            %%
            %% INPUTS 
            % nRls number of realisations
            % A [index of super bins] bins to simulate
            % nA number of super bins
            % RspCond  [nA x nBoot x nRtr]
            %
            %% OUTPUTS
            % Sml.A    [nA x nRls x nRtr]
            % Sml.I    [nA x nRls x nRtr]
            % Sml.Unf  [nA x nRls x nRtr]
            % Sml.Org  [nA x nRls x nRtr]
            
            %Input A and nA super bins definitions
            if nargin<=2
                A=ones(obj.Bn.nBin,1); %omni situation
            end
            uA=unique(A(A>0)); %find unique bins (handles rescricted domain case where not all may exist!!)
            
            %Input 4 number of super bins
            if nargin<=3
                nA=max(A);
            end
            Sml.nRls=nRls;
            %decide which bootstrap sample to use over all bootstraps HT.nBoot
            Sml.I=randi(obj.nBoot,1,nRls);          
            tnRtr=size(RspCond,3);
            %preallocate
            [Sml.Unf,Sml.A]=deal(NaN(nA,nRls,tnRtr));
            Sml.Org=RspCond(:,Sml.I,:); %assign response 
            
            %% Sample bin with right rate condition on A
            for iA=1:numel(uA)
                tA=find(A==uA(iA));% local index within the loop
                if ~isempty(tA)
                    tP=NaN(numel(tA),obj.nBoot,tnRtr);
                    tf=NaN(numel(tA),obj.nBoot,tnRtr);
                    for iRtr=1:tnRtr
                        tP(:,:,iRtr)=obj.CDF(RspCond(uA(iA),:,iRtr),tA); %probability of the response
                        tf(:,:,iRtr)=obj.PDF(RspCond(uA(iA),:,iRtr),tA); %density of the response
                    end
                    W=cumsum(tf,1)./sum(tf,1); %relative weights across each bin
                    tJ=sum(rand(1,nRls)>W(:,Sml.I,:),1)+1; %sample according to rate.
                    
                    Sml.A(uA(iA),:,:)=tA(tJ);%bin allocation
                    for iRls=1:nRls
                        Sml.Unf(uA(iA),iRls,:)=tP(tJ(iRls),Sml.I(iRls),:);
                    end
                end
            end %iA
            
            Sml.I=repmat(Sml.I,nA,1,tnRtr);
        end %sample_MC_CondX
        
        function P=CDF(obj,X,A,I)
            %CDF function for marginal model  (empirical below threshold - GP above)
            %CDF(obj,X) compute CDF for all bins and bootstraps using original data X is the locations to compute cdf at
            %CDF(obj,X,A) compute CDF for all bins and bootstraps using original data at specific
            %bins A
            %CDF(obj,X,A,I) compute CDF for all bins and bootstraps using original data at specific
            %bins and bootstraps indexes I
            if nargin==3  %Y not specified but A specifed
                I=1:obj.nBoot;
            end
            
            if nargin<=2  %cdf for everything
                P=MarginalModel.gamgpcdf(X,obj.Shp',obj.Scl,obj.Thr,obj.Omg,obj.Kpp,obj.GmmLct,obj.NEP');
            else %specified specfic bins and bootstraps
                if numel(A)==numel(I) && numel(A)>1   %cdf for subset defined by bins (A) and bootstraps (I), where output is vector
                    if obj.Bn.nBin == 1  %single covariate bins
                        P=MarginalModel.gamgpcdf(X,obj.Shp(I),obj.Scl(I)',obj.Thr(I)',obj.Omg(I)',obj.Kpp(I)',obj.GmmLct(A),obj.NEP(I));
                    else  %multiple covariate bins
                        J=sub2ind([obj.Bn.nBin,obj.nBoot],A,I);
                        P=MarginalModel.gamgpcdf(X,obj.Shp(I),obj.Scl(J),obj.Thr(J),obj.Omg(J),obj.Kpp(J),obj.GmmLct(A),obj.NEP(I));
                    end
                else  % cdf for subset defined by bins (A) and bootstraps (I)
                    P=MarginalModel.gamgpcdf(X,obj.Shp(I)',obj.Scl(A,I),obj.Thr(A,I),obj.Omg(A,I),obj.Kpp(A,I),obj.GmmLct(A),obj.NEP(I)');
                end
            end
            
        end %CDF
        
        function P=SURVIVOR(obj,X,A,I)
            %SURVIVOR function for marginal model  (empirical below threshold - GP above)
            %SURVIVOR(obj,X) compute survivor function for all bins and bootstraps using original data X is the locations to compute cdf at
            %SURVIVOR(obj,X,A) compute survivor function for all bins and bootstraps using original data at specific
            %bins A
            %SURVIVOR(obj,X,A,I) compute survivor function for all bins and bootstraps using original data at specific
            %bins and bootstraps indexes I
            if nargin==3  %Y not specified but A specifed
                I=1:obj.nBoot;
            end
            
            if nargin<=2  %cdf for everything
                P=MarginalModel.gamgpsurvivor(X,obj.Shp',obj.Scl,obj.Thr,obj.Omg,obj.Kpp,obj.GmmLct,obj.NEP');
            else %specified specfic bins and bootstraps
                if numel(A)==numel(I) && numel(A)>1   %cdf for subset defined by bins (A) and bootstraps (I), where output is vector
                    if obj.Bn.nBin == 1  %single covariate bins
                        P=MarginalModel.gamgpsurvivor(X,obj.Shp(I),obj.Scl(I)',obj.Thr(I)',obj.Omg(I)',obj.Kpp(I)',obj.GmmLct(A),obj.NEP(I));
                    else  %multiple covariate bins
                        J=sub2ind([obj.Bn.nBin,obj.nBoot],A,I);
                        P=MarginalModel.gamgpsurvivor(X,obj.Shp(I),obj.Scl(J),obj.Thr(J),obj.Omg(J),obj.Kpp(J),obj.GmmLct(A),obj.NEP(I));
                    end
                else  % cdf for subset defined by bins (A) and bootstraps (I)
                    P=MarginalModel.gamgpsurvivor(X,obj.Shp(I)',obj.Scl(A,I),obj.Thr(A,I),obj.Omg(A,I),obj.Kpp(A,I),obj.GmmLct(A),obj.NEP(I)');
                end
            end
            
        end %SURVIVOR
        
        function P=PDF(obj,X,A,I)
            %PDF function for marginal model  (empirical below threshold - PDF above)
            %PDF(obj,X) compute PDF for all bins and bootstraps using original data X is the locations to compute cdf at
            %PDF(obj,X,A) compute PDF for all bins and bootstraps using original data at specific
            %bins A
            %PDF(obj,X,A,I) compute PDF for all bins and bootstraps using original data at specific
            %bins and bootstrap indexes I
            if nargin==3  %Y not specified but A specifed
                I=1:obj.nBoot;
            end
            
            if nargin<=2  %cdf for everything
                P=MarginalModel.gamgppdf(X,obj.Shp',obj.Scl,obj.Thr,obj.Omg,obj.Kpp,obj.GmmLct,obj.NEP');
            else %specified specifc bins and bootstraps
                if numel(A)==numel(I) && numel(A)>1
                    if obj.Bn.nBin == 1  %single covariate bins
                        P=MarginalModel.gamgppdf(X,obj.Shp(I),obj.Scl(I)',obj.Thr(I)',obj.Omg(I)',obj.Kpp(I)',obj.GmmLct(A),obj.NEP(I));
                    else  %multiple covariate bins
                        J=sub2ind([obj.Bn.nBin,obj.nBoot],A,I);
                        P=MarginalModel.gamgppdf(X,obj.Shp(I),obj.Scl(J),obj.Thr(J),obj.Omg(J),obj.Kpp(J),obj.GmmLct(A),obj.NEP(I));
                    end
                else  % cdf for subset defined by bins (A) and bootstraps (I)
                    P=MarginalModel.gamgppdf(X,obj.Shp(I)',obj.Scl(A,I),obj.Thr(A,I),obj.Omg(A,I),obj.Kpp(A,I),obj.GmmLct(A),obj.NEP(I)');
                end
            end
            
        end %PDF
        
        function logF=LogPDF(obj,X,A,I)
            %PDF function for marginal model  (empirical below threshold - PDF above)
            %PDF(obj,X) compute PDF for all bins and bootstraps using original data X is the locations to compute cdf at
            %PDF(obj,X,A) compute PDF for all bins and bootstraps using original data at specific
            %bins A
            %PDF(obj,X,A,I) compute PDF for all bins and bootstraps using original data at specific
            %bins and bootstrap indexes I
            if nargin==3  %Y not specified but A specifed
                I=1:obj.nBoot;
            end
            
            if nargin<=2  %cdf for everything
                logF=MarginalModel.loggamgppdf(X,obj.Shp',obj.Scl,obj.Thr,obj.Omg,obj.Kpp,obj.GmmLct,obj.NEP');
            else %specified specifc bins and bootstraps
                if numel(A)==numel(I) && numel(A)>1
                    if obj.Bn.nBin == 1  %single covariate bins
                        logF=MarginalModel.loggamgppdf(X,obj.Shp(I),obj.Scl(I)',obj.Thr(I)',obj.Omg(I)',obj.Kpp(I)',obj.GmmLct(A),obj.NEP(I));
                    else  %multiple covariate bins
                        J=sub2ind([obj.Bn.nBin,obj.nBoot],A,I);
                        logF=MarginalModel.loggamgppdf(X,obj.Shp(I),obj.Scl(J),obj.Thr(J),obj.Omg(J),obj.Kpp(J),obj.GmmLct(A),obj.NEP(I));
                    end
                else  % cdf for subset defined by bins (A) and bootstraps (I)
                    logF=MarginalModel.loggamgppdf(X,obj.Shp(I)',obj.Scl(A,I),obj.Thr(A,I),obj.Omg(A,I),obj.Kpp(A,I),obj.GmmLct(A),obj.NEP(I)');
                end
            end
            
        end %LogPDF
        
        function X=INV(obj,P,I,A)
            %Inverse CDF for marginal model  (empirical below threshold - GP above)
            %% INPUTS
            %P probability
            %I index of bootstraps to use
            %A index of bins to use
            %if I scalar --> case where finding inverse CDF in single bin
            %if I vector --> case where inverse CDF in across sampled bins and bootstraps
            %if P matrix --> case where finding inverse CDF for all bootstraps and bins
            %% OUTPUTS
            %X data on original margins 
            X=NaN(size(P));
            p=size(P);
            if numel(I)==1
                Cs=1; %scalar
            elseif prod(p(2:end))>1
                Cs=3; %matrix
            else
                Cs=2; %vector
            end
            if nargin<4
                Cs=3;
            end
            
            switch Cs
                case 1 %I scalar --> case where finding inverse CDF in single bin
                    X=MarginalModel.gamgpinv(P,obj.Shp(I),obj.Scl(A,I),obj.Thr(A,I),obj.Omg(A,I),obj.Kpp(A,I),obj.GmmLct(A,I),obj.NEP(I));
                    
                case 2  %I vector --> case where inverse CDF in across sampled bins and bootstraps                  
                    if obj.Bn.nBin==1
                        X=MarginalModel.gamgpinv(P,obj.Shp(I),obj.Scl(I)',obj.Thr(I)',obj.Omg(I)',obj.Kpp(I)',obj.GmmLct,obj.NEP(I));
                    else
                        J=sub2ind([obj.Bn.nBin,obj.nBoot],A,I);
                        X=MarginalModel.gamgpinv(P,obj.Shp(I),obj.Scl(J),obj.Thr(J),obj.Omg(J),obj.Kpp(J),obj.GmmLct(A),obj.NEP(I));
                    end
                    
                case 3 %I matrix --> case where finding inverse CDF for all bootstraps and bins
                    X=MarginalModel.gamgpinv(P,obj.Shp(I)',obj.Scl(:,I),obj.Thr(:,I),obj.Omg(:,I),obj.Kpp(:,I),obj.GmmLct,obj.NEP(I)');
            end
        end %INV
        
        function X=INV_survivor(obj,Q,I,A)
            %Inverse CDF for marginal model  (empirical below threshold - GP above)
            %% INPUTS
            %P probability
            %I index of bootstraps to use
            %A index of bins to use
            %if I scalar --> case where finding inverse CDF in single bin
            %if I vector --> case where inverse CDF in across sampled bins and bootstraps
            %if P matrix --> case where finding inverse CDF for all bootstraps and bins
            %% OUTPUTS
            %X data on original margins 
            X=NaN(size(Q));
            p=size(Q);
            if numel(I)==1
                Cs=1; %scalar
            elseif prod(p(2:end))>1
                Cs=3; %matrix
            else
                Cs=2; %vector
            end
            if nargin<4
                Cs=3;
            end
            
            switch Cs
                case 1 %I scalar --> case where finding inverse CDF in single bin
                    X=MarginalModel.gamgpinvsurvivor(Q,obj.Shp(I),obj.Scl(A,I),obj.Thr(A,I),obj.Omg(A,I),obj.Kpp(A,I),obj.GmmLct(A,I),obj.NEP(I));
                    
                case 2  %I vector --> case where inverse CDF in across sampled bins and bootstraps                  
                    if obj.Bn.nBin==1
                        X=MarginalModel.gamgpinvsurvivor(Q,obj.Shp(I),obj.Scl(I)',obj.Thr(I)',obj.Omg(I)',obj.Kpp(I)',obj.GmmLct,obj.NEP(I));
                    else
                        J=sub2ind([obj.Bn.nBin,obj.nBoot],A,I);
                        X=MarginalModel.gamgpinvsurvivor(Q,obj.Shp(J),obj.Scl(J),obj.Thr(J),obj.Omg(J),obj.Kpp(J),obj.GmmLct(I),obj.NEP(I));
                    end
                    
                case 3 %I matrix --> case where finding inverse CDF for all bootstraps and bins
                    X=MarginalModel.gamgpinvsurvivor(Q,obj.Shp(I)',obj.Scl(:,I),obj.Thr(:,I),obj.Omg(:,I),obj.Kpp(:,I),obj.GmmLct,obj.NEP(I)');
            end
        end %INV
        
        function X=INV_Standard(obj,P)
            %transform from uniform to standard margins
            %using inverse CDF
            switch obj.MarginType
                case 'Gumbel'
                    X = -log(-log(P));
                case 'Laplace'
                    X = sign(0.5-P).*log(2*min(1-P,P));
                otherwise
                    error('Margin Type not recognised')
            end
            
        end %INV_Standard
        
        function X=INV_Standard_survivor(obj,Q)
            %transform (1-uniform) to standard margins
            %using inverse CDF
            switch obj.MarginType
                case 'Gumbel'
                    X = -log(-log1p(-Q));
                case 'Laplace'
                    X = sign(Q-0.5).*log(2*min(Q,1-Q));
                otherwise
                    error('Margin Type not recognised')
            end
            
        end %INV_Standard
        
        function P=CDF_Standard(obj,X)
            %transform from standard to uniform margins using CDF
            switch obj.MarginType
                case 'Gumbel'
                    P = exp(-exp(-X));
                case 'Laplace'
                    P = (X>0)-.5*sign(X).*exp(-abs(X));
                otherwise
                    error('Margin Type not recognised')
            end
        end %CDF_Standard
        
        function P=Survivor_Standard(obj,X)
            %transform from standard to uniform margins using CDF
            switch obj.MarginType
                case 'Gumbel'
                    P = -expm1(-exp(-X));
                case 'Laplace'
                    P = (X<=0) + 0.5*sign(X).*exp(-abs(X));
                otherwise
                    error('Margin Type not recognised')
            end
        end %CDF_Standard
        
        function F=PDF_Standard(obj,X)
            %get probability density on standard margins
            switch obj.MarginType
                case 'Gumbel'
                    F = exp(-(X +exp(-X)));
                case 'Laplace'
                    F =0.5.*exp(-abs(X));
                otherwise
                    error('Margin Type not recognised')
            end
            F(isnan(F))=0;
        end %PDF_Standard
        
        function F=LogPDF_Standard(obj,X)
            %get log probability density on standard margins
            switch obj.MarginType
                case 'Gumbel'
                    F = -(X +exp(-X));
                case 'Laplace'
                    F =log(0.5)-abs(X);
                otherwise
                    error('Margin Type not recognised')
            end
            F(isnan(F))=0;
        end %LogPDF_Standard
        
        function Plot(obj)
            %Plot results of PPC marginal fit
            
            %% Transform Orignal Sample to Standard Margins
            figure(1);
            clf;
            obj.PlotMarginCheck;
            
            %% Diagnostic plots for GP fit
            %Piecewise Sigma: fitted GP scale (black)
            figure(2)
            clf;
            obj.PlotNonStatParam;
            
            %Constant Xi: Fitted GP shape (black)  and true shape (green)
            figure(3);
            clf
            obj.PlotGPShape;
            
            %% Lack of fit against different smoothness Lambda
            if ~all(isnan(obj.CVLackOfFit(:))) && obj.Bn.nBin > 1
                figure(4);
                clf;
                obj.PlotSmoothness;
            end
            
            %% Q-Q plots for each bin
            obj.PlotQQ;
            
            %% Threshold uncertianty
            if obj.nBoot>1 && numel(unique(obj.NEP))>1
                figure(7);
                clf;
                obj.PlotThresholdStability;
            end
            
            %% Return Value CDF plot
            figure(8);
            clf;
            obj.PlotRV;
            
        end %Plot
        
        function nRtr=get.nRtr(obj)
            nRtr=numel(obj.RtrPrd);
        end %get.nRtr
        
    end %methods
    
    methods (Access = private)
        
        function PlotMarginCheck(obj)
            %Plot of data on standard and uniform margin as check of the
            %goodness of fit
            
            %Plot raw data on all margins you have calculated so far
            mThr=median(obj.Thr,2);
            IExc=(obj.Y>mThr(obj.Bn.A));   %Index of exceedenses
            
            [YMrg,YUnif]=Margins(obj,1);
            nSubPlt = 1+2*(~isempty(YUnif)); %number of subplots; if have transformed raw data onto standard margins,
            
            for iC=1:obj.Bn.nCvr
                subplot(obj.Bn.nCvr,nSubPlt,1+3*(iC-1))
                plot(obj.X(IExc,iC),obj.Y(IExc),'k.')
                hold on
                plot(obj.X(~IExc,iC),obj.Y(~IExc),'.','color',[1,1,1]*0.7)
                if iC==1
                    title(sprintf('%s: Raw data', obj.RspLbl))
                end
                
                PlotBinEdge(obj.Bn,iC);
                PlotParameter(obj.Bn,obj.GmmLct,iC,'color','r','linewidth',2);
                PlotParameter(obj.Bn,obj.Thr,iC,'color','b','linewidth',2);
                xlabel(obj.CvrLbl(iC))
                ylabel(obj.RspLbl)
                axis tight
                grid on

                if nSubPlt > 1
                    %Data on uniform margins
                    subplot(obj.Bn.nCvr,nSubPlt,2+3*(iC-1))
                    plot(obj.X(:,iC),YUnif,'k.')
                    hold on
                    if iC==1
                        title(sprintf('On uniform margins \n NEP = %.2f', obj.NEP(1)))
                    end
                    xlabel(obj.CvrLbl(iC))
                    ylabel(sprintf('%s-uniform scale',obj.RspLbl))
                    
                    PlotBinEdge(obj.Bn,iC);
                    
                    set(gca,'xtick',0:45:360)
                    ylim([0,1])
                    xlim([0,360])
                    grid on;
                    
                    %Data on standard margins
                    subplot(obj.Bn.nCvr,nSubPlt,3*iC)
                    plot(obj.X(:,iC),YMrg,'k.')
                    hold on
                    if iC==1
                        title(sprintf('On %s margins \n NEP = %.2f',obj.MarginType, obj.NEP(1)))
                    end
                    xlabel(obj.CvrLbl(iC))
                    ylabel(sprintf('%s-%s scale',obj.RspLbl,obj.MarginType))
                    
                    PlotBinEdge(obj.Bn,iC)
                    
                    xlim([0,360])
                    set(gca,'xtick',0:45:360)
                    grid on;
                end
            end
            savePics(fullfile(obj.FigureFolder,sprintf('Stg3_%s_1_DataTransform',obj.RspSavLbl)))
        end %Plot Margin Check
        
        function PlotNonStatParam(obj)
            %Piecewise Sigma: fitted GP scale (black)
            for iC=1:obj.Bn.nCvr
                %GP Scale
                subplot(obj.Bn.nCvr,3,(iC-1)*3+1)
                if obj.Bn.nBin > 1  %if non-stationary, plot as function of covariate
                    PlotParameter(obj.Bn,obj.Scl,iC,'color','k','linewidth',2);
                    PlotBinEdge(obj.Bn,iC);
                    PlotParameter(obj.Bn,obj.Scl,iC,'color','k','linewidth',2);
                    xlabel(obj.CvrLbl(iC))
                    ylabel('\sigma')
                else   %if stationary, histogram
                    histogram(obj.Scl,'edgecolor','none','facecolor','k','normalization','pdf')
                    xlabel('\sigma')
                    ylabel('Empirical density');
                end
                title(sprintf('%s: GP scale',obj.RspLbl))
                grid on;
                
                %Gamma shape
                subplot(obj.Bn.nCvr,3,(iC-1)*3+2)
                if obj.Bn.nBin > 1
                    PlotParameter(obj.Bn,obj.Omg,iC,'color','k','linewidth',2);
                    PlotBinEdge(obj.Bn,iC);
                    PlotParameter(obj.Bn,obj.Omg,iC,'color','k','linewidth',2);
                    xlabel(obj.CvrLbl(iC))
                    ylabel('\omega')
                else
                    histogram(obj.Omg,'edgecolor','none','facecolor','k','normalization','pdf')
                    xlabel('\omega')
                    ylabel('Empirical density');
                end
                title(sprintf('%s: Gamma shape',obj.RspLbl))
                grid on;
                 
                %Gam scale
                subplot(obj.Bn.nCvr,3,(iC-1)*3+3)
                if obj.Bn.nBin > 1
                    PlotParameter(obj.Bn,obj.Kpp,iC,'color','k','linewidth',2);
                    PlotBinEdge(obj.Bn,iC);
                    PlotParameter(obj.Bn,obj.Kpp,iC,'color','k','linewidth',2);
                    xlabel(obj.CvrLbl(iC))
                    ylabel('\kappa')
                else
                    histogram(obj.Kpp,'edgecolor','none','facecolor','k','normalization','pdf')
                    xlabel('\kappa')
                    ylabel('Empirical density');
                end
                title(sprintf('%s: Gamma scale ',obj.RspLbl))
                grid on;
                
            end
            
            savePics(fullfile(obj.FigureFolder,sprintf('Stg3_%s_2_Parameters',obj.RspSavLbl)))
        end %PlotNonStatParam
        
        function PlotGPShape(obj)
            %Constant Xi: Fitted GP shape (black)
            histogram(obj.Shp,'edgecolor','none','facecolor','k','normalization','pdf')
            xlabel('\xi')
            ylabel('Empirical density');      
            title(sprintf('%s: GP shape',obj.RspLbl))
            grid on;
            
            savePics(fullfile(obj.FigureFolder,sprintf('Stg3_%s_3_ParametersShape',obj.RspSavLbl)))
        end %PlotGPShape
        
        function PlotSmoothness(obj)
            %Smoothness parameter plot
            if obj.nBoot>1 && obj.CVMth==1
                plot(obj.SmthSet,median(obj.CVLackOfFit,2,'nanflag','omitnan'),'k-','linewidth',2)
                hold on
                plot(obj.SmthSet,quantile(obj.CVLackOfFit,0.025,2),'k--','linewidth',2)
                plot(obj.SmthSet,quantile(obj.CVLackOfFit,0.975,2),'k--','linewidth',2)
            else
                plot(obj.SmthSet,obj.CVLackOfFit(:,1),'k-','linewidth',2)
            end
            axis tight
            hold on
            grid on
            plot(median(obj.OptSmth)*[1,1],ylim,'r--','linewidth',2)
            ylabel('Lack of fit')
            xlabel('\lambda')
            set(gca,'xscale','log');
            title(sprintf('%s: GP cross-validation lack of fit',obj.RspLbl))
            savePics(fullfile(obj.FigureFolder,sprintf('Stg3_%s_4_CV',obj.RspSavLbl)))
        end %PlotSmoothness
        
        function PlotQQ(obj)
            %QQ plot of goodness of fit of the model
            nQ=100;
            Q=permute(linspace(min(obj.Y),max(obj.Y)*1.2,nQ),[1,3,2]);
            C=MarginalModel.gamgpcdf(Q,obj.Shp',obj.Scl,obj.Thr,obj.Omg,obj.Kpp,obj.GmmLct,obj.NEP');
            Q=squeeze(Q);
            C=permute(C,[3,1,2]);
            
            if obj.Bn.nBin > 1
                figure(5);
                clf
                
                nPlt1=ceil(sqrt(obj.Bn.nBin)); %max size [nPlt x nPlt]
                nPlt2=ceil(obj.Bn.nBin./nPlt1);
                
                if obj.nBoot>1
                    qC=quantile(C,[0.025,0.5,0.975],3);
                end
                
                for iB=1:obj.Bn.nBin

                    subplot(nPlt2,nPlt1,iB)
                    
                    I=obj.Bn.A==iB;
                    if any(I)
                        P = (0:sum(I)-1)./sum(I);                        
                        plot(sort(obj.Y(I)),log10(1-P),'r.') %data
                        axis tight
                        hold on
                        if obj.nBoot>1
                            plot(Q,log10(1-qC(:,iB,2)),'k-')   %fitted CDF
                            plot(Q,log10(1-qC(:,iB,1)),'k--')
                            plot(Q,log10(1-qC(:,iB,3)),'k--')
                        else
                            plot(Q,log10(1-C(:,iB)),'k-')
                        end
                        xlim([min(obj.Y(I)),max(obj.Y(I))])
                        
                        title(sprintf('%s: Bin %s, nExc %g',obj.RspLbl,obj.Bn.BinLbl{iB},sum(I)))
                        
                        ylabel('log(1-p)')
                        xlabel(obj.RspLbl)
                        grid on;
                    else %if no data in bin, leave plot empty
                        title(sprintf('%s: Bin %s, nExc %g',obj.RspLbl,obj.Bn.BinLbl{iB},sum(I)))
                        box on
                    end
                end
                hold off;
                savePics(fullfile(obj.FigureFolder,sprintf('Stg3_%s_5_SectorGoodnessOfFit',obj.RspSavLbl)))
            end
            
            figure(6)
            clf;
            P=(0:obj.nDat-1)./obj.nDat;
            plot(sort(obj.Y),log10(1-P),'r.')
            axis tight
            hold on
            w=bsxfun(@rdivide,obj.Rat,sum(obj.Rat,1));  %propn of exceedence data falling in each bin
            COmni=sum(bsxfun(@times,shiftdim(w,-1),C),2); %sum (prob in bin .* CDF(:,iBin))
            
            COmni(COmni>=1)=1;
            if obj.nBoot>1
                qCOmni=squeeze(quantile(COmni,[0.025,0.5,0.975],3));
                plot(Q,log10(1-qCOmni(:,2)),'k-')
                hold on
                plot(Q,log10(1-qCOmni(:,[1,3])),'k--')
            else
                plot(Q,log10(1-COmni),'k-')
            end
            xlim([0,max(obj.Y)])
            title(sprintf('%s: Overall goodness of fit',obj.RspLbl))
            ylabel('log(1-p)')
            grid on
            xlabel(obj.RspLbl)
            savePics(fullfile(obj.FigureFolder,sprintf('Stg3_%s_6_OverallGoodnessOfFit',obj.RspSavLbl)))
        end %PlotQQ
        
        function PlotThresholdStability(obj)
            %plot of threshold stability
            if obj.nBoot>50
                nNEP = 10; %number of bins for threshold diagnostic plot
                NEPEdg=linspace(min(obj.NEP),max(obj.NEP),nNEP+1);
                NEPBin=(NEPEdg(1:end-1)+NEPEdg(2:end))/2;
                
                A = discretize(obj.NEP,NEPEdg);
                
                S=accumarray(A,obj.Shp,[nNEP,1],@median,NaN);
                Slb=accumarray(A,obj.Shp,[nNEP,1],@(x)quantile(x,0.025),NaN);
                Sub=accumarray(A,obj.Shp,[nNEP,1],@(x)quantile(x,0.975),NaN);
                NanI=isnan(S);
                plot(NEPBin(~NanI),S(~NanI),'r-','linewidth',2)
                hold on
                plot(NEPBin(~NanI),Slb(~NanI),'r--','linewidth',2)
                plot(NEPBin(~NanI),Sub(~NanI),'r--','linewidth',2)
            end
            plot(obj.NEP,obj.Shp,'k.','markersize',20)
            xlabel('NEP')
            ylabel('GP shape \xi')
            grid on
            title(sprintf('%s: GP shape stability by threshold',obj.RspLbl))
            savePics(fullfile(obj.FigureFolder,sprintf('Stg3_%s_7_ThresholdStability',obj.RspSavLbl)))
        end %PlotThresholdStability
        
        function PlotRV(obj)
            %Plot Return Value CDF's
            ColMat=hsv(obj.Bn.nBin);
            grid on
            for iRtr=1:obj.nRtr
                subplot(1,obj.nRtr,iRtr)
                hold on
                if obj.Bn.nBin > 1
                    for iBin = 1:obj.Bn.nBin
                        tX=sort(obj.RVSml.Org(iBin,:,iRtr),'missingplacement','first');
                        plot(tX,linspace(0,1,numel(tX)),'color',ColMat(iBin,:),'linewidth',2)
                        %hold on
                    end
                end
                XOmni=sort(max(obj.RVSml.Org(:,:,iRtr),[],1)); %take max over bins
                plot(XOmni,linspace(0,1,numel(XOmni)),'k','linewidth',2)

                %20230426 Phil tidy up xlim (for consistency with conditional return value plot)
                talPrb=1e-3;
                lb=min(quantile(obj.RVSml.Org(:,:,iRtr),talPrb,2));
                ub=max(quantile(obj.RVSml.Org(:,:,iRtr),1-talPrb,2));
                xlim([lb ub]);

                plot(xlim,[0.5,0.5],'--k')
                plot(xlim,[exp(-1),exp(-1)],'--k')

                ylabel('Cumulative probability')
                xlabel(obj.RspLbl)
                
                if obj.Bn.nBin >1 && iRtr==1
                    legend([obj.Bn.BinLbl(:);'Omni'],'location','best','fontsize',6);
                end
                
                title(sprintf('Maximum over %g years',obj.RtrPrd(iRtr)))
                grid on;
                box on;
            end
            savePics(fullfile(obj.FigureFolder,sprintf('Stg3_%s_8_ReturnValueCDF',obj.RspSavLbl)))
        end %PlotRV
        
        function [Y,A,I]=GetBootstrapSample(obj,I)
            % if iBt==1 use original sample
            Y = obj.Y(I);
            A = obj.Bn.A(I);
        end %GetBootstrapSample
        
        function obj=GPCrossValidation(obj,Rsd,ObsMax,AExc,iBt)
            % fit piecewise constant GP to threshold exceedances
            % INPUTS:
            % Rsd - Y-u threshold exceedances
            % AExc - bin allocation of exceedances
            % ObsMax - observation max used in the upper end point of the GP 
            % iBt - bootstrap index
            % OUTPUT:
            % fitted GP model with optimal smoothness 
            
            %% Constant Starting Solution (method of moments)
            xbar=mean(Rsd);
            s2=var(Rsd);
            xi0 =0;
            sigma0 = .5 .* xbar .* (xbar.^2 ./ s2 + 1)*2;
            p0=[xi0;log(ones(obj.Bn.nBin,1)*sigma0)];
            
            nExc=length(Rsd);
            opts=optimset('display','off');
            
            %% Flag for do Cross Validation
            if (obj.CVMth==0 &&  iBt>1) || obj.nSmth==1
                %% Cross validation off
                if obj.nSmth==1  %only one smoothness
                    obj.OptSmth(iBt)=obj.SmthSet;
                else  %use first bootstrap instead of cross validating every time
                    obj.OptSmth(iBt)=obj.OptSmth(1);
                end
            else
                fprintf('Starting Cross Validation to estimate GP Scale smoothness:\n');
                LackOfFit=NaN(obj.nCV,obj.nSmth); %lack of fit
                ICV=randi(obj.nCV,nExc,1); %split nExc observations into nCV cross validation groups
                
                for iCV=1:obj.nCV  %cross validation loop
                    fprintf('.')
                    %Pull out fit set
                    ObsMaxFit=  ObsMax(ICV~=iCV);
                    RsdFit=Rsd(ICV~=iCV);
                    AExcFit=AExc(ICV~=iCV);
                    %Pull out prediction set
                    ObsMaxPrd=  ObsMax(ICV==iCV);
                    RsdPrd=Rsd(ICV==iCV);
                    AExcPrd=AExc(ICV==iCV);
                    for iLam=1:obj.nSmth %Smoothness penalties
                        InitialNLL=MarginalModel.GPLikeNonStationary(p0,RsdFit,AExcFit,ObsMaxFit,obj.SmthSet(iLam));
                        if isinf(InitialNLL)
                            error('bad starting value')
                        end
                        PrmHat=fminsearch(@(p)MarginalModel.GPLikeNonStationary(p,RsdFit,AExcFit,ObsMaxFit,obj.SmthSet(iLam)),p0,opts);
                        LackOfFit(iCV,iLam)=MarginalModel.GPLikeNonStationary(PrmHat,RsdPrd,AExcPrd,ObsMaxPrd);
                    end %SigPen
                end %CV
                fprintf('\n')
                obj.CVLackOfFit(:,iBt)=sum(LackOfFit,1)';
                [~,tI]=min(obj.CVLackOfFit(:,iBt));
                obj.OptSmth(iBt)=obj.SmthSet(tI)';
            end
            
            %% Fit Model using optimal smoothness
            InitialNLL=MarginalModel.GPLikeNonStationary(p0,Rsd,AExc,ObsMax,obj.OptSmth(iBt));
            if isinf(InitialNLL)
                error('bad starting value')
            end
            PrmHat = fminsearch(@(p)MarginalModel.GPLikeNonStationary(p,Rsd,AExc,ObsMax,obj.OptSmth(iBt)),p0,opts);
            obj.Shp(iBt)=PrmHat(1);
            obj.Scl(:,iBt)=exp(PrmHat(2:end));
            
        end %GPCrossValidation
        
    end %methods private
    
    methods (Static)
        function PLOGL=GPLikeNonStationary(PARAMS,Dat,BinAlc,ObsMax,SigPen)
            %INPUTS:
            % - [(nBins+1) x 1] PARAMS, piecewise constant marginal model
            % parameters(= [xi, sig_1 ,..., sig_nBins])
            % - [nDat x 1] Dat, GP data (= threshold exceedences - threshold)
            % - [nDat x 1] BinAlc, bin allocation index on data
            % - [1 x 1] SigPen, penalty imposing smoothness in GP scale across covariate
            % - max seen in  each bin used for global upper end point constraint
            %OUTPUT:
            % - [1 x 1] PLOGL, Penalised negative log likelihood
            
            if size(Dat,1) ~= size(BinAlc,1)
                error('Input dimension mismatch: dmn Dat =?= dmn BinAllocation')
            end
            
            if nargin<=4
                SigPen=0;  %unpenalised
            end
            Xi=PARAMS(1);
            Sig=exp(PARAMS(2:end));
            nBins =numel(Sig);
            
            NLOGL=0;
            for iBin=1:nBins %Additive likelihood over covariate bins
                I=BinAlc==iBin;
                if any(I)
                    tGPLike = gplike([Xi,Sig(iBin)],Dat(I));
                    if isinf(tGPLike) || any(Xi<-0.5) || (Xi<0 && any(ObsMax(I)>-Sig(iBin)./Xi))
                        NLOGL=Inf;
                    else
                        NLOGL=NLOGL+tGPLike;
                    end
                end
            end
            
            PLOGL=NLOGL+SigPen.*sum((Sig-mean(Sig)).^2);
            
        end %GPLikeNonStationary
        
        function X=gpinv(P,Xi,Sgm,Thr)
            %     X=gpinv(P,Xi,Sgm,Thr) returns the inverse of generalized Pareto (GP)
            %     cdf with tail index (shape) parameter Xi, scale parameter Sgm,
            %     and threshold (location) parameter Thr, evaluated at the values in X.
            %     The size of P is the common size of the input arguments.
            
            t1=bsxfun(@rdivide,bsxfun(@power,1-P,-Xi)-1,Xi);
            
            %%Gumbel case
            I=abs(Xi)<1e-5;
            if any(I(:))
                t1(I)=-log(1-P(I));
            end
            
            X=bsxfun(@plus,bsxfun(@times,Sgm,t1), Thr);
        end %gpinv
        
        function X=gpinvsurvivor(Q,Xi,Sgm,Thr)
            %     X=gpinv(P,Xi,Sgm,Thr) returns the inverse of generalized Pareto (GP)
            %     cdf with tail index (shape) parameter Xi, scale parameter Sgm,
            %     and threshold (location) parameter Thr, evaluated at the values in X.
            %     The size of P is the common size of the input arguments.
            
            t1=bsxfun(@rdivide,bsxfun(@power,Q,-Xi)-1,Xi);
            
            %%Gumbel case
            I=abs(Xi)<1e-5;
            if any(I(:))
                t1(I)=-log(Q(I));
            end
            
            X=bsxfun(@plus,bsxfun(@times,Sgm,t1), Thr);
        end %gpinv
        
        function F=gppdf(X,Xi,Sgm,Thr)
            %     F = gppdf(X,Xi,Sgm,Thr) returns the pdf of the generalized Pareto (GP)
            %     distribution with tail index (shape) parameter Xi, scale parameter Sgm,
            %     and threshold (location) parameter THETA, evaluated at the values in X.
            %     The size of P is the common size of the input arguments.
            
            Z=bsxfun(@rdivide,bsxfun(@minus,X,Thr),Sgm);
            t1=1+bsxfun(@times,Xi,Z);
            t1(t1<=0)=0;  %deal with upper end point of distribution
            
            F=bsxfun(@times,1./Sgm,bsxfun(@power,t1,-1./Xi-1));
            
            %% Gumbel case
            I=abs(Xi)<1e-5;
            if any(I(:))
                Xi=bsxfun(@times,Xi,ones(size(Z)));
                I=abs(Xi)<1e-5;
                Sgm=bsxfun(@times,Sgm,ones(size(Z)));
                
                F(I)=bsxfun(@times,1./Sgm(I),exp(-Z(I)));
            end
            
            F(bsxfun(@lt,X,Thr))=0; %density of zero below the threshold
            F(Z<0)=0;
            
        end %gppdf
        
        function logF=loggppdf(X,Xi,Sgm,Thr)
            %     logF = loggppdf(X,Xi,Sgm,Thr) returns the log pdf of the generalized Pareto (GP)
            %     distribution with tail index (shape) parameter Xi, scale parameter Sgm,
            %     and threshold (location) parameter THETA, evaluated at the values in X.
            %     The size of P is the common size of the input arguments.
            
            Z=(X-Thr)./Sgm;
            t1=1+bsxfun(@times,Xi,Z);
            t1(t1<=0)=0;  %deal with upper end point of distribution
            
            logF=-log(Sgm)-(1./Xi+1).*log(t1);
            
            %% Gumbel case
            I=abs(Xi)<1e-5;
            if any(I(:))
                Xi=Xi.*ones(size(Z));
                I=abs(Xi)<1e-5;
                Sgm=Sgm.*ones(size(Z));
                
                logF(I)=-log(Sgm(I))-Z(I);
            end
            
            logF(Z<0)=-Inf; %density of zero below the threshold
            
        end %loggppdf
        
        function P=gpcdf(X,Xi,Sgm,Thr)
            %     P = CumulativeDensityFunction(obj,X) returns the cdf of the generalized Pareto (GP)
            %     distribution with tail index (shape) parameter Xi, scale parameter Sgm,
            %     and threshold (location) parameter THEThrTA, evaluated at the values in X.
            %     The size of P is the common size of the input arguments.
            Z=bsxfun(@rdivide,bsxfun(@minus,X,Thr),Sgm);
            t1=1+bsxfun(@times,Xi,Z);
            t1(t1<=0)=0;  %deal with upper end point of distribtion
            P=1-bsxfun(@power,t1,-1./Xi);
            %% Gumbel case
            I=abs(Xi)<1e-5;
            if any(I(:))
                Xi=bsxfun(@times,Xi,ones(size(Z)));
                I=abs(Xi)<1e-5;
                
                P(I)=1-exp(-Z(I));
            end
            P(bsxfun(@le,X,Thr))=NaN; %density of zero below the threshold
            P(Z<=0)=NaN;
            
        end %gpcdf
        
        function Q=gpsurvivor(X,Xi,Sgm,Thr)
            %     Q=gpsurvivor(X,Xi,Sgm,Thr) returns the survivor functoin of the generalized Pareto (GP)
            %     distribution with tail index (shape) parameter Xi, scale parameter Sgm,
            %     and threshold (location) parameter Thr, evaluated at the values in X.
            %     The size of P is the common size of the input arguments.
            Z=bsxfun(@rdivide,bsxfun(@minus,X,Thr),Sgm);
            t1=1+bsxfun(@times,Xi,Z);
            t1(t1<=0)=0;  %deal with upper end point of distribtion
            Q=bsxfun(@power,t1,-1./Xi);
            %% Gumbel case
            I=abs(Xi)<1e-5;
            if any(I(:))
                Xi=bsxfun(@times,Xi,ones(size(Z)));
                I=abs(Xi)<1e-5;
                
                Q(I)=exp(-Z(I));
            end
            Q(bsxfun(@le,X,Thr))=NaN; %density of zero below the threshold
            Q(Z<=0)=NaN;
        end %gpsurvivor
        
        function F=gampdf(X,Alp,Bet,GmmLct)
            %     F = gampdf(X,Alp,Bet,GmmLct) returns the pdf of the gamma distribution using orthog0nal
            %     parameterisation density
            Z=bsxfun(@minus,X,GmmLct);
            Z(Z<0)=0;
            
            F = bsxfun(@times,(Bet./Alp).^(-Alp)./gamma(Alp),bsxfun(@power,Z,(Alp-1))).*...
                exp(-bsxfun(@times,Alp./Bet,Z));
        end %gampdf
        
        function logF=loggampdf(X,Alp,Bet,GmmLct)
            %     logF = loggampdf(X,Alp,Bet,GmmLct) returns the log pdf of the gamma distribution using orthog0nal
            %     parameterisation
            %density
            Z=bsxfun(@minus,X,GmmLct);
            Z(Z<0)=0;
            
            logF = -Alp.*(log(Bet)-log(Alp))-gammaln(Alp)+(Alp-1).*log(Z)-Alp.*Z./Bet;
        end %loggampdf
        
        function P=gamcdf(X,Alp,Bet,GmmLct)
            %     P = gpcdf(X,Alp,Bet,GmmLct) returns the cdf of the gamma distribution using orthognal
            %     parameterisation
            Z=bsxfun(@minus,X,GmmLct);
            Z(Z<0)=0;
            Z = bsxfun(@times,Alp./Bet,Z);
            P = bsxfun(@gammainc,Z, Alp);
        end %gamcdf
        
        function Q=gamsurvivor(X,Alp,Bet,GmmLct)
            %     Q=gamsurvivor(X,Alp,Bet,GmmLct) returns the survivor 
            %     function of the gamma distribution
            Z=bsxfun(@minus,X,GmmLct);
            Z(Z<0)=0;
            Z = bsxfun(@times,Alp./Bet,Z);
            Q = bsxfun(@(x,y)gammainc(x,y,'upper'),Z, Alp);
        end %gamsurvivor
        
        function X=gaminv(P,Alp,Bet,GmmLct)
            %     X = gaminv(P,Alp,Bet,GmmLct)returns the inverse the gamma distribution using orthognal
            %     parameterisation
            
            P=bsxfun(@times,P,ones(size(Alp)));
            Alp=bsxfun(@times,Alp,ones(size(P)));
            
            
            q=NaN(size(P));
            
            I=~isnan(P) & P>=1e-2;
            q(I) = gammaincinv(P(I),Alp(I));
            I=~isnan(P) & P<1e-2;
            q(I)= gammaincinv(1-P(I),Alp(I),'upper');
            
            X = bsxfun(@plus,bsxfun(@times,q,Bet./Alp),GmmLct);
        end %gaminv
        
        function X=gaminvsurvivor(Q,Alp,Bet,GmmLct)
            %     X = gaminv(P,Alp,Bet,GmmLct)returns the inverse the gamma distribution using orthognal
            %     parameterisation

            Q=bsxfun(@times,Q,ones(size(Alp)));
            Alp=bsxfun(@times,Alp,ones(size(Q)));
            
            
            q=NaN(size(Q));
            
            I=~isnan(Q) & Q<=1-1e-2;
            q(I) = gammaincinv(1-Q(I),Alp(I));
            I=~isnan(Q) & Q>1-1e-2;
            q(I)= gammaincinv(Q(I),Alp(I),'upper');
            
            X = bsxfun(@plus,bsxfun(@times,q,Bet./Alp),GmmLct);
        end %gaminv
        
        function P=gamgpcdf(X,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau)
            %     P=gamgpcdf(X,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau) returns the cdf of the gamma-gp distribution
            
            %threshold
            %           Thr=MarginalModel.gaminv(Tau,Alp,Bet,GmmLct);
            IBlw=X<Thr;
            %gamma part
            P1=MarginalModel.gamcdf(X,Alp,Bet,GmmLct);
            % mjj 08/11/18: I think tau scaling for Gamma was missing:
            % added it.
            %gp part
            P2=MarginalModel.gpcdf(X,Xi,Sgm,Thr).*(1-Tau)+Tau;
            %combine
            P=NaN(size(P1));
            P(IBlw)=P1(IBlw);
            P(~IBlw)=P2(~IBlw);
            
        end %gamgpcdf
        
        function Q=gamgpsurvivor(X,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau)
            %     Q=gamgpsurvivor(X,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau) returns 
            %     the survivor function of the gamma-gp distribution
            
            %threshold
            %           Thr=MarginalModel.gaminv(Tau,Alp,Bet,GmmLct);
            IBlw=X<Thr;
            %gamma part
            P1=MarginalModel.gamsurvivor(X,Alp,Bet,GmmLct);
            % mjj 08/11/18: I think tau scaling for Gamma was missing:
            % added it.
            %gp part
            P2=MarginalModel.gpsurvivor(X,Xi,Sgm,Thr).*(1-Tau);
            %combine
            Q=NaN(size(P1));
            Q(IBlw)=P1(IBlw);
            Q(~IBlw)=P2(~IBlw);
            
        end %gamgpsurvivor
        
        function f=gamgppdf(X,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau)
            % f=gamgpcdf(X,Xi,Sgm,Alp,Bet,GmmLct,Tau) returns the pdf of the gamma-gp distribution         
            %threshold
            IBlw=bsxfun(@le,X,Thr);
            %gamma part
            f1=MarginalModel.gampdf(X,Alp,Bet,GmmLct);
            %gp part
            f2=MarginalModel.gppdf(X,Xi,Sgm,Thr).*(1-Tau);
            %combine
            f=NaN(size(f1));
            f(IBlw)=f1(IBlw);
            f(~IBlw)=f2(~IBlw);
            
        end %gamgppdf
        
        function f=loggamgppdf(X,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau)
            % f=loggamgppdf(X,Xi,Sgm,Alp,Bet,GmmLct,Tau) returns the pdf of the gamma-gp distribution           
            %threshold
            IBlw=X<=Thr;
            %gamma part
            f1=MarginalModel.loggampdf(X,Alp,Bet,GmmLct);
            %gp part
            f2=MarginalModel.loggppdf(X,Xi,Sgm,Thr)+log((1-Tau));
            %combine
            f=NaN(size(f1));
            f(IBlw)=f1(IBlw);
            f(~IBlw)=f2(~IBlw);
            
        end %loggamgppdf
        
        function X=gamgpinv(P,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau)
            %X=gamgpinv(P,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau) returns the inv of the gamma-gp distribution
            
            %gamma part
            X1=MarginalModel.gaminv(P,Alp,Bet,GmmLct);
            IBlw=P<(Tau.*ones(size(X1)));
            
            %gp part (adjust probabilities)
            PGP=(P-Tau)./(1-Tau);
            PGP(PGP<0)=0;
            X2=MarginalModel.gpinv(PGP,Xi,Sgm,Thr);
            
            X=NaN(size(X1));
            X(IBlw)=X1(IBlw);
            X(~IBlw)=X2(~IBlw);
        end %gamgpinv
        
        function X=gamgpinvsurvivor(Q,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau)
            %X=gamgpinv(P,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau) returns the inv of the gamma-gp distribution
            
            P=1-Q;
            %gamma part
            X1=MarginalModel.gaminvsurvivor(Q,Alp,Bet,GmmLct);
            IBlw=P<(Tau.*ones(size(X1)));
            
            %gp part (adjust probabilities)
            QGP=Q./(1-Tau);
            X2=MarginalModel.gpinvsurvivor(QGP,Xi,Sgm,Thr);
            
            X=NaN(size(X1));
            X(IBlw)=X1(IBlw);
            X(~IBlw)=X2(~IBlw);
        end %gamgpinv
        
        
    end %methods
end %class
