classdef HeffernanTawn
    
    %functions for fitting the conditional extremes model of Heffernan and Tawn (2004) "A conditional approach to modelling multivariate extreme values)"
    
    properties        
        Alp   %nBin x nDmn x nBoot or 1 x nBoot for slope parameter of Hefferenan
        Bet   %nBin x nDmn x nBoot  
        Mu    %nBin x nDmn x nBoot 
        Sig   %nBin x nDmn x nBoot       
        Rsd   % nBoot x 1 (cell array), residuals from each bootstrap
        Thr   % nBoot x 1, H&T conditional threshold
        NEP   % 1 x 1, non exceedence probability
        nBoot % 1 x 1, number of bootstraps (from marginal)
        X     % n x nBoot, Conditioned variable on Standard Margins
        Y     % n x nBoot, conditioning variable on Standard Margins
        A     % n x 1    %bin allocation
        RV    % Conditional Return Value Y | X
        n     % 1 x 1, number of observations
        nDmn  % 1 x 1, number of dimensions
        SmpLclRsdOn=true;   % 1 x 1, flag for sampling residuals locally from cvr bin (1) or globally (0)
        Sml   %Importance sampled simulation
        Alg= 'fit_fminsearch';  %algorithm used to fit HT model either fit_fminsearch or fit_newtonraphson
        PrmInd   %index of parameters, alpha=1, beta=2, mu=3, sig=4;
    end
    
    properties(Hidden=true)
        FigureFolder='Figures';  %folder to save figures into (default Figures)
        
    end
    
    properties(SetAccess=protected ) 
        nPrm     %total number of parameter
       
        Prm_A    %n x 4 link between parameter
        nBin     %number of model nBins
        nRtr     %number of return periods
        NonStat=false(4,1);   %4 x 1 flag for non stationary parameters
        
        RsdInd   % nBoot x 1 (cell array), index of residuals from each bootstrap
        %% Cross Validation stuff TODO given options to change these defaults
        CVMth=0; %0 Only Cross Validate smoothness for original dataset (fast);
        %1 Cross Validate smoothness for every bootstrap resample (slow),
        nCV=10; %no. cross-validation groups
        nSmth=10;   %no. smoothnesses tried in CV
        SmthLB=-4;   %lower bound (log10)  for smmothness range
        SmthUB=4;   %upper bound (log10)  for smmothness range
        
        SmthSet     %Set of Candidate smoothness params in Cross Validation
        OptSmth     %Optimal smoothness
        CVLackOfFit %Lack of Fit of candidate smoothness params in CV
        MarginType = 'Laplace'; %Laplace or Gumbel
        TwoParamFlag = false;%false;  %2P (alpha, beta) Hefferenan and Tawn model or 4P (alpha, beta , mu, sigma)
    end
    
    methods
        function HT=HeffernanTawn(Mrg,HTNEP,NonStationary,CV,SmpLclRsdOn)
            %HT=HeffernanTawn(Mrg,HTNEP,NonStationary,CV,SmpLclRsdOn)
            %INPUT
            % - Mrg 2 x 1, marignal model structure *output from stage 3
            % - HTNEP (scalar), Non-exceedance probability for conditional
            % - NonStationary flag for non stationary alpha
            %(Optional)
            %CV cross validation structure with control parameter
            
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
            HT.nDmn = numel(Mrg);
            
            %check have same number of observations for each margin
            for i = 2:HT.nDmn % loop over margins
                if size(Mrg(1).Y,1)~=size(Mrg(i).Y,1)
                    error('Marginal 1 and %d should have the same number of observations. Mrg1 nObs = %d; Mrg%d nObs = %d',i,size(Mrg(1).Y,1),i,size(Mrg(i).Y,1));
                end
            end
            %NEP must lead to postive threshold
            switch Mrg(1).MarginType
                case 'Gumbel' %Gumbel margin
                    validateattributes(HTNEP, {'numeric'},{'<=',1,'>=',exp(-exp(0))},'HeffernanTawn','HTNEP',2);  %0<=Tau<=1
                case 'Laplace'%Laplace margin
                    validateattributes(HTNEP, {'numeric'},{'<=',1,'>=',0.5},'HeffernanTawn','HTNEP',2);  %0<=Tau<=1
            end
            
            validateattributes(NonStationary,{'logical','numeric'},{'numel',4},'HeffernanTawn','NonStationary',3)
            HT.NonStat=logical(NonStationary(:));  %flag for nonstationary HT model
            if nargin>=5
                validateattributes(SmpLclRsdOn,{'logical','numeric'},{'integer','scalar'},'HeffernanTawn','SmpLclRsdOn',5)
                HT.SmpLclRsdOn=logical(SmpLclRsdOn); %flag for sampling residuals locally from bin or globally
            end
                      
            %check have same number of bootstrap resamples for each margin
            %HT.NonEmptyBins = []; %preallocate
            for i = 2:HT.nDmn % loop over margins
                %check have same number of bootstrap resamples for each margin
                if Mrg(1).nBoot~=Mrg(i).nBoot
                    error('Marginal 1 and %d should have the same number of bootstrap resamples. Mrg1 nBoot = %d; Mrg%d nBoot = %d',i,Mrg(1).nBoot,i,Mrg(i).nBoot);
                end
                %Check fitted each margin to same bootstrap resamples in
                if any(Mrg(1).BSInd(:)~=Mrg(i).BSInd(:))
                    error('Marginal 1 and %d should have the same bootstrap resamples',i);
                end
                %       HT.NonEmptyBins=unique([HT.NonEmptyBins;Mrg(i).NonEmptyBins;Mrg(j).NonEmptyBins]);
            end
            
            %% Smoothness parameters
            if nargin>=4  && ~isempty(CV) %% Cross Validation parameters (used for generalied Pareto fitting
                if isfield(CV,'CVMth')
                    validateattributes(CV.CVMth, {'numeric'},{'binary'},'MarginalModel','CV.CVMth',7);
                    HT.CVMth=CV.CVMth;
                end
                if isfield(CV,'nCV')
                    validateattributes(CV.nCV, {'numeric'},{'scalar','positive'},'MarginalModel','CV.nCV',7);
                    HT.nCV=CV.nCV;      %number cross-validation groups
                end
                if isfield(CV,'nSmth')
                    validateattributes(CV.nSmth, {'numeric'},{'scalar','positive'},'MarginalModel','CV.nSmth',7);
                    HT.nSmth=CV.nSmth;    %number smoothnesses tried in CV                                        
                end
                if isfield(CV,'SmthLB')
                    validateattributes(CV.SmthLB, {'numeric'},{'scalar'},'MarginalModel','CV.SmthLB',7);
                    HT.SmthLB=CV.SmthLB;   %lower bound (log10)  for smmothness range
                end
                if isfield(CV,'SmthUB')
                    validateattributes(CV.SmthUB, {'numeric'},{'scalar'},'MarginalModel','CV.SmthUB',7);
                    HT.SmthUB=CV.SmthUB;   %upper bound (log10)  for smmothness range
                end
            end
            
            HT.FigureFolder=Mrg(1).FigureFolder;
            
            %% Pre-allocation
            HT.n=size(Mrg(1).BSInd,1);
            HT.nBoot=Mrg(1).nBoot;
            
            HT.NEP = rand(HT.nBoot,1).*(HTNEP(2)-HTNEP(1))+HTNEP(1);
            HT.nBin=Mrg(1).Bn.nBin;
            
            Prm={'Alp','Bet','Mu','Sig'};
            
            if HT.TwoParamFlag  
               HT.NonStat(3:4)=false; 
            end
                
            HT.PrmInd=[]; 
           
            if HT.nBin == 1 %all stationary since only 1 bin
                HT.NonStat=false(4,1);
                HT.nPrm=ones(4,1);
                HT.PrmInd=(1:4)';
            else
                for iP=1:4 %loop oer parameters
                    if HT.NonStat(iP) %get number of model bins (number of alpha parameters).
                        HT.(Prm{iP})=NaN(HT.nBin,HT.nDmn-1,HT.nBoot);
                        HT.nPrm(iP)=HT.nBin;
                        HT.PrmInd=[HT.PrmInd;ones(HT.nBin,1).*iP];                       
                    else  %stationary model case
                        HT.(Prm{iP})=NaN(1,HT.nDmn-1,HT.nBoot);
                        HT.nPrm(iP)=1;
                        HT.PrmInd=[HT.PrmInd;iP];                        
                    end                                                                                    
                end
            end                 
            HT.Thr=NaN(HT.nBoot,HT.nDmn-1);
            HT.Rsd=cell(HT.nBoot,1); %different numbers of residuals for each bootstrap so store in cell
            HT.RsdInd=cell(HT.nBoot,1); %bin index of exceedences used in local sampling of residual;
            
            HT.Y=NaN(HT.n,HT.nDmn-1,HT.nBoot);
            HT.A=NaN(HT.n,HT.nBoot);
            HT.Prm_A=ones(HT.n,HT.nBoot,4);
            HT.OptSmth=NaN(HT.nBoot,1);
            
            if any(HT.NonStat)
                HT.SmthSet=logspace(HT.SmthLB,HT.SmthUB,HT.nSmth); %try range smoothness penalties for sigma varying by bin
            else
                HT.nSmth=1;
                HT.SmthSet=0;    %Switch off smoothness
            end
            HT.CVLackOfFit=NaN(HT.nSmth,HT.nBoot);
                                             
        end %HeffernanTawn constructor
        
        function HT=Fit(HT,Mrg,nSml)
            %Fit hefferenan and tawn model and compute return values
            %
            %INPUT
            %HT model object
            %Marginal model object
            %nSml number of simulation of importance sampler
            %
            %OUPTUT
            %HT fitted HT model with return values
            
            
            if nargin<3
                nSml=1e6; %number of importance smaples for diagnostic and contour plots
            else
                validateattributes(nSml,{'numeric'},{'integer','positive','scalar'},'HeffernanTawn','nSml',6)
            end
            
            HT.X=Margins(Mrg(1));
            
            %% Fit H&T Model
            for iBt=1:HT.nBoot           %loop over bootstrap resamples
                fprintf('Fitting for bootstrap sample %d of %d\n',iBt,HT.nBoot);
                %transform conditioned variable to Standard margins for iBt'th bootstrap
                HT.Thr(iBt)=Mrg(1).INV_Standard(HT.NEP(iBt));
                %HT.Thr(iBt)=quantile( HT.X(:,iBt),HT.NEP(iBt)); %find the H&T threshold
                IExc= HT.X(:,iBt)>HT.Thr(iBt);   %threshold exceedences
                %transform conditioning variable to Standard margins for iBt'th bootstrap
                
                INotNaN=Mrg(1).BSInd(:,iBt)>0;
                J=Mrg(1).BSInd(INotNaN,iBt);
                HT.A(INotNaN,iBt)=Mrg(1).Bn.A(J);    
                for iP=1:4
                    if HT.NonStat(iP)
                        HT.Prm_A(INotNaN,iBt,iP)=Mrg(1).Bn.A(J);                                       
                    end
                end
                        
                
                HT.RsdInd{iBt}=HT.A(IExc,iBt);
                for iDmn=2:HT.nDmn
                    HT.Y(:,iDmn-1,iBt)=Margins(Mrg(iDmn),iBt);
                end
                %% Fit Model
                HT=Fit_private(HT,IExc,iBt);
            end
            
            %% Compute conditional Return Value
            fprintf('Computing Conditional return value\n')
            HT=ConditionalReturnValue(HT,Mrg);
            
            %% simulate
            fprintf('Simulating under Heffernan and Tawn model using importance sampling\n')            
            HT.Sml=HT.SimulateIS(Mrg,nSml);   
        end
        
        function Sml=Simulate(HT,Mrg,nRls)
            %Sml=Simulate(HT,Mrg,nRls)
            %simulate realisations under HT model
            %INPUTS:
            % - nDmn x 1 Mrg, MarginalModel class Mrg (used to transform back to orginal scale)
            % - 1 x 1 nRls, number of realisations under the model
            %OUTPUT:
            % - Data structure Sml, data simulated from H&T model on Org,
            % Unif and Standard margins
            
            I=randi(HT.nBoot,nRls,1);%decide which bootstrap sample to use over all bootstraps HT.nBoot
            
            %% simulate covariate
            
            if Mrg(1).Bn.nBin>1
                Rat=Mrg(1).Rat(:,I);  %get rate of all observations
                RatCdf=cumsum(Rat)/sum(Rat); %get rate cdf
                
                %Simulate covariate with right rate
                Sml.A=sum(bsxfun(@gt,rand(1,nRls),RatCdf),1)'+1;  %bin allocation
                Sml.X=SampleCovariateFromBin(Mrg(1).Bn,Sml.A);  %Sample covariate value uniformly from  within bin
            else
                Sml.A=ones(nRls,1);
                Sml.X=rand(nRls,1)*360; %todo think about non periodic case here.
            end
                                         
            %% simulate on standard scale            
            U= rand(nRls,1); %sample uniform value
            G=INV_Standard(Mrg(1),U);
            
            Sml.StnMrg=NaN(nRls,HT.nDmn);    
            
            APrm=Sml.A.*HT.NonStat';
            APrm(APrm==0)=1; %nRls x 4
            
            for iRls=1:nRls  %loop over realisations
                
                if mod(iRls,1e3)==0
                    fprintf('HT Simulation %d of %d\n',iRls,nRls) %progress printer
                end
                
                iBt=I(iRls); %bootstrap index for curent realisation
                if G(iRls)>HT.Thr(iBt)  %X value above the threshold
                    tA=Sml.A(iRls);     
                    Z=HT.Rsd{iBt}(HT.RsdInd{iBt}==tA,:);  %nRsd x nDmn-1 joint residuals in each dimension
                    if size(Z,1)>0 && HT.SmpLclRsdOn %sample resid from bin if at least 1                     
                        ISmp= randi(size(Z,1),1);  %sample residual
                    else %too many bins -
                        Z=HT.Rsd{iBt};  %nRsd x nDmn-1 joint residuals in each dimension
                        ISmp= randi(size(Z,1),1);  %sample residual
                    end
                    Z=Z(ISmp,:);
                    
                    Sml.StnMrg(iRls,1)=G(iRls); %use random gumbel/laplace
                    %get conditioning value from H&T model
                    
                    Sml.StnMrg(iRls,2:end)=HT.Alp(APrm(iRls,1),:,iBt).* G(iRls) + G(iRls).^HT.Bet(APrm(iRls,2),:,iBt).* ( HT.Mu(APrm(iRls,3),:,iBt)+ HT.Sig(APrm(iRls,4),:,iBt).*Z);
                    
                else %resample below threshold
                    
                    IBlw=find(HT.X(:,iBt)<=HT.Thr(iBt)); %threshold on the gumbel scale
                    if any(IBlw)
                        if numel(IBlw)==1
                            ISmp=IBlw; %TODO check in nDmn
                        else
                            ISmp=randsample(IBlw,1); %TODO check in nDmn
                        end
                        Sml.StnMrg(iRls,1)=HT.X(ISmp,iBt);
                        Sml.StnMrg(iRls,2:end)=HT.Y(ISmp,:,iBt);
                    else
                        Sml.StnMrg(iRls,1)=HT.Thr(iBt);
                        Sml.StnMrg(iRls,2:end)=-Inf;
                    end
                end
                
            end
            
            %% Transform to uniform
            Sml.Unf=CDF_Standard(Mrg(1),Sml.StnMrg);
            
            %% Transform back to original margin
            Sml.Org=NaN(size(Sml.Unf));
            for iDmn=1:HT.nDmn %loop over dimension
                Sml.Org(:,iDmn)=Mrg(iDmn).INV(Sml.Unf(:,iDmn),I,Sml.A);
            end
            
        end %simulate
        
        function Sml=SimulateIS(HT,Mrg,nRls)
            %Sml=SimulateIS(HT,Mrg,nRls)
            %simulate uniformly space realisations under HT model (used in importance sampling
            %contours)
            %INPUTS:
            % - nDmn x 1 Mrg, MarginalModel class Mrg (used to transform back to orginal scale)
            % - 1 x 1 nRls, number of realisations under the model
            %OUTPUT:
            % - Data structure Sml, data simulated from H&T model on Org,
            % Unif and Standard margins
            
            Sml.I=randi(HT.nBoot,nRls,1);%decide which bootstrap sample to use over all bootstraps HT.nBoot
          
            %% simulate covariate
            RatNrm=Mrg(1).Rat./sum(Mrg(1).Rat,1); %nBin x nBoot
            %Simulate covariate with right rate
            Sml.A=randi(HT.nBin,nRls,1);  %bin allocation
            Sml.X=SampleCovariateFromBin(Mrg(1).Bn,Sml.A);  %Sample covariate value uniformly from  within bin
            
            if HT.nBin==1
                J=Sml.I;
                Sml.W=RatNrm(J)'; %nRls x 1 %probability weight for chosen bin;
            else
                J=sub2ind([HT.nBin,HT.nBoot],Sml.A,Sml.I);  %common index across bin allocation and bootstrap
                Sml.W=RatNrm(J); %nRls x 1 %probability weight for chosen bin;
            end
                        
            %% simulate on original scale
            %gp upper endpoint in each dimension
            
            Sml.Org= rand(nRls,HT.nDmn); %simulation resluts on original scale
            Sml.Unf= NaN(nRls,HT.nDmn);   %simulation resluts on uniform scale
            Sml.StnMrg= NaN(nRls,HT.nDmn); %simulation resluts on standard margins            
            logf_Stn=NaN(nRls,HT.nDmn);  %log importance weights on standard margins
            logf=NaN(nRls,HT.nDmn);  %%log importance weights on original margins
            logg_Stn=NaN(nRls,HT.nDmn);  %log importance weights on standard margins
            logg=NaN(nRls,HT.nDmn);  %%log importance weights on original margins                       
            
            %% simulated data
            for iDmn=1:HT.nDmn
                %% Find sensible range to draw values over
                [UL,LL,Rng] = Mrg(iDmn).makeRange(Sml.I,Sml.A);
                Sml.Org(:,iDmn)= Sml.Org(:,iDmn).*Rng+LL; %sample uniform value
                
                logg(:,iDmn)=-log(Rng); %importance weights                    
                
                %transform to uniform scale;
                Sml.Unf(:,iDmn)=Mrg(iDmn).CDF(Sml.Org(:,iDmn),Sml.A,Sml.I);
                %transform to standard margins;
                Sml.StnMrg(:,iDmn)=INV_Standard(Mrg(iDmn),Sml.Unf(:,iDmn));
                                
                logfx = Mrg(iDmn).LogPDF_Standard(Sml.StnMrg(:,iDmn)); %density of values on standard margins
                
                %% Compute Jacobians for transformation to orignal marigns
                J_Unf_Stn=logfx; %Uniform --> Gumbel transform Jacobian
                J_Unf_GamGP=Mrg(iDmn).LogPDF(Sml.Org(:,iDmn),Sml.A,Sml.I); %Uniform --> GamGP transform Jacobian
                
                logg_Stn(:,iDmn)=logg(:,iDmn)-J_Unf_GamGP+J_Unf_Stn;
                
                %% Density of computed points
                if iDmn==1  %first dimension standard marginal denisty
                    %compute density of chosen points (other dimensions use the HT density)
                    logf_Stn(:,1)=logfx;
                    
                    IExc=HT.Thr(Sml.I)<Sml.StnMrg(:,1); %index of exceedences
                    nE=sum(IExc); %number of exceedences
                    logfx1=logfx;
                    
                else %conditional density of each subsequent dimension.
                    iAsc=iDmn-1;                                                           
                    DimSet=[1,iDmn];
                    
                    logfyx=NaN(nRls,1);
                    %% below threshold 2D kernel density estimation
                    Q=[HT.X(:,1),HT.Y(:,iAsc,1)]; %don't use bootrap uncertainty
                    Z=Sml.StnMrg(:,DimSet);
                    logfyx(~IExc) = log(ksdensity(Q,Z(~IExc,:))); %joint empirical density of data.
                    logfyx(~IExc)= logfyx(~IExc)-logfx1(~IExc); %conditional density                    
                    %% above threshold  HT density
                    %HT parameters;
                    IBtExc=Sml.I(IExc); %bootstrap index of exceedences only                                     
                    
                    % Y= aX+X^b (m+s Z);
                    %  ((Y- aX) -m X ^b)./ s X^b =Z
                     %% HT parameter
                    tPrm=NaN(nRls,4);
                    PrmLbl={'Alp','Bet','Mu','Sig'};
                    for iP=1:4
                        if HT.NonStat(iP) %non stationary alpha
                            t1=permute(HT.(PrmLbl{iP})(:,iAsc,:),[1,3,2]);  %nBin x nBoot
                            tPrm(:,iP)=t1(J); %nRls x 1 
                        %%
                        else
                            tPrm(:,iP)=HT.(PrmLbl{iP})(1,iAsc,Sml.I); % nRls x nDmn-1
                        end
                    end
                    tPrm=tPrm(IExc,:);
                    
                    Xpb=Z(IExc,1).^tPrm(:,2);
                    mu=Z(IExc,2)-tPrm(:,1).*Z(IExc,1)-Xpb.*tPrm(:,3);
                    sgm=Xpb.*tPrm(:,4);
                    ZStn=mu./sgm;
                    %residual density comes from kernel deinsty estimation
                    if 1
                        fygx_HT=NaN(nE,1);
                        
                        for iBt=1:HT.nBoot %loop over bootstraps
                            if HT.SmpLclRsdOn
                                for iA=1:HT.nBin
                                    tRsd=HT.Rsd{iBt}(HT.RsdInd{iBt}==iA,iAsc);  %nRsd x nDmn-1 joint residuals in each dimension
                                    if size(ZStn,1)>0 %sample resid from bin if at least 1
                                        ISmp= randi(size(tRsd,1),1);  %sample residual
                                    else
                                        tRsd=HT.Rsd{iBt}(:,iAsc);  %nRsd x nDmn-1 joint residuals in each dimension
                                        ISmp= randi(size(tRsd,1),1);  %sample residual
                                    end
                                    tRsd=tRsd(ISmp,:);
                                    I=IBtExc==iBt & Sml.A(IExc)==iA; %
                                    if any(I)
                                        fygx_HT(I)=ksdensity(tRsd,ZStn(I)); %f(y|x) on standard scale
                                    end
                                end
                            else
                                I=IBtExc==iBt; %
                                if any(I)
                                    fygx_HT(I)=ksdensity(HT.Rsd{iBt}(:,iAsc),ZStn(I)); %f(y|x) on standard scale
                                end
                            end
                        end
                        %f(x,y) = f(y|x)f(x)
                        logfyx(IExc)=log(fygx_HT); %conditional density f(y|x) on laplace scale
                        
                    else
                        % use gaussian density for residuals
                        fygx_HT=-log(2*pi)/2-log(sgm)-0.5*ZStn.^2; 
                        fygx_HT(isnan(fygx_HT))=-Inf;
                        logfyx(IExc)=fygx_HT;
                    end
                   
                    logf_Stn(:,iDmn)=logfyx+log(Rng); %importance weights
                                       
                end
                               
                logf(:,iDmn)=logf_Stn(:,iDmn)+J_Unf_GamGP-J_Unf_Stn; 
                               
            end %loop over dimensions
            
            Sml.logfog=logf-logg;
            Sml.logfog_Stn=logf_Stn-logg_Stn;
            
            %Handle bad cases set density to zero
            if any(isnan(Sml.logfog(:)))
                Sml.logfog(isnan(Sml.logfog))=-Inf;
            end
            if any(isnan(Sml.logfog_Stn(:)))
                Sml.logfog_Stn(isnan(Sml.logfog_Stn))=-Inf;
            end                         
            
        end %SimulateIS
        
        function HT=ConditionalReturnValue(HT,Mrg,nRls)
            %HT=ConditionalReturnValue(HT)
            %compute conditional return value Y | X using monte carlo simuation use return period defined in marginal model.
            if nargin<=2
                nRls=1000; %default number of MC draws
            end
            
            HT.nRtr=numel(Mrg(1).RtrPrd);
            %preallocate
            nAsc=HT.nDmn-1;
            RV.X_Stn=NaN(HT.nBin,nRls,HT.nRtr);
            RV.X=NaN(HT.nBin,nRls,HT.nRtr);
            RV.Y_Stn=NaN(HT.nBin,nAsc,nRls,HT.nRtr);
            RV.Y=NaN(HT.nBin,nAsc,nRls,HT.nRtr);
            
            I=randi(HT.nBoot,nRls,1);  %bootstrap samples to use
            RV.I=I; %store bootstraps for use in threshold plot
            %draw random resdiuals
            
            if ~HT.SmpLclRsdOn   %if too many bins (few obs per bin), sample residuals globally
                Z=cell2mat(cellfun(@(x)x(randi(numel(x),HT.nBin,1))',HT.Rsd(I,:),'uniformoutput',false))'; %TODO check multivariate cases
                Z=permute(Z,[1,3,2]);
                
            else   %when your bins are big enough (decent no. of obs per bin), sample residuals locally from bin
                Z=NaN(nRls,nAsc,HT.nBin);
                for iBin=1:HT.nBin
                    tRsd=cellfun(@(x,y)y(x==iBin,:),HT.RsdInd(I,:),HT.Rsd(I,:),'uniformoutput',false);
                    J=cellfun(@length,tRsd)>0;
                    Z(J,:,iBin)=cell2mat(cellfun(@(x)x(randi(size(x,1)),:),tRsd(J,:),'uniformoutput',false));
                end
                Z=permute(Z,[3,2,1]);
            end
            
            for iRtr=1:HT.nRtr %loop over return periods
                
                %% Sample from return value distribution within each bin
                rho=Mrg(1).Rat(:,I); %annual rate of occurence
                LT=rho*Mrg(1).RtrPrd(iRtr); %poisson Rate
                
                UX=rand(HT.nBin,nRls); %U should be in the range [ P0, 1] where P0 is the non occurence rate.
                P=1+log(UX)./(LT);    %adjust for return period  (inv of exp(-L*T*(1-C));
                P(bsxfun(@lt,P,HT.NEP(I)'))=NaN;  %this is the non-exceedence rate on standard scale
                %P(P<0)=NaN;  %this is the non-exceedence rate on standard scale
                RV.X(:,:,iRtr)=Mrg(1).INV(P,I); %RVX values on original scale
                
                %transform form uniform to standard margins using CDF
                RV.X_Stn(:,:,iRtr)=INV_Standard(Mrg(1),P);
                
                %compute Y using conditional model given X
                tX=permute(RV.X_Stn(:,:,iRtr),[1,3,2]);
                RV.Y_Stn(:,:,:,iRtr)=bsxfun(@times,HT.Alp(:,:,I),tX) + bsxfun(@power,tX,HT.Bet(:,:,I)).*bsxfun(@plus,HT.Mu(:,:,I),bsxfun(@times,HT.Sig(:,:,I),Z));
                
                for iAsc=1:nAsc
                    %transform Y from standard to uniform margins using CDF
                    UY=permute(Mrg(iAsc+1).CDF_Standard(RV.Y_Stn(:,iAsc,:,iRtr)),[1,3,2]); %nBin x nAsc x nRls x nRtr
                    %Transform Y to back original margins
                    RV.Y(:,iAsc,:,iRtr)=permute(Mrg(iAsc+1).INV(UY,I),[1,3,2]);
                end
            end
            
            %% Get Omni Value (covariate free) return value and its associated conditions
            if HT.nBin > 1
                
                XOmni_Stn=max(RV.X_Stn,[],1);
                [XOmni,J]=max(RV.X,[],1);    %J is index of location of max
                
                YOmni=NaN(1,nAsc,nRls,HT.nRtr);
                YOmni_Stn=NaN(1,nAsc,nRls,HT.nRtr);
                for iRtr=1:HT.nRtr
                    I=sub2ind([HT.nBin,nRls],J(:,:,iRtr),(1:nRls)); %need composite index of realisation and max location
                    %original margins
                    tYOmni=reshape(permute(RV.Y(:,:,:,iRtr),[1,3,2]),HT.nBin*nRls,nAsc);
                    YOmni(1,:,:,iRtr)=permute(tYOmni(I,:),[3,2,1]);
                    %standard margins
                    tYOmni_Stn=reshape(permute(RV.Y_Stn(:,:,:,iRtr),[1,3,2]),HT.nBin*nRls,nAsc);
                    YOmni_Stn(1,:,:,iRtr)=permute(tYOmni_Stn(I,:),[3,2,1]);
                end
                
                RV.Y=cat(1,RV.Y,YOmni);
                RV.Y_Stn=cat(1,RV.Y_Stn,YOmni_Stn);
                RV.X=cat(1,RV.X,XOmni);
                RV.X_Stn=cat(1,RV.X_Stn,XOmni_Stn);
                
            end
            
            %% Store Return Value simulation
            HT.RV=RV;
            HT.RV.nRls=nRls;
            
        end %ConditionalReturnValue
        
        function Plot(HT,Mrg)
            %% Store bin start and end points for plot labels
            %             BinSt=Mrg(1).DrcEdg;
            %             BinEnd=circshift(Mrg(1).DrcEdg,-1);
            %             if Mrg(1).DrcEdg(1) > 0 %if there is a bin which straddles 0, make it bin 1
            %                 BinSt=circshift(BinSt,1);
            %                 BinEnd=circshift(BinEnd,1);
            %             end
            %% simulate under model
            nRls=length(Mrg(1).Y)*20;
            MCSml=Simulate(HT,Mrg,nRls); %simulate 10 times length of original data
            
            %% Plot  data and simulation
            figure;
            clf;
            nG=50;
            Edg=linspace(min(HT.Sml.Org(:,1)),max(Mrg(1).RVMed(:).*1.3),nG+1);
            EdgStn=linspace(Mrg(1).INV_Standard(0.01),max(HT.Sml.StnMrg(:,1)),nG+1);            
            Md=(Edg(2:end)+Edg(1:end-1))./2;
            MdStn=(EdgStn(2:end)+EdgStn(1:end-1))./2;
            Dst=(Md-HT.Sml.Org(:,1))./(Edg(end)-Edg(1));
            Dst_Stn=(MdStn-HT.Sml.StnMrg(:,1))./(EdgStn(end)-EdgStn(1));
            nu=0.01;
            
            for iDmn = 2:HT.nDmn
                
                %% compute quantiles of regression line from importance sampling
                [Q,QStn]=deal(NaN(nG,3));                
                
                %Original margins
                [tY,srtI]=sort(HT.Sml.Org(:,iDmn));                                              
                tW=exp(sum(HT.Sml.logfog(srtI,:),2)+log(HT.Sml.W)-Dst(srtI,:).^2./(nu.^2));
                tP=cumsum(tW,1)./sum(tW,1);
                
                %Standard margins
                [tY_Stn,srtI]=sort(HT.Sml.StnMrg(:,iDmn));               
                tW_Stn=exp(sum(HT.Sml.logfog_Stn(srtI,:),2)+log(HT.Sml.W)-Dst_Stn(srtI,:).^2./(nu.^2));
                tP_Stn=cumsum(tW_Stn,1)./sum(tW_Stn,1);
 
                for iA=1:nG
                    %% Original Margins
                    [sP,J]=unique(tP(:,iA));                    
                    if numel(sP)>1
                        Q(iA,:)=interp1(sP,tY(J),[0.1,0.5,0.9]);
                    end
                    %% Original Margins
                    [sP,J]=unique(tP_Stn(:,iA));
                    if numel(sP)>1
                        QStn(iA,:)=interp1(sP,tY_Stn(J),[0.1,0.5,0.9]);
                    end
                end
                
                subplot(2,HT.nDmn-1,iDmn-1)                
                plot(MCSml.StnMrg(:,1),MCSml.StnMrg(:,iDmn),'.','color',[1,1,1]*0.8)
                hold on
                plot(HT.X(:,1),HT.Y(:,iDmn-1,1),'k.','markersize',10)
                hold on
                plot(MdStn,QStn,'r--','linewidth',2);
              
                xlabel(sprintf('%s: Conditioning variable',Mrg(1).RspLbl))
                ylabel(sprintf('%s: Conditioned variable',Mrg(iDmn).RspLbl));
                title(sprintf('%s|%s: H&T simulation on %s margins',Mrg(iDmn).RspLbl,Mrg(1).RspLbl,Mrg(1).MarginType))
                legend('Simulation','Data','location','NorthWest')
                axis tight
                box on
                grid on
                              
                ylim([Mrg(2).INV_Standard(0.01),max(HT.Sml.StnMrg(:,2))])
                xlim([Mrg(1).INV_Standard(0.01),max(HT.Sml.StnMrg(:,1))])
                hold on
                plot([min(HT.Thr),max(HT.Thr)].*[1,1]',ylim,'k--')
                plot(xlim,xlim,'k--');
                
                subplot(2,HT.nDmn-1,HT.nDmn-1+iDmn-1)
                plot(MCSml.Org(:,1),MCSml.Org(:,iDmn),'.','color',[1,1,1]*0.8)
                hold on
                plot(Mrg(1).Y,Mrg(iDmn).Y,'k.','markersize',10)                            
                hold on
                plot(Md,Q,'r--','linewidth',2);
                hold on
                plot([1,1]'.*Mrg(1).RVMed(end,:),ylim,'g--','linewidth',2)
                xlabel(sprintf('%s: Conditioning variable',Mrg(1).RspLbl))
                ylabel(sprintf('%s: Conditioned variable',Mrg(iDmn).RspLbl));
                title(sprintf('%s|%s: H&T simulation on original scale',Mrg(iDmn).RspLbl,Mrg(1).RspLbl))
                axis tight
                grid on
            end
            savePics(fullfile(HT.FigureFolder,'Stg4_HT_1_SmlvsData'))
            
            %% Residual Diagnostics orignal sample
            figure;
            clf;
            tRsd=cell2mat(HT.Rsd);
            
            %histogram of redisuals in ith variable
            for iDmn = 1:(HT.nDmn-1)
                subplot((HT.nDmn-1),2,2*iDmn-1)
                if verLessThan('Matlab','8.5')
                    hist(tRsd(:,iDmn));
                    h = findobj(gca,'Type','patch');
                    set(h,'FaceColor',[1 1 1]*0.5);
                    set(h,'linestyle','none')
                else
                    histogram(tRsd(:,iDmn),'edgecolor','none','facecolor','k')
                end
                axis tight
                title(sprintf('%s|%s: resid hist',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl))
                xlabel('Residual')
            end
            %qq plot of redisuals in ith variable
            for iDmn = 1:(HT.nDmn-1)
                
                subplot(HT.nDmn-1,2,2*iDmn);
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
                title(sprintf('%s|%s: resid QQ-plot',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl))
            end
            savePics(fullfile(HT.FigureFolder,'Stg4_HT_2_ResidualDiagnostics1'))
            
            
            %% Residual diagnostics: residuals against covariates
            figure;clf;
            iP = 0; %counter on subplot number
            for iDmn = 1:(HT.nDmn-1)
                for iC=1:Mrg(1).nCvr
                    iP = iP +1;
                    subplot(HT.nDmn-1,Mrg(1).nCvr,iP)
                    IExc=HT.X(:,1)>HT.Thr(1);   %threshold exceedences
                    IExc=IExc(Mrg(1).BSInd(:,1)>0);
                    plot(Mrg(1).X(IExc,iC),HT.Rsd{1}(:,iDmn),'k.')
                    axis tight
                    box on
                    grid on
                    xlabel(Mrg(1).CvrLbl(iC))
                    title(sprintf('%s|%s: Residuals on %s',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl,Mrg(1).CvrLbl{iC}))
                end
            end
            
            
            savePics(fullfile(HT.FigureFolder,'Stg4_HT_2_ResidualDiagnostics2'))
            
            %% Parameter Plot (bootstrap uncertainty)
            LgnLbl=cell(HT.nDmn-1,1);
            for iDmn = 1:(HT.nDmn-1)
                LgnLbl{iDmn} = sprintf('%s|%s',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl);
            end
            Lbl={'\alpha','\beta','m','s'};
            
            for iDmn = 1:(HT.nDmn-1) %different parameters figure for each dimension
                figure;
                clf;
                
                PrmLbl={'Alp','Bet','Mu','Sig'};
                for iP=1:4          %loop over parameters                  
                    if HT.NonStat(iP)
                        for iC=1:Mrg(1).nCvr
                            subplot(4,Mrg(1).nCvr,(iC-1)*4+iP)
                            hold on                            
                            PlotParameter(Mrg(1).Bn,squeeze(HT.(PrmLbl{iP})(:,iDmn,:)),iC,'color',[0,0,0],'linewidth',2)                            
                            xlabel(Mrg(1).CvrLbl(iC))
                            ylabel(Lbl{iP})
                            box on
                            grid on
                        end
                    else
                        subplot(4,1,iP)                        
                        Prm=squeeze(HT.(PrmLbl{iP})(1,iDmn,:));                        
                                               
                        histogram(Prm,'edgecolor','none','facecolor','k')
                    
                        hold on
                        xlabel(Lbl{iP})                                                                    
                    end
                    title(LgnLbl{iDmn})
                end
                axes('position',[0.1300    0.1100    0.7750    0.8150]);
                axis off
                title(sprintf('%s|%s: Histograms of H&T parameter uncertainty',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl))
                
                savePics(fullfile(HT.FigureFolder,sprintf('Stg4_HT_3_Parameters_%s',Mrg(iDmn+1).RspLbl)))
            end
            
            %% Threshold plot  (alpha as function of NEP)
            %Alpha varies by covariate bin, so to avoid over-crowded plots; have a different version of this plot for each dimension.
            nPlt1=ceil(sqrt(HT.nPrm(1))); %max size nPlt x nPlt
            nPlt2=ceil(HT.nPrm(1)./nPlt1);                        
            
            for iDmn = 1:(HT.nDmn-1)
                figure;
                clf;
                for iAlp = 1:HT.nPrm(1)
                    subplot(nPlt2,nPlt1,iAlp)
                    tAlp=squeeze(HT.Alp(iAlp,iDmn,:));
                    plot(HT.NEP,tAlp,'k.','markersize',20)
                                        
                    if Mrg(1).nBoot>=20
                        hold on
                        nNEP=10; %number of bins for threshold diagnostic plot
                        NEPEdg=linspace(min(HT.NEP),max(HT.NEP),nNEP+1);
                        NEPBin=(NEPEdg(1:end-1)+NEPEdg(2:end))/2;
                        [~,tA] = histc(HT.NEP,NEPEdg);
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
                        xlabel('HT NEP')
                        ylabel('\alpha')
                        if HT.nPrm(1)==1 %stationary case
                            title(sprintf('%s|%s: H&T parameter stability by threshold',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl))
                        else
                            title(sprintf('%s|%s: H&T parameter stability by threshold \n Bin %s',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl,Mrg(1).Bn.BinLbl{iAlp}))
                        end
                    else
                        title(sprintf('Bin %s',Mrg(1).Bn.BinLbl{iAlp}))
                    end
                end
                hold off;
                savePics(fullfile(HT.FigureFolder,sprintf('Stg4_HT_4_AlphaThresholdStability_%s',Mrg(iDmn+1).RspLbl)))
            end
            
                                   
            %% Lack of fit plot
            if any(HT.NonStat)
                figure;
                clf;
                if HT.nBoot>1 && HT.CVMth==1
                    plot(log10(HT.SmthSet),nanmedian(HT.CVLackOfFit,2),'k-','linewidth',2)
                    hold on
                    plot(log10(HT.SmthSet),quantile(HT.CVLackOfFit,0.025,2),'k--','linewidth',2)
                    plot(log10(HT.SmthSet),quantile(HT.CVLackOfFit,0.975,2),'k--','linewidth',2)
                else
                    plot(log10(HT.SmthSet),HT.CVLackOfFit(:,1),'k-','linewidth',2)
                end
                axis tight
                hold on
                grid on
                plot(log10(median(HT.OptSmth))*[1,1],ylim,'r--','linewidth',2)
                ylabel('Lack of fit')
                xlabel('$\log_{10}(\tilde{\lambda)}$','Interpreter','latex')
                title(sprintf('HT cross-validation lack of fit'))
                savePics(fullfile(HT.FigureFolder,'Stg4_HT_5_SectorGoodnessOfFit'))
            end
            
            %% Cond. RV CDFs
            ColMat=hsv(HT.nBin);
            
            figure;
            clf;
            c=0;
            for iDmn=2:HT.nDmn %loop over associated variables
                for iRtr=1:HT.nRtr
                    c=c+1;
                    
                    subplot(HT.nDmn-1,HT.nRtr,c)
                    grid on
                    % Per bin: -----------
                    if HT.nBin<16  && HT.nBin > 1  %plot directional CDFs when there are no more than 15 bins
                        try %sort() different in matlab versions
                            tX=sort(permute(HT.RV.Y(1:HT.nBin,iDmn-1,:,iRtr),[1,3,2]),2,'MissingPlacement','first');
                        catch
                            tX=sort(permute(HT.RV.Y(1:HT.nBin,iDmn-1,:,iRtr),[1,3,2]),2);
                        end
                        for iBin = 1:HT.nBin
                            plot(tX(iBin,:),linspace(0,1,HT.RV.nRls),'linewidth',2,'color',ColMat(iBin,:))
                            hold on
                        end
                    end
                    
                    %Omni: --------------
                    try
                        tX=sort(permute(HT.RV.Y(end,iDmn-1,:,iRtr),[1,3,2]),2,'MissingPlacement','first');
                    catch
                        tX=sort(permute(HT.RV.Y(end,iDmn-1,:,iRtr),[1,3,2]),2);
                    end
                    if HT.nBin == 1
                        plot(tX,linspace(0,1,HT.RV.nRls),'k-','linewidth',2)
                    else
                        plot(tX,linspace(0,1,HT.RV.nRls),'--k','linewidth',2)
                    end
                    
                    ylabel('Cumulative Probability')
                    xlabel(sprintf('%s | %s',Mrg(iDmn).RspLbl, Mrg(1).RspLbl))
                    grid on
                    
                    if (HT.nBin<16  && HT.nBin >1 ) && (c == (HT.nDmn-1)*HT.nRtr)
                        legend([Mrg(1).Bn.BinLbl(:);'Omni'],'location','best');
                    end
                    
                    if iDmn==2
                        title(sprintf('Cnd RV CDF %g Years',Mrg(iDmn).RtrPrd(iRtr)))
                    end
                end
            end
            savePics(fullfile(HT.FigureFolder,'Stg4_HT_6_ConditionalReturnValueCDF'))
        end %plot
        
    end %methods
    methods% (Access = private)
        function HT=Fit_private(HT,IExc,iBt)
            % HT=Fit_private(HT,IExc,iBt)
            % Fit single bootstrap of HT model
            % INPUTS:
            % - n x 1 IExc logical array of exceedances
            % - 1 x 1 iBt bootstrap number
            % OUTPUTS:
            % - HT.Alp(:,:,iBt) fitted H&T parameters for iBt'th bootstrap
            % - HT.Bet(:,:,iBt) fitted H&T parameters for iBt'th bootstrap
            % - HT.Mu(:,:,iBt) fitted H&T parameters for iBt'th bootstrap
            % - HT.Sig(:,:,iBt) fitted H&T parameters for iBt'th bootstrap
            % - HT.Rsd{:,iBt} residuals for iBt'th bootstrap
                        
            %% Get exceedances
            nExc=sum(IExc);
            tX=HT.X(IExc,iBt); %exceedances for iBt'th bootstrap
            tY=HT.Y(IExc,:,iBt);  %conditioning value given exceedance in conditioned value            
            tA=permute(HT.Prm_A(IExc,iBt,:),[1,3,2]);  %nExc x 4
            if sum(isnan(tA))>0
                warning('NaNs in allocation')
            end
                                    
            %% Flag for do Cross Validation
            if (HT.CVMth==0 &&  iBt>1) || HT.nSmth==1
                %% Cross validation off
                if HT.nSmth==1  %only one smoothness
                    HT.OptSmth(iBt)=HT.SmthSet;
                else  %use first bootstrap instead of cross validating every time
                    HT.OptSmth(iBt)=HT.OptSmth(1);
                end
            else
                fprintf('Starting Cross Validation to estimate HT alpha smoothness:\n');
                LackOfFit=NaN(HT.nCV,HT.nSmth); %lack of fit
                ICV=randi(HT.nCV,nExc,1); %split nExc observations into nCV cross validation groups
                for iCV=1:HT.nCV  %cross validation loop
                    fprintf('.')
                    %Pull out fit set
                    tXFit=tX(ICV~=iCV); %main response
                    tYFit=tY(ICV~=iCV,:); %associated variable
                    tAFit=tA(ICV~=iCV,:);
                    %Pull out prediction set
                    tXPrd=tX(ICV==iCV);
                    tYPrd=tY(ICV==iCV,:);
                    tAPrd=tA(ICV==iCV,:);
                    for iLam=1:HT.nSmth %Smoothness penalties                        
                        PrmHat=HT.(HT.Alg)(tXFit,tYFit,tAFit,HT.SmthSet(iLam));                       
                        LackOfFit(iCV,iLam)=HT.likelihood(tXPrd,tYPrd,tAPrd,PrmHat);
                    end %SigPen
                end %CV
                fprintf('\n')
                HT.CVLackOfFit(:,iBt)=sum(LackOfFit,1)';
                [~,tI]=min(HT.CVLackOfFit(:,iBt));
                HT.OptSmth(iBt)=HT.SmthSet(tI)';
            end
            
            %% fit model: find optimal parameter values for heffernan and tawn model           
            tPrm=HT.(HT.Alg)(tX,tY,tA,HT.OptSmth(iBt));           
            %% keep residuals
            HT.Alp(:,:,iBt)=tPrm(HT.PrmInd==1,:);
            HT.Bet(:,:,iBt)=tPrm(HT.PrmInd==2,:);
            HT.Mu(:,:,iBt)=tPrm(HT.PrmInd==3,:);
            HT.Sig(:,:,iBt)=sqrt(tPrm(HT.PrmInd==4,:));
            tXpB=tX.^HT.Bet(tA(:,2),:,iBt);
            mn=tY-HT.Alp(tA(:,1),:,iBt).*tX-HT.Mu(tA(:,3),:,iBt).*tXpB;
            std= HT.Sig(tA(:,4),:,iBt).*tXpB;
            HT.Rsd{iBt}=mn./std;
            if any(~isreal(HT.Rsd{iBt}(:)))
                warning('Complex Residuals found!!\n')
            end
            
        end %fit
        
        function PrmHat=fit_fminsearch(HT,X,Y,A,Lam)
            %Fit Hefferenan and Tawn model using simple fminsearch algorithm, works well for small
            %number of bins
            %INPUT
            %X    %conditioning data
            %Y    %conditioned data
            %A    %bin allocation
            %Lam  %smoothness parameter
            
            %% Starting Solution
            p0=HT.startingSolution(X,Y,A);
            
            opts=optimset('display','off');
            
            PrmHat=fminsearch(@(p)HT.likelihood(X,Y,A,p,Lam),p0,opts);
            
            if HT.TwoParamFlag
                PrmHat=[PrmHat;zeros(HT.nDmn-1,1)';ones(HT.nDmn-1,1)'];
            end
        end
        
        function PrmHat=fit_fmincon(HT,X,Y,A,Lam)
            %Fit Hefferenan and Tawn model using simple fit_fmincon algorithm, works well for small
            %number of bins
            %INPUT
            %X    %conditioning data
            %Y    %conditioned data
            %A    %bin allocation
            %Lam  %smoothness parameter
            
            %% Starting Solution
            p0=HT.startingSolution(X,Y,A);
            
            %opts=optimoptions(@fmincon,'algorithm','interior-point','SpecifyObjectiveGradient',true);
            opts=optimoptions(@fmincon,'Display','notify','algorithm','interior-point','SpecifyObjectiveGradient',true,'CheckGradients',false);
            
             %Alp in [-1,1], Bet in [-Inf, 1], Mu[-Inf,Inf], Sig in [0,Inf];
            LB=[ones(HT.nAlp,1)*-1;-Inf;-Inf;0];  %
            UB=[ones(HT.nAlp,1)*1;1;Inf;Inf];
                        
            %assumes mu at the end!!
            A_con=[blkdiag([eye(HT.nAlp);eye(HT.nAlp)],1,-1),zeros(2*HT.nAlp+2,1) ];
            A_con=A_con(:,[1:end-2,end,end-1]); %swap mu and sig around to make sure in the right order!!
            %(2p+2) x (p+3)  (p+3) x 1    
            b=[ones(2.*HT.nAlp,1);1;0];
            %  (2p+2) x 1
            
            PrmHat=fmincon(@(p)HT.myobj(X,Y,A,p,Lam),p0,A_con,b,[],[],LB,UB,[],opts);
            
            if HT.TwoParamFlag
                PrmHat=[PrmHat;zeros(HT.nDmn-1,1)';ones(HT.nDmn-1,1)'];
            end
        end
        
        function p=fit_newtonraphson(HT,X,Y,A,Lam)
            %Fit Hefferenan and Tawn model using newton rapshon algorithm, works well for
            %bigger cases
            %INPUT
            %X    %conditioning data
            %Y    %conditioned data
            %A    %bin allocation
            %Lam  %smoothness parameter
            p=HT.startingSolution(X,Y,A);
            FitInd=[1,2,1,2];
            if HT.TwoParamFlag                
                FitInd=FitInd(HT.PrmInd(HT.PrmInd<=2));
            else
                FitInd=FitInd(HT.PrmInd);
            end
            
            maxIter=1000; %maximum number of iterartions of solver
            tol=1e-3; %tolerance for convergence
            Dlt0=[0.2,0.2]; %step size on [0,1])
            Dlt=Dlt0; %step size on [0,1])
            
            for iIt=1:maxIter %iteration loop
                Chg=0;
                % iterate between pairs (alpha,mu) & (beta,sigma)
                for iP=1:2
                    %G p x nDmn, H p x p nDmn
                    [G,H]=HT.Gradient(iP,X,Y,A,p,Lam);
                    
                    newp=p;
                    for iDmn=1:HT.nDmn-1
                        newp(FitInd==iP,iDmn)= p(FitInd==iP,iDmn)-Dlt(iP).*(H(:,:,iDmn)\G(:,iDmn));
                    end
                    
                    if HT.checkParam(newp)
                        %parameters are bad so reset and reduce step size
                        Dlt(iP)=Dlt(iP).*0.1; %
                        Chg=Inf; %need another loop
                    else %good case
                        Dlt(iP)=Dlt0(iP); %reset step size on good step                        
                        Chg=Chg+sum(abs((newp(:)-p(:))));
                        p=newp; %set new parameters
                    end
                    
                end %loop over parameters
                
                %fprintf('%g\n',Chg);
                %% Check parameters              
                if Chg<tol %stopping criteria
                 %   fprintf('Converged in %g iterations\n',iIt)
                    break
                end
            end %iterations
                                    
            if HT.TwoParamFlag
                p=[p;zeros(1,HT.nDmn-1);ones(1,HT.nDmn-1)];
            end                        
        end                        
        
        function p0=startingSolution(HT,X,Y,A)
            
            %% Assume b=0 for starting solution 
            b0=0;
            
            %% Alp and Mu | b
            if any(HT.NonStat([1,3,4]))
                B=A(:,find(HT.NonStat,1))==1:HT.nBin;
               
            end
            if HT.NonStat(1)
                BX=B.*X;                 
            else
                BX=X;                
            end
            if HT.NonStat(3)
                BXpb=B.*X.^b0;
               
            else
                BXpb=X.^b0;
                
            end
            Q=[BX,BXpb];
            p=(Q'*Q)\Q'*Y;
            a0=min(p(1:HT.nPrm(1),:),1-1e-6);
            m0=p(HT.nPrm(1)+1:end,:);

            %% variance | alp,mu,b
            R=(Y-Q*p)./(X.^b0);
            b0=ones(HT.nPrm(2),1).*b0;
            v0=ones(HT.nPrm(4),1).*var(R); %TODO make non stationary!!
           
            
            p0=[a0;b0;m0;v0];
                                          
            if HT.TwoParamFlag
                p0=p0(HT.PrmInd<=2);                                               
            end
        end
        
        function [Alp,Bet,M,V]=p2Prm(HT,p)
            %convert optimization p to separate parameters
            nAsc=HT.nDmn-1;
            if HT.TwoParamFlag                
                Alp=p(HT.PrmInd==1,:);
                Bet=p(HT.PrmInd==2,:);
                M=zeros(1,nAsc-1);
                V=ones(1,nAsc-1); %variance
            else                
                Alp=p(HT.PrmInd==1,:); %nBin x nD-1
                Bet=p(HT.PrmInd==2,:); %1 xnD-1
                M=p(HT.PrmInd==3,:);  %1 x nD -1
                V=p(HT.PrmInd==4,:); % 1 x nD-1
            end            
        end
        
        function PLOGL=likelihood(HT,X,Y,A,p,L)
            %function PLOGL=likelihood(HT,X,Y,A,p,L)
            %compute penalised Heffernan and Tawn likelihood for data on Gumbel scale.
            %INPUT
            %-X         n x 1 conditioned value
            %-Y         n x 1 x (nDmn-1) conditioning value
            %-A         n x 1 x (nDmn-1) bin allocation
            %-p         nBin+3 x (nDmn-1) parameter values
            %-L         smoothness parameter (alpha)
            %-MrgTp     Gumbel or Laplace
            %-TwoPrmFlag TwoParameters HT Model if true, else 4 parameter
            if nargin<=5
                L=0; %no penalty
            end
                                      
            if HT.checkParam(p)
                PLOGL=Inf;
                return
            end
             
            [Alp,Bet,M,V]=HT.p2Prm(p);
                            
            Alp_Dat = Alp(A(:,1),:);                          
            Bet_Dat = Bet(A(:,2),:);                      
            M_Dat = M(A(:,3),:);                                  
            V_Dat = V(A(:,4),:);  %B*Alp
                                    
            Xb=bsxfun(@power,X,Bet_Dat);  %X^b
            Std=bsxfun(@times,sqrt(V_Dat),Xb); %standard deviation
            
            Y = reshape(Y,[],HT.nDmn-1);           
            
            NLOGL=sum(sum(0.5*((Y-bsxfun(@times,Alp_Dat,X)-bsxfun(@times,M_Dat,Xb)).^2./(Std).^2)+log(Std)));
            
            if HT.nBin>1
                PLOGL=NLOGL+L.*sum(sum((Alp-mean(Alp,1)).^2+sum((Bet-mean(Bet,1)).^2)...
                    +sum((M-mean(M,1)).^2)+sum((V-mean(V,1)).^2),1),2);     %Penalised Likelhiood
            else
                PLOGL=NLOGL;
            end
        end %likelihood
        
        %
        function [PLOGL,G]=myobj(HT,X,Y,A,p,L)
            %my objective function for fmincon
            PLOGL=likelihood(HT,X,Y,A,p,L);
            G=HT.mygradfmincon(X,Y,A,p,L);
        end
               
        function Bad=checkParam(HT,varargin)
            %Check parameter values are in domain
            %TODO add Keef constraints??
            Bad=false;
            
            [Alp,Bet,~,V]=HT.p2Prm(varargin{:});
            
            %nD checks need unwrap Alp(:)
            if any(Alp(:)>1) || any(Alp(:)<-1) || any(V(:)<0) || any(Bet(:)>1)  %constraints
                Bad=true;
                return
            end
            
            if strcmp(HT.MarginType,'Gumbel')
                if any(Alp(:)<0)
                    Bad=true;
                    return
                end
            end
        end %checkParam
        
        function [G,H]=Gradient(HT,iP,X,Y,A,p,L)
            %[G,H]=Gradient(iP,X,Y,A,p,L,TwoPrmFlg)
            %compute gradients for the Hefferenan and Tawn model
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Title:    Gradient and Expected Fisher information. Vectorized and
            %           suitable for splines.
            % Author:   Thijs Willems (graduate intern TU Delft)
            % Date:     17-10-2016
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin<=5
                L=0; %no penalty
            end
            
            
            [Alp,Bet,M,V]=HT.p2Prm(p);
            
            if HT.nBin==1
                L=0; %switch off penalty
            end
            
            % Get paramters at data
            Alp_Dat = Alp(A(:,1),:);
            Bet_Dat = Bet(A(:,2),:);         
            if HT.TwoParamFlag
                M_Dat=0;
                V_Dat=1;
            else
                M_Dat = M(A(:,3),:);
                V_Dat = V(A(:,4),:);  %B*Alp
            end
            
            %% First Derivatives
            switch iP
                case 1 %alpha and mu
                    G = NaN(sum(HT.nPrm([1,3])),HT.nDmn-1);
                    %dL/dalp
                    t1=-1./V_Dat.*X.^(1-2.*Bet_Dat).*(Y - Alp_Dat.*X - X.^Bet_Dat.*M_Dat);
                    
                    RsdAlp=(Alp-mean(Alp));%residual alpha difference (penalty term)
                    for iDmn=1:HT.nDmn-1
                        G(1:HT.nPrm(1),iDmn) = accumarray(A(:,1),t1(:,iDmn),[HT.nPrm(1),1],@sum,0)+2.*L.*RsdAlp(:,iDmn);
                    end
                    
                    %dL/dmu
                    if HT.TwoParamFlag
                        G(HT.nPrm(1)+1,:)=[];
                    else
                        t1=-(1./V_Dat).*X.^(-Bet_Dat).*(Y - Alp_Dat.*X - X.^Bet_Dat.*M_Dat);
                        
                        RsdM=(M-mean(M)); %residual mu difference (penalty term)
                        for iDmn=1:HT.nDmn-1
                            G(HT.nPrm(1)+1:end,iDmn) = accumarray(A(:,3),t1(:,iDmn),[HT.nPrm(3),1],@sum,0)+2.*L.*RsdM(:,iDmn);
                        end
                        
                    end
                case 2 %beta and sigma
                    G = NaN(sum(HT.nPrm([2,4])),HT.nDmn-1);
                    %dL/dbeta
                    t1= -(1./V_Dat).*X.^(-2.*Bet_Dat).*log(X).*(Y - Alp_Dat.*X).*(Y - Alp_Dat.*X - X.^Bet_Dat.*M_Dat) + log(X);
                    RsdBet=(Bet-mean(Bet));%residual alpha difference (penalty term)
                    for iDmn=1:HT.nDmn-1
                        G(1:HT.nPrm(2),iDmn) = accumarray(A(:,2),t1(:,iDmn),[HT.nPrm(2),1],@sum,0)+2.*L.*RsdBet(:,iDmn);
                    end
                    
                    %dL/dsigma
                    if HT.TwoParamFlag
                        G(HT.nPrm(2)+1,:)=[];
                    else
                        % G(2,:) = sum(-1./S.*(1./(2.*S).*X.^(-2*Bet).*(Y - Al.*X - X.^Bet.*M).^2 - 1/2));
                        t1=(1./(2.*V_Dat)).*(1-(1./V_Dat).*((Y - Alp_Dat.*X - X.^Bet_Dat.*M_Dat)./(X.^Bet_Dat)).^2);
                        RsdV=(V-mean(V)); %residual V difference (penalty term)
                        for iDmn=1:HT.nDmn-1
                            G(HT.nPrm(2)+1:end,iDmn) = accumarray(A(:,4),t1(:,iDmn),[HT.nPrm(4),1],@sum,0)+2.*L.*RsdV(:,iDmn);
                        end
                    end
            end
            
            %% Expected Fisher information matrix
            if nargout==2
                H=NaN(size(G,1),size(G,1),size(G,2));
                switch iP
                    case 1 %alpha and mu
                        t1=1./V_Dat.*X.^(2-2*Bet_Dat);
                        for iDmn=1:HT.nDmn-1 %loop iDmn
                            E11 = diag(accumarray(A(:,1),t1(:,iDmn),[HT.nPrm(1),1],@sum,0))+2.*L.*(HT.nPrm(1)-1).*(HT.nPrm(1)-2)./(HT.nPrm(1).^2);
                            if  HT.TwoParamFlag
                                H(:,:,iDmn)=E11;
                            else
                                E33 =diag(accumarray(A(:,3),1./V_Dat(:,iDmn),[HT.nPrm(3),1],@sum,0))+2.*L.*(HT.nPrm(3)-1).*(HT.nPrm(3)-2)./(HT.nPrm(3).^2);
                                t1=(1./V_Dat(:,iDmn)).*X.^(1-Bet_Dat(:,iDmn));
                                if HT.NonStat(1) &&  HT.NonStat(3) %both non stationary
                                    E13 = diag(accumarray(A(:,1),t1,[HT.nPrm(1),1],@sum,0)); % nBin  x nBin
                                elseif HT.NonStat(1)
                                    E13 = accumarray(A(:,1),t1,[HT.nPrm(1),1],@sum,0); % nBin  x 1
                                elseif HT.NonStat(3)
                                    E13 = accumarray(A(:,3),t1,[HT.nPrm(3),1],@sum,0)'; % 1 x nBin
                                else %both stationary
                                    E13 = sum(t1); % 1 x 1
                                end
                                H(:,:,iDmn)=[E11,E13;E13',E33]; %[nBin +1 x nBin +1]
                            end
                        end %iDmn
                    case 2 %beta and sigma
                        %d2L/DLbeta2
                        t1=(log(X).^2).*(2 + M_Dat.^2./V_Dat);
                        for iDmn=1:HT.nDmn-1 %loop iDmn
                            E22 = diag(accumarray(A(:,2),t1(:,iDmn),[HT.nPrm(2),1],@sum,0))+2.*L.*(HT.nPrm(2)-1).*(HT.nPrm(2)-2)./(HT.nPrm(2).^2);
                            
                            if  HT.TwoParamFlag
                                H(:,:,iDmn)=E22;
                            else
                                
                                %d2L/DLsigma2
                                E44 =diag(accumarray(A(:,4),(1./(2.*V_Dat(:,iDmn).^2)),[HT.nPrm(4),1],@sum,0))+2.*L.*(HT.nPrm(4)-1).*(HT.nPrm(4)-2)./(HT.nPrm(4).^2);
                                
                                %d2L/DLsigma2
                                t1=(1./V_Dat).*log(X);
                                if HT.NonStat(2) &&  HT.NonStat(4) %both non stationary
                                    E24 = diag(accumarray(A(:,2),t1,[HT.nPrm(2),1],@sum,0)); % nBin  x nBin
                                elseif HT.NonStat(2)
                                    E24 = accumarray(A(:,2),t1,[HT.nPrm(2),1],@sum,0); % nBin  x 1
                                elseif HT.NonStat(4)
                                    E24 = accumarray(A(:,4),t1,[HT.nPrm(4),1],@sum,0)'; % 1 x nBin
                                else %both stationary
                                    E24 = sum(t1); % 1 x 1
                                end
                                
                                H(:,:,iDmn)=[E22,E24;E24',E44]; %[2 x 2]
                            end
                        end %iDmn
                end
            end
        end%Gradient
    end %Methods private
    methods(Static)

        
         function G=mygradfmincon(X,Y,A,p,L)
            %G=mygradfmincon(X,Y,A,p,L)
            %compute gradients for the Hefferenan and Tawn model
            %INPUT
            %-X         n x 1 conditioned value
            %-Y         n x 1 x (nDmn-1) conditioning value
            %-A         n x 1 x (nDmn-1) bin allocation
            %-p         nBin+3 x (nDmn-1) parameter values
            %-L         smoothness parameter (alpha)
            %
            %OUTPUT
            % G       p x (nDmn-1)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Title:    Gradient and Expected Fisher information. Vectorized and
            %           suitable for splines.
            % Author:   Thijs Willems (graduate intern TU Delft)
            % Date:     17-10-2016
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin<=4
                L=0; %no penalty
            end
            n = size(p,1);
            nDmn = size(p,2)+1;
            nBin=n-3;
            Alp=p(1:nBin,:); %nBin x nD-1
            Bet=p(nBin+1,:); %1 xnD-1
            M=p(nBin+2,:);  %1 x nD -1
            V=p(nBin+3,:); % 1 x nD-1  %this the variance!! in the gradient
            
            if nBin==1
                L=0; %switch off penalty
            end
            
            % Initliatlize
            
            Al = Alp(A,:); %expand alpha for all bins
            
            AlpRsd=(Alp-mean(Alp,1)); %residual alpha difference (penalty term)
            
            %% First Derivatives
            G = NaN(nBin+3,nDmn-1);
            %dL/dalp
            t1=(-1./V).*X.^(1-2.*Bet).*(Y - Al.*X - X.^Bet.*M);
            for iDmn=1:nDmn-1
                G(1:nBin,iDmn) = accumarray(A,t1(:,iDmn),[nBin,1],@sum,0)+2.*L.*AlpRsd(:,iDmn);
            end            
            %dL/dbeta
            G(nBin+1,:) = sum((-1./V).*X.^(-2.*Bet).*log(X).*(Y - Al.*X).*(Y - Al.*X - X.^Bet.*M) + log(X));
            %dL/dmu
            G(nBin+2,:) = sum((-1./V).*X.^(-Bet).*(Y - Al.*X - X.^Bet.*M));            
            %dL/dsigma            
            %G(nBin+3,:) = sum(-1./S.*(1./(2.*S).*X.^(-2*Bet).*(Y - Al.*X - X.^Bet.*M).^2 - 1/2));
            G(nBin+3,:) = (1./(2.*V)).*sum(1-(1./V).*((Y - Al.*X - X.^Bet.*M)./(X.^Bet)).^2);
                                                            
        end%mygradfmincon
        
         function [c,ceq]=nonlincon(p)
            n = size(p,1);
            nBin=n-3;
            Alp=p(1:nBin,:); %nBin x nD-1
            Bet=p(nBin+1,:); %1 xnD-1            
            S=p(nBin+3,:); % 1 x nD-1
            
            c=[Alp-1;
                -1-Alp;
                Bet-1;
                -S];
            ceq=0;
        end
    end %methods static
end %class