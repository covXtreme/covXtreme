% Copyright 2023 covXtreme
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%      http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

classdef CovariateBinning
    %Allocate data to covariate bins
    %% INPUTS
    % Dat
    %    X [n x nCvr] covariate values (periodic covairates must be on on 0,360)
    %    Edg  cell array [nCvr x 1]  Covaraite Bin edges [0,360]; periodic bins will wrap around 0 degrees
    %    IsPrd [nCvr x 1] vector if covairte is periodic or not
    %% OPTIONAL INPUTS USED cFOR PLOTTING
    %    Y [n x nDmn] responses
    %    RspLbl  cell array of reponse variable labels
    %    CvrLbl  cell array of covariate variable labels
    %% OUTPUTS
    %    A  n x 1 bin allocation
    %    Plot of bin edges on data
    
    properties
        nCvr      %[1 x 1]      number of covariates
        n         %[1 x 1]      number of observations
        nBin      %[1 x 1]      total nunber of covariate bins
        nBinCvr   %[nCvr x 1]  number of covariate bins in each covariate
        Edg       %[nCvr x 1]   cell array Covaraite Bin edges [0,360]; periodic bins will wrap around 0 degrees        
        A         %[n x 1]      bin allocation of each observation to a bin
        Cnt       %[nBin x 1]   number of observations per bin
        BinLbl    %[nBin x 1]   Bin label
    end %properties
    
    properties (Hidden=true)
        LeftEdge   %[nBin x nCvr] Center of each bin
        Length     %[nBin x nCvr] Length of each bin
        EdgExd     %Extended padded with [0,360]
        IsPrd      %[nCvr x 1] logical of whether the covariate periodic        
        APrj       %[nBin x nCrv] Index matrix to go from full nBin x nDrc matrix over all covariates to indiviual covariate
        RspLbl     %[nRsp x 1] response label (used in plotting)
        CvrLbl     %[nCvr x 1] covariate label (used in plotting)
    end %properties(hidden)
        
    methods %methods
        
        function Bn=CovariateBinning(X,Edg,IsPrd,Y,RspLbl,CvrLbl)
            %Allocate data to covariate bins
            %% INPUTS
            % Dat
            %    X [n x nCvr] covariate values (periodic covairates must be on on 0,360)
            %    Edg  cell array [nCvr x 1]  Covariate Bin edges [0,360]; periodic bins will wrap around 0 degrees
            %    IsPrd [nCvr x 1] vector if covairte is periodic or not
            %Optional used for plotting
            %    Y [n x nDmn] responses
            %    RspLbl  cell array of reponse variable labels
            %    CvrLbl  cell array of covariate variable labels
            %% OUTPUTS
            %    Bn [n x 1] bin allocation
            %    Plot of bin edges on data
            
            if nargin==0
               return 
            end
            %% Input Check
            [Bn.n,Bn.nCvr]=size(X); %number of covariates
            % validation of attributes 
            validateattributes(Edg,{'cell'},{'numel',Bn.nCvr},'BinAllocation','Edg',2);
            Bn.Edg=Edg;
            validateattributes(IsPrd,{'numeric','logical'},{'numel',Bn.nCvr,'integer'},'BinAllocation','IsPrd',3);
            Bn.IsPrd=IsPrd;                        
            nDmn=size(Y,2); %number of reponses
            validateattributes(Y,{'numeric'},{'nrows',Bn.n},'BinAllocation','Y',4);
            validateattributes(RspLbl,{'cell'},{'numel',nDmn},'BinAllocation','RspLbl',6);
            Bn.RspLbl=RspLbl;
            validateattributes(CvrLbl,{'cell'},{'numel',Bn.nCvr},'BinAllocation','CvrLbl',7);
            Bn.CvrLbl=CvrLbl;
                                    
            %% Get bin allocation
            Bn=BinAllocation(Bn,X);
            
            %% Plot bins
            PlotBins(Bn,X,Y)
        end %CovariateBinning
        
        function Bn=BinAllocation(Bn,X)
            %% Loop over covariates and bin in each dimension
            %% INPUTS
            %    Bn [n x 1] bin allocation
            %% OUTPUTS
            %    Bn populated Bn object
            Bn.nBinCvr=NaN(Bn.nCvr,1); %number of bins
            ADmn=NaN(Bn.n, Bn.nCvr);  %bin allocation for each dimension
            for iC=1:Bn.nCvr %loop over covariates
                Bn.Edg{iC}=sort(Bn.Edg{iC}); %check edges are sorted
                
                if Bn.IsPrd(iC) %periodic
                    if any(X(:,iC)>360) || any(X(:,iC)<0)
                        error('Periodic X(:,%g) expected on range [0,360]',iC);
                    end
                    if any(Bn.Edg{iC}>360) || any(Bn.Edg{iC}<0)
                        error('Periodic Edg{%g} expected on range [0,360]',iC);
                    end
                    if numel(unique(mod(Bn.Edg{iC},360))) < numel(Bn.Edg{iC})
                        error('You have supplied 2 bin edges which correspond to same point on a circle. For a single bin, please provide 1D DrcEdg (any value will do). For 2 bins, please provide 2D DrcEdg containing values which differ on 360deg circle.');
                    end
                    Bn.nBinCvr(iC)=numel(Bn.Edg{iC}); %number of bins
                    if verLessThan('Matlab','8.5')
                        Bn.EdgExd{iC}=unique([0-1;Bn.Edg{iC};360+1]); %pad with 0, 360 if needed
                    else
                        Bn.EdgExd{iC}=unique([0;Bn.Edg{iC};360]); %pad with 0, 360 if needed
                    end
                else %non periodic
                    Bn.EdgExd{iC}=Bn.Edg{iC};
                    Bn.nBinCvr(iC)=numel(Bn.Edg{iC})-1; %number of bins
                    if any(X(:,iC)<Bn.Edg{iC}(1) | X(:,iC)>Bn.Edg{iC}(end))
                        error('Non Periodic X(:,%g) outside of range defined by bins',iC)
                    end
                end
                
                %find indices of directions X which fall in each directional bin
                ADmn(:,iC)=discretize(X(:,iC),Bn.EdgExd{iC});
                if Bn.IsPrd(iC)
                    ADmn(ADmn(:,iC)==Bn.nBinCvr(iC)+1,iC)=1; %wrap periodic bins around 0
                end
            end %iC
            
            %% combine bins over dimensions
            k=cumprod(Bn.nBinCvr);
            Bn.A=ADmn(:,1);
            for iC=2:Bn.nCvr
                Bn.A = Bn.A + (ADmn(:,iC)-1)*k(iC-1);
            end %iC
            
            Bn.nBin=prod(Bn.nBinCvr);
            Bn.APrj=cell(1,Bn.nCvr);
            [Bn.APrj{:}]=ind2sub(Bn.nBinCvr,(1:Bn.nBin)');
            Bn.APrj=cell2mat(Bn.APrj);
            
            Bn.BinLbl=cell(Bn.nCvr,1);
            Bn.LeftEdge= NaN(Bn.nBin,Bn.nCvr); %LeftEdge of the bins
            Bn.Length= NaN(Bn.nBin,Bn.nCvr); %length of bins
            for iC=1:Bn.nCvr                
                if Bn.IsPrd(iC) %period need to wrap bins
                    BinSt=Bn.Edg{iC};
                    BinEnd=circshift(Bn.Edg{iC},-1);
                    
                    if Bn.Edg{iC}(1) > 0 %if there is a bin which straddles 0, make it bin 1
                        BinSt=circshift(BinSt,1);
                        BinEnd=circshift(BinEnd,1);
                    end
                                        
                else %non-periodic
                    BinSt=Bn.Edg{iC}(1:end-1);
                    BinEnd=Bn.Edg{iC}(2:end);
                end                                
                
                for iB=1:Bn.nBin  %loop over bins to make label
                    
                    I=Bn.APrj(iB,iC);
                    
                    if BinEnd(I)<BinSt(I) && Bn.IsPrd(iC)  %periodic case
                        Bn.Length(iB,iC)=BinEnd(I)+360-BinSt(I);
                    else
                        Bn.Length(iB,iC)=BinEnd(I)-BinSt(I);     %nBin x nCvr Length of each bin
                    end
                    Bn.LeftEdge(iB,iC)=BinSt(I);    %nBin x nCvr LeftEdge of each bin
                    
                                        
                    if iC==1
                        Bn.BinLbl{iB}=sprintf('%s[%g, %g)',Bn.CvrLbl{iC}(1),BinSt(I),BinEnd(I));
                    else
                        Bn.BinLbl{iB}=[Bn.BinLbl{iB},' x ',sprintf('%s[%g, %g)',Bn.CvrLbl{iC}(1),BinSt(I),BinEnd(I))];
                    end
                    %% Store bin start and end points for plot labels
                end  %iB
            end  %iC          
            %% Ensure have no empty bins (or bins with <30 obs)          
            Bn.Cnt=accumarray(Bn.A,Bn.A,[Bn.nBin,1],@numel); %
            if any(Bn.Cnt < 30)
                warning('Too few exceedances in one or more bins. Consider reducing number of bins.');
                Bn.Cnt
            end
        end %BinAllocation
        
        function PlotBins(Bn,X,Y)
            %% Plotting function      
            
            %% Marginal plot
            figure(1)
            clf;
            nDmn=size(Y,2);
            c=0;
            for i=1:nDmn
                for iC=1:Bn.nCvr
                    c=c+1;
                    subplot(nDmn,Bn.nCvr,c)
                    plot(X(:,iC),Y(:,i),'k.')
                    axis tight
                    hold on
                    grid on
                    plot(repmat(Bn.Edg{iC},1,2)',ylim,'r--','linewidth',2)
                    xlabel(Bn.CvrLbl{iC})
                    if Bn.IsPrd(iC)
                        set(gca,'xtick',0:45:360,'xlim',[0,360])
                    end
                    ylabel(Bn.RspLbl{i})
                end  %iC
            end %i
            savePics('Figures/Stg2_Data_BinEdges')        
                        
            %% Joint plot
                  
            for iD=2:nDmn                                                    
                figure(iD+1)    
                clf;
                for iC=1:Bn.nBin 
                    if Bn.nBin > 1  %if more than one bin, use subplots
                        subplot(2,ceil(Bn.nBin/2),iC)
                    end
                    plot(Y(Bn.A==iC,1),Y(Bn.A==iC,iD),'k.')
                    axis tight
                    hold on
                    grid on
                    xlabel(Bn.RspLbl{1})
                    ylabel(Bn.RspLbl{iD})
                    if Bn.nBin > 1
                        title(Bn.BinLbl{iC})
                    end
                end  %iC
                savePics(sprintf('Figures/Stg2_Data_BinScatterPlot_Y%g_Y1',iD))
            end  %iD       
            
        end %PlotBins
        
        function PlotBinEdge(Bn,iC)
            %add vertical lines to plot with bin boundarys for covariate iC
            ty=ylim;
            hold on
            for i=1:Bn.nBinCvr(iC)  %loop over bins
                plot(Bn.Edg{iC}(i)*[1,1],ylim,'r--','linewidth',1)
            end %i     
            ylim(ty);
        end %PlotBinEdge
                     
        function X=SampleCovariateFromBin(Bn,A)
            %given bin number sample covariate uniformly from within bin  
            X=rand(numel(A),Bn.nCvr).*Bn.Length(A,:)+ Bn.LeftEdge(A,:);  %move covariate value within bin
            for iC=1:Bn.nCvr  %loop over covariates
                if Bn.IsPrd(iC)
                   X=mod(X,360);
                end
            end %iC
        end %SampleCovariateFromBin
        
        function PlotParameter(Bn, P,iC,varargin)
            %Project and plot parameter for ith Covariates
            %% INPUTS
            %P parameters
            %iC covariate dimension
            % optional plotting arguements                        
            
            %% Take P and project onto iC dimension    
            if Bn.nBin==1
                PPrj=mean(P);
            else
                PPrj=grpstats(P,Bn.APrj(:,iC),{'mean'});
            end
            
            X=linspace(min(Bn.EdgExd{iC}),max(Bn.EdgExd{iC}),360);
            
            if Bn.IsPrd(iC)
               tX=mod(X,360); 
            else
                tX=X;
            end
                     
            %find indices of directions X which fall in each directional bin
            tA=discretize(tX,Bn.EdgExd{iC});
            if Bn.IsPrd(iC)
                tA(tA==Bn.nBinCvr(iC)+1)=1; %wrap periodic bins around 0
            end
                         
            if size(P,2)>1
                qP=quantile(PPrj,[0.025,0.5,0.975],2);
                plot(X,qP(tA,2),'-',varargin{:})
                hold on
                plot(X,qP(tA,[1,3]),'--',varargin{:})
            else
                plot(X,PPrj(tA,1),'-',varargin{:})
            end
                
            if Bn.IsPrd(iC)
               xlim([0,360]);
               set(gca,'xtick',0:90:360);
            end
                               
        end %PlotParameter
   
    end %methods
    
end %class
