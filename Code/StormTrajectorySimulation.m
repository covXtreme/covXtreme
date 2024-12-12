% Copyright © [2023] Shell Global Solutions International B.V. All Rights Reserved.
% SPDX-License-Identifier: Apache-2.0

classdef StormTrajectorySimulation
    properties
        RA % nExc x nDmn (cell array) associated variables storm trajectory
        Cvr % nExc x nDmn (cell array) covariate storm trajectory
        Sml % nSml x 1 (cell array) set of simulations of the event set
        A % nExc x nAsc (cell array) associated variables storm trajectory
        
        nCvr      % 1 x 1, number of covariates
        nDat      % 1 x 1, number of data obs
        nAsc      % 1 x 1, number of associated variables
        
        nSml=10; % 1 x 1, number of simulations
        nNgh=10; %Number of near neighbour historical storm trajectories to consider
        RtrPrd=100; % nRtr x 1 Return Period  in years
        DistanceMetric='Square'; % distance metric either Square or Absolute 
    end %properties
    
    methods
        function obj=StormTrajectorySimulation(Rsp,Asc,Cvr,Dat,Bn,Mrg,Opts)
            % wrapper function for StormTrajectorySimulation class
            %INPUT
            % Rsp     [n x 1] vector of main response 
            % Cvr     [n x nCvr] matrix of covariates (e.g. direction, season)
            % Asc     [n x nAsc] matrix of other associated variables 
            % Dat structure  (from stage 1)
            %     - Dat.X     nDat x nCvr  covariate values
            %     - Dat.Y     nDat x nDmn  response data (Y)
            %     - Dat.IsPrd   nCvr x 1 vector if covairte is periodic or not
            %     - Dat.CvrLbl    char string for response label
            %     - Dat.RspLbl    char string for response label
            % Bn bin allocation structure from stage 2
            % Mrg 2 x 1, marginal model structure *output from stage 3
            % Opts Options for the storm trajectory class object 
            if nargin == 0
                return
            end
            if nargin>=7  && ~isempty(Opts)
                if isfield(Opts,'nSml')
                    validateattributes(Opts.nSml, {'numeric'},{'scalar','positive'},'StormTrajectorySimulation','Opts.nSml',7);
                    obj.nSml=Opts.nSml;
                end %nSml
                if isfield(Opts,'nNgh')
                    validateattributes(Opts.nNgh, {'numeric'},{'scalar','positive'},'StormTrajectorySimulation','Opts.nNgh',7);
                    obj.nNgh=Opts.nNgh;
                end %nNgh
                if isfield(Opts,'RtrPrd')
                    validateattributes(Opts.RtrPrd, {'numeric'},{'scalar','positive'},'StormTrajectorySimulation','Opts.RtrPrd',7);
                    obj.RtrPrd=Opts.RtrPrd;
                end %nNgh
                if isfield(Opts,'DistanceMetric')
                    validateattributes(DistanceMetric,{'string','char'},{},'StormTrajectorySimulation','Opts.DistanceMetric',7)
                    obj.DistanceMetric = validatestring(Opts.DistanceMetric,{'Square','Absolute'});
                end
            end
            [obj.nDat,obj.nCvr]=size(Dat.X);
            if ~isa(Mrg,'MarginalModel')
                error('Mrg should be a nDmn x 1 Marginal Model')
            end
            if numel(Mrg)<=1
                error('Mrg should be a nDmn x 1 Marginal Model')
            end
            
            obj=IdentifyStormTrajectories(obj,Dat,Rsp,Asc,Cvr);
                       
            obj=CreateStormTrajectoryBinAllocation(obj,Dat,Bn);
                      
            obj=SimulateEventSet(obj,Mrg);
            
            obj=StormMatching(obj,Dat,Bn);
        end %StormTrajectorySimulation
        
        function obj=IdentifyStormTrajectories(obj,Dat,Rsp,Asc,Cvr)
            % obj=IdentifyStormTrajectories(obj,Dat)
            % Function to separate out each storm trajectory
            % INPUT
            % Dat structure from stage 1
            % Rsp     [n x 1] vector of main response 
            % Cvr     [n x nCvr] matrix of covariates (e.g. direction, season)
            % Asc     [n x nAsc] matrix of other associated variables 
            obj.nAsc=size(Asc,2);
            obj.RA=cell(obj.nDat,1+obj.nAsc); %Initialise empty cell array for storm trajectories
            obj.Cvr=cell(obj.nDat, obj.nCvr); %Initialise empty cell array for covariate trajectories.
            
            %% For each storm peak, get its trajectory.
            for iDat = 1:obj.nDat
                % Find the start and end of each storm
                startIdx = Dat.Prd(iDat, 1); % start of the storm
                endIdx = Dat.Prd(iDat, 2);   % end of the storm
                
                % Extract the storm trajectory within the storm period
                stormRsp = Rsp(startIdx:endIdx);  % Response trajectory
                stormAsc = Asc(startIdx:endIdx, :);  % Associated variables trajectory
                stormCvr = Cvr(startIdx:endIdx, :); % Covariates
                
                % Store the storm trajectory
                obj.RA{iDat, 1} = stormRsp; % Store response trajectory
                for iA = 1:obj.nAsc
                    obj.RA{iDat, iA+1} = stormAsc(:, iA);  % Store associated variables trajectory
                    obj.Cvr{iDat, iA} = stormCvr(:, iA);
                end %iAsc
            end %iDat
            
            TrajectoryPlot(obj, Dat, [], 'Stg6_StormTrajectory');
                       
        end %IdentifyStormTrajectories
        
        function obj=CreateStormTrajectoryBinAllocation(obj,Dat,Bn)
            % CreateStormTrajectoryBinAllocation(obj)
            % from the bin allocation determine the bins of the historical
            % trajectories 
            % INPUTS
            % Dat structure from stage 1
            % Bn bin allocation structure from stage 2
            fprintf(1,'Calculating bin allocations for historical trajectories.\n');
            obj.A=cell(Bn.n,obj.nAsc);
            for iBn=1:Bn.n
                tCvr=zeros(Bn.n,1);
                tStrTrjCvr=obj.Cvr{iBn,:};
                tCvr(1:size(tStrTrjCvr,1))=tStrTrjCvr;
                tBn=BinAllocation(Bn,tCvr, 0);
                obj.A{iBn,:}=tBn.A(1:size(tStrTrjCvr,1));
            end %iS
            
            TrajectoryPlot(obj, Dat, Bn, 'Stg6_StormTrajectory');
        end %CreateStormTrajectoryBinAllocation
        
        function obj=SimulateEventSet(obj, Mrg)
            %obj=SimulateEventSet(obj, Mrg)
            % simulate a number of trajectories
            % INPUT
            % Mrg 2 x 1, marginal model structure *output from stage 3
            %poissrnd(Mrg(1).nDat*obj.RtrPrd/Mrg(1).Yrs,10,1)
            % TODO move from cell structure to tabular 
            nStrVec=poissrnd(Mrg(1).nDat*obj.RtrPrd/Mrg(1).Yrs,10,1);
            nStr=sum(nStrVec);
            obj.Sml=struct;
            obj.Sml.SmlIndex=repelem(1:10,nStrVec)';
            [obj.Sml.A,obj.Sml.Org]=deal(NaN(nStr,1));
            for iS=1:obj.nSml
                tSml=sample_MC(Mrg(1),nStrVec(iS));
                obj.Sml.A(SmlIndex==iS,:)=tSml.A;
                obj.Sml.Org(SmlIndex==iS,:)=tSml.Org;
            end %iS
            for iS=1:obj.nSml
                nStr = poissrnd(Mrg(1).nDat*obj.RtrPrd/Mrg(1).Yrs);
                obj.Sml{iS}=sample_MC(Mrg(1),nStrV); %This is input to storm matching
            end %iS
        end %SimulateEventSet
        
        function obj=StormMatching(obj,Dat,Bn)
            %% Storm matching to allocate historical trajectories to simulated storms
            % Match using the dominant variate only, and in the covariate bin
            % corresponding to the dominant variate only.
            % INPUT
            % Dat peak picked data structure from stage 1
            % Bn bin allocation structure from stage 2
            fprintf(1,'Storm matching to allocate historical trajectories to simulated storms.\n');
            iE=0;
            for iS=1:obj.nSml
                for iR=1:obj.Sml{iS}.nRls
                    iE=iE+1;
                    tBn=obj.Sml{iS}.A(iR); %bin for realisation iR from simulation iS
                    tMB=find(Bn.A==tBn);%all historical matched with correct bin allocation
                    difference_matched=Dat.Y(tMB,1)-obj.Sml{iS}.Org(iR);
                    switch obj.DistanceMetric
                        case 'Square'
                            [~,tD2]=sort((difference_matched).^2);
                        case 'Absolute'
                            [~,tD2]=sort(abs(difference_matched));
                        otherwise
                            error('Distance metric not recognised')
                    end
                    tMV=tMB(tD2(1:obj.nNgh));
                    tMtc=tMV(randi(obj.nNgh,1));
                    obj.Sml{iS}.StrTrj.RA{iR,:}=obj.RA{tMtc,:};
                    obj.Sml{iS}.StrTrj.Cvr{iR,:}=obj.Cvr{tMtc,:};
                    obj.Sml{iS}.StrTrj.A{iR,:}=obj.A{tMtc,:};
                end %iR
            end %iS
            
           %TrajectoryPlot(obj,Dat, [], 'Stg6_Data_StormTrajectory_Simulated');
            
        end %StormMatching
        
        function obj=TrajectoryPlot(obj, Dat, Bn, FilNam)
            %% Storm trajectory plot
            % Match using the dominant variate only, and in the covariate bin
            % corresponding to the dominant variate only.
            % INPUT
            % Dat structure from stage 1
            % Bn structure from stage 2
            % FlNm character string for the filename
            fprintf(1,'Creating trajectory plot.\n');
            if isempty(Bn)
                nAscp1 = size(obj.RA, 2);               
                for iA = 1:nAscp1
                    % Plot without covariate
                    tFilNam = sprintf('%s_%s_Time', FilNam, Dat.RspLbl{iA});
                    obj.plotStormTrajectories(Dat.Y, obj.RA(:, iA), [], tFilNam, Dat.RspLbl{iA}, '');
                    
                    for iC = 1:obj.nCvr
                        tFilNam = sprintf('%s_%s_%s', FilNam, Dat.RspLbl{iA}, Dat.CvrLbl{iC});
                        obj.plotStormTrajectories(Dat.Y, obj.RA(:, iA), obj.Cvr(:, iC), tFilNam, Dat.RspLbl{iA}, Dat.CvrLbl{iC});
                    end %iC
                end %iA
            else
                for iB=1:Bn.nBin
                    nAscp1 = size(obj.RA, 2);                   
                    for iA = 1:nAscp1
                        % Plot without covariate
                        tFilNam = sprintf('%s_%s_Bin%g_Time', FilNam, Dat.RspLbl{iA}, iB);
                        tRsp = sprintf('%s:%s', Dat.RspLbl{iA}, Bn.BinLbl{iB});
                        obj.plotStormTrajectories(Dat.Y(Bn.A==iB,1), obj.RA(Bn.A==iB, iA), [], tFilNam, tRsp, '');
                        
                        for iC = 1:obj.nCvr
                            tFilNam = sprintf('%s_%s_Bin%g_%s', FilNam, Dat.RspLbl{iA}, iB, Dat.CvrLbl{iC});
                            tRsp = sprintf('%s:%s', Dat.RspLbl{iA}, Bn.BinLbl{iB});
                            obj.plotStormTrajectories(Dat.Y(Bn.A==iB,1), obj.RA(Bn.A==iB, iA), obj.Cvr(Bn.A==iB, iC), tFilNam, tRsp, Dat.CvrLbl{iC});
                        end%iC
                    end %iA
                end %iB
            end
        end
        function obj=plotStormTrajectories(obj,Y, TrjRA, TrjCvr, tFilNam, LblRA, LblCvr)
            % Wrapper function to plot the storm trajectories
            [peakOrder, colorMap] = obj.getPeakOrderAndColorMap(Y(:, 1));
            clf;
            nPk=size(Y,1);
            for isNormalised = [true, false]
                subplot(2, 1, 1 + ~isNormalised)
                hold on;
                
                for iPk = 1:nPk
                    peakIndex = peakOrder(iPk);
                    stormRsp = TrjRA{peakIndex, 1};
                    if any(isnan(stormRsp))
                        continue
                    else
                        peakValue = max(stormRsp);
                        normalisedRsp = stormRsp / peakValue;
                        
                        if ~isempty(TrjCvr)
                            stormCvr = TrjCvr{peakIndex, 1};
                            peakCvr = stormCvr(find(stormRsp == peakValue, 1));
                            CntStrCvr = obj.adjustCovariate(stormCvr - peakCvr);
                            xData = CntStrCvr;
                        else
                            lenBefore = find(stormRsp == peakValue, 1) - 1;
                            lenAfter = length(stormRsp) - lenBefore - 1;
                            xData = (-lenBefore:lenAfter)';
                        end
                        
                        yData = normalisedRsp;
                        if ~isNormalised, yData = stormRsp; end
                        plot(xData, yData, 'Color', colorMap(iPk, :), 'LineWidth', 0.5, 'LineStyle', ':', 'Marker', '.', 'MarkerSize', 10);
                    end
                end%iPk
                
                % Axis labels and title
                [titleStr, xlabelStr, ylabelStr, xlimits] = obj.getPlotProperties(isNormalised, TrjCvr, LblRA, LblCvr);
                xlabel(xlabelStr);
                ylabel(ylabelStr);
                title(titleStr);
                grid on;
                box on;
                axis tight;
                if ~isempty(xlimits)
                    xlim(xlimits);
                end
                hold off;
            end %isNormalised
            savePics(fullfile('Figures', tFilNam));
        end %plotStormTrajectories
        
    end %methods
    
    methods (Static)
        function [peakOrder, colorMap] = getPeakOrderAndColorMap(peakValues)
            % Function to select colour ranges depending on the number of
            % peaks
            % INPUTS
            % peakValues  nDat x 1 array of storm peak values 
            [~, peakOrder] = sort(peakValues, 'descend');
            colorMap = jet(size(peakValues,1));
            colorMap = colorMap(end:-1:1, :);
        end %getPeakOrderAndColorMap
        
        function CntStrCvr = adjustCovariate(CntStrCvr)
            % Function to center covariate values for a storm peak
            % INPUTS
            % Covariate centered around the storm peak covariate
            CntStrCvr(CntStrCvr < -180) = CntStrCvr(CntStrCvr < -180) + 360;
            CntStrCvr(CntStrCvr > 180) = CntStrCvr(CntStrCvr > 180) - 360;
        end %adjustCovariate
        
        function [titleStr, xlabelStr, ylabelStr, xlimits] = getPlotProperties(isNormalised, TrjCvr, LblRA, LblCvr)
            % Generate properties for storm trajectory plotting 
            if isNormalised
                titleStr = sprintf('%s: Normalised Storm Trajectories', LblRA);
                ylabelStr = 'Normalised Response (Max=1)';
            else
                titleStr =sprintf('%s: Unnormalised Storm Trajectories', LblRA);
                ylabelStr = 'Unnormalised Response';
            end
            if isempty(TrjCvr)
                xlabelStr = 'Time relative to storm peak (t=0)';
                xlimits = [-48, 48];
            else
                xlabelStr = sprintf('Angle relative to storm peak (at 0) for covariate %s', LblCvr);
                xlimits = [];
            end
        end % getPlotProperties
        
    end %methods(Static)
    
end %classdef