% Copyright © [2023] Shell Global Solutions International B.V. All Rights Reserved.
% SPDX-License-Identifier: Apache-2.0

classdef StormTrajectorySimulation
    properties
        RA % nExc x nDmn (cell array) associated variables storm trajectory
        Cvr % nExc x nDmn (cell array) covariate storm trajectory
        RtrPrd=100; % nRtr x 1 Return Period  in years
        Sml % nSml x 1 (cell array) set of simulations of the event set
        nSml=10; % 1 x 1, number of simulations
        nNgh=10; %Number of near neighbour historical storm trajectories to consider
        A % nExc x nAsc (cell array) associated variables storm trajectory
        
        nCvr      % 1 x 1, number of covariates
        nDat      % 1 x 1, number of data obs
        nAsc      % 1 x 1, number of associated variables
    end %properties
    
    methods
        function obj=StormTrajectorySimulation(Rsp,Asc,Cvr,Dat,Bn,Mrg)
            % wrapper function for StormTrajectorySimulation class
            %INPUT
            % Rsp
            % Asc
            % Cvr
            % Dat structure  (from stage 1)
            %     - Dat.X     nDat x nCvr  covariate values
            %     - Dat.Y     nDat x nDmn  response data (Y)
            %     - Dat.IsPrd   nCvr x 1 vector if covairte is periodic or not
            %     - Dat.CvrLbl    char string for response label
            %     - Dat.RspLbl    char string for response label
            % Mrg 2 x 1, marginal model structure *output from stage 3
            if nargin == 0
                return
            end
            [obj.nDat,obj.nCvr]=size(Dat.X);
            if ~isa(Mrg,'MarginalModel')
                error('Mrg should be a nDmn x 1 Marginal Model')
            end
            if numel(Mrg)<=1
                error('Mrg should be a nDmn x 1 Marginal Model')
            end
            
            obj=IdentifyStormTrajectories(obj,Dat,Rsp,Asc,Cvr);
            
            obj=CreateStormTrajectoryBinAllocation(obj,Bn);
            
            obj=SimulateEventSet(obj,Mrg);
            
            obj=StormMatching(obj,Dat,Bn);
        end %StormTrajectorySimulation
        
        function obj=IdentifyStormTrajectories(obj,Dat,Rsp,Asc,Cvr)
            % obj=IdentifyStormTrajectories(obj,Dat)
            % Function to separate out each storm trajectory
            % INPUT
            % Dat structure from stage 1
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
            
            TrajectoryPlot(obj,Dat.Y, obj, Dat.RspLbl, Dat.CvrLbl, [], 'Stg6_Data_StormTrajectory')
            
        end %IdentifyStormTrajectories
        
        function obj=CreateStormTrajectoryBinAllocation(obj,Bn)
            % CreateStormTrajectoryBinAllocation(obj)
            fprintf(1,'Calculating bin allocations for historical trajectories\n');
            obj.A=cell(Bn.n,obj.nAsc);
            for iBn=1:Bn.n
                tCvr=zeros(Bn.n,1);
                tStrTrjCvr=obj.Cvr{iBn,:};
                tCvr(1:size(tStrTrjCvr,1))=tStrTrjCvr;
                tBn=BinAllocation(Bn,tCvr, 0);
                obj.A{iBn,:}=tBn.A(1:size(tStrTrjCvr,1));
            end %iS
        end %CreateStormTrajectoryBinAllocation
        
        function obj=SimulateEventSet(obj, Mrg)
            %obj=SimulateEventSet(obj, Mrg)
            % INPUT
            % Mrg 2 x 1, marginal model structure *output from stage 3
            obj.Sml=cell(obj.nSml,1);
            for iS=1:obj.nSml
                nStr = poissrnd(Mrg(1).nDat*obj.RtrPrd/Mrg(1).Yrs);
                obj.Sml{iS}=sample_MC(Mrg(1),nStr); %This is input to storm matching
            end %iSml
        end %SimulateEventSet
        
        function obj=StormMatching(obj,Dat,Bn)
            %% Storm matching to allocate historical trajectories to simulated storms
            % Match using the dominant variate only, and in the covariate bin
            % corresponding to the dominant variate only.
            % INPUT
            % Dat structure from stage 1
            % Bn structure from stage 2
            fprintf(1,'Storm matching to allocate historical trajectories to simulated storms:\n');
            iE=0;
            for iS=1:obj.nSml
                for iR=1:obj.Sml{iS}.nRls
                    iE=iE+1;
                    tBn=obj.Sml{iS}.A(iR); %bin for realisation iR from simulation iS
                    tMB=find(Bn.A==tBn); %all historical matched with correct bin allocation
                    [~,tD2]=sort((Dat.Y(tMB,1)-obj.Sml{iS}.Org(iR)).^2);
                    tMV=tMB(tD2(1:obj.nNgh));
                    tMtc=tMV(randi(obj.nNgh,1));
                    obj.Sml{iS}.StrTrj.RA{iR,:}=obj.RA{tMtc,:};
                    obj.Sml{iS}.StrTrj.Cvr{iR,:}=obj.Cvr{tMtc,:};
                    obj.Sml{iS}.StrTrj.A{iR,:}=obj.A{tMtc,:};
                    %Dgn(iE,:)=[max(StrTrj.RA{tMtc,1}) obj.Sml{iS}.Org(iR) obj.Sml{iS}.A(iR)];
                end %iR
            end %iS
            
        end %StormMatching
        
        function obj=TrajectoryPlot(obj, Y, Trj, RspLbl, CvrLbl, Bn, FilNam)
            
            if isempty(Bn)
                nAscp1 = size(Trj.RA, 2);
                nCvr = size(Trj.Cvr, 2);
                
                for iA = 1:nAscp1
                    % Plot without covariate
                    tFilNam = sprintf('%s_%s_Time', FilNam, RspLbl{iA});
                    obj.plotStormTrajectories(Y, Trj.RA(:, iA), [], tFilNam, RspLbl{iA}, '');
                    
                    for iC = 1:nCvr
                        tFilNam = sprintf('%s_%s_%s', FilNam, RspLbl{iA}, CvrLbl{iC});
                        obj.plotStormTrajectories(Y, Trj.RA(:, iA), Trj.Cvr(:, iC), tFilNam, RspLbl{iA}, CvrLbl{iC});
                    end
                end
            else
                
                for iB=1:Bn.nBin
                    nAscp1 = size(Trj.RA, 2);
                    nCvr = size(Trj.Cvr, 2);
                    
                    for iA = 1:nAscp1
                        % Plot without covariate
                        tFilNam = sprintf('%s_%s_Bin%g_Time', FilNam, RspLbl{iA}, iB);
                        tRsp = sprintf('%s:%s', RspLbl{iA}, Bn.BinLbl{iB});
                        obj.plotStormTrajectories(Y(Bn.A==iB,1), Trj.RA(Bn.A==iB, iA), [], tFilNam, tRsp, '');
                        
                        for iC = 1:nCvr
                            tFilNam = sprintf('%s_%s_Bin%g_%s', FilNam, RspLbl{iA}, iB, CvrLbl{iC});
                            tRsp = sprintf('%s:%s', RspLbl{iA}, Bn.BinLbl{iB});
                            obj.plotStormTrajectories(Y(Bn.A==iB,1), Trj.RA(Bn.A==iB, iA), Trj.Cvr(Bn.A==iB, iC), tFilNam, tRsp, CvrLbl{iC});
                        end
                    end
                end
            end
        end   
        function obj=plotStormTrajectories(obj,Y, TrjRA, TrjCvr, tFilNam, LblRA, LblCvr)
            nPk = size(Y, 1);
            [peakOrder, colorMap] = obj.getPeakOrderAndColorMap(Y(:, 1), nPk);
            
            clf;
            for isNormalised = [true, false]
                subplot(2, 1, 1 + ~isNormalised)
                hold on;
                for i = 1:nPk
                    peakIndex = peakOrder(i);
                    stormRsp = TrjRA{peakIndex, 1};
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
                    plot(xData, yData, 'Color', colorMap(i, :), 'LineWidth', 0.5, 'LineStyle', ':', 'Marker', '.', 'MarkerSize', 10);
                end
                
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
            end
            savePics(fullfile('Figures', tFilNam));
        end %plotStormTrajectories
        
    end %methods
    
    methods (Static)
        function [peakOrder, colorMap] = getPeakOrderAndColorMap(peakValues, nPk)
            [~, peakOrder] = sort(peakValues, 'descend');
            colorMap = jet(nPk);
            colorMap = colorMap(end:-1:1, :);
        end %getPeakOrderAndColorMap
        
        function CntStrCvr = adjustCovariate(CntStrCvr)
            CntStrCvr(CntStrCvr < -180) = CntStrCvr(CntStrCvr < -180) + 360;
            CntStrCvr(CntStrCvr > 180) = CntStrCvr(CntStrCvr > 180) - 360;
        end %adjustCovariate
        
        function [titleStr, xlabelStr, ylabelStr, xlimits] = getPlotProperties(isNormalised, TrjCvr, LblRA, LblCvr)
            if isNormalised
                titleStr = sprintf('%s: Normalised Storm Trajectories', LblRA);
                ylabelStr = 'Normalised Response (Max=1)';
            else
                titleStr = 'Unnormalised Storm Trajectories';
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