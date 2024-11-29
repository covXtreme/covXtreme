% Copyright © [2023] Shell Global Solutions International B.V. All Rights Reserved.
% SPDX-License-Identifier: Apache-2.0

classdef StormTrajectorySimulation
    properties
    end %properties
    
    methods
        function obj=StormTrajectorySimulation(Dat,Mrg, HT)
            
            obj.IdentifyStormTrajectories(Dat);
        end %StormTrajectorySimulation
        
        function obj=IdentifyStormTrajectories(obj,Dat)
            Rsp=Dat.Y(:,1);
            Asc=Dat.Y(:,2:end);
            Cvr=Dat.X;
            nExc=size(Dat.Y,1);
            
            obj.RA=cell(nExc,1+nAsc); %Initialise empty cell array for storm trajectories
            obj.Cvr=cell(nExc, nCvr); %Initialise empty cell array for covariate trajectories.
            
            %% For each storm peak, get its trajectory.
            for i = 1:nExc
                % Find the start and end of each storm
                startIdx = Prd(i, 1); % start of the storm
                endIdx = Prd(i, 2);   % end of the storm
                
                % Extract the storm trajectory within the storm period
                stormRsp = Rsp(startIdx:endIdx);  % Response trajectory
                stormAsc = Asc(startIdx:endIdx, :);  % Associated variables trajectory
                stormCvr = Cvr(startIdx:endIdx, :); % Covariates
                
                % Store the storm trajectory
                obj.RA{i, 1} = stormRsp; % Store response trajectory
                for iA = 1:nAsc
                    obj.RA{i, iA+1} = stormAsc(:, iA);  % Store associated variables trajectory
                    obj.Cvr{i, iA} = stormCvr(:, iA);
                end
            end
            
            TrajectoryPlot(Dat.Y, obj, Dat.RspLbl, Dat.CvrLbl, [], 'Stg6_Data_StormTrajectory')
            
            
        end %IdentifyStormTrajectories
        
        function obj=CreateStormTrajectoryBinAllocation(obj, StrTrj)
            fprintf(1,'Calculating bin allocations for historical trajectories\n');
            for iS=1:Bn.n
                tCvr=zeros(Bn.n,1);
                tStrTrjCvr=Dat.StrTrj.Cvr{iS,:};
                tCvr(1:size(tStrTrjCvr,1))=tStrTrjCvr;
                iVrb = 0;
                tBn=BinAllocation(Bn,tCvr, iVrb);
                StrTrj.A{iS,:}=tBn.A(1:size(tStrTrjCvr,1));
            end %iS
        end %CreateStormTrajectoryBinAllocation
        
        function obj=StormMatching(obj,StrTrj)
            %% Storm matching to allocate historical trajectories to simulated storms
            % Match using the dominant variate only, and in the covariate bin
            % corresponding to the dominant variate only.
            nNgh=10; %Number of near neighbour historical storm trajectories to consider
            fprintf(1,'Storm matching to allocate historical trajectories to simulated storms:\n');
            iE=0;
            for iS=1:nSml
                for iR=1:Sml{iS}.nRls
                    
                    iE=iE+1;
                    
                    tBn=Sml{iS}.A(iR); %bin for realisation iR from simulation iS
                    
                    tMB=find(Bn.A==tBn); %all historical matched with correct bin allocation
                    
                    [~,tD2]=sort((Dat.Y(tMB,1)-Sml{iS}.Org(iR)).^2);
                    
                    tMV=tMB(tD2(1:nNgh));
                    
                    tMtc=tMV(randi(nNgh,1));
                    
                    Sml{iS}.StrTrj.RA{iR,:}=StrTrj.RA{tMtc,:};
                    Sml{iS}.StrTrj.Cvr{iR,:}=StrTrj.Cvr{tMtc,:};
                    Sml{iS}.StrTrj.A{iR,:}=StrTrj.A{tMtc,:};
                    Dgn(iE,:)=[max(StrTrj.RA{tMtc,1}) Sml{iS}.Org(iR) Sml{iS}.A(iR)];
                    
                end; %iR
            end; %iS
            
        end %StormMatching
        
        function TrajectoryPlot(Y, Trj, RspLbl, CvrLbl, Bn, FilNam)
            
            if isempty(Bn)
                nAscp1 = size(Trj.RA, 2);
                nCvr = size(Trj.Cvr, 2);
                
                for iA = 1:nAscp1
                    % Plot without covariate
                    tFilNam = sprintf('%s_%s_Time', FilNam, RspLbl{iA});
                    plotStormTrajectories(Y, Trj.RA(:, iA), [], tFilNam, RspLbl{iA}, '');
                    
                    for iC = 1:nCvr
                        tFilNam = sprintf('%s_%s_%s', FilNam, RspLbl{iA}, CvrLbl{iC});
                        plotStormTrajectories(Y, Trj.RA(:, iA), Trj.Cvr(:, iC), tFilNam, RspLbl{iA}, CvrLbl{iC});
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
                        plotStormTrajectories(Y(Bn.A==iB,1), Trj.RA(Bn.A==iB, iA), [], tFilNam, tRsp, '');
                        
                        for iC = 1:nCvr
                            tFilNam = sprintf('%s_%s_Bin%g_%s', FilNam, RspLbl{iA}, iB, CvrLbl{iC});
                            tRsp = sprintf('%s:%s', RspLbl{iA}, Bn.BinLbl{iB});
                            plotStormTrajectories(Y(Bn.A==iB,1), Trj.RA(Bn.A==iB, iA), Trj.Cvr(Bn.A==iB, iC), tFilNam, tRsp, CvrLbl{iC});
                        end
                    end
                end
            end
        end
        
        function plotStormTrajectories(Y, TrjRA, TrjCvr, tFilNam, LblRA, LblCvr)
            nPk = size(Y, 1);
            [peakOrder, colorMap] = getPeakOrderAndColorMap(Y(:, 1), nPk);
            
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
                        CntStrCvr = adjustCovariate(stormCvr - peakCvr);
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
                [titleStr, xlabelStr, ylabelStr, xlimits] = getPlotProperties(isNormalised, TrjCvr, LblRA, LblCvr);
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
        
    end %methods
    
end %classdef