function TrajectoryPlot2(Y, Trj, RspLbl, CvrLbl, Bn, FilNam)

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
end

function [peakOrder, colorMap] = getPeakOrderAndColorMap(peakValues, nPk)
    [~, peakOrder] = sort(peakValues, 'descend');
    colorMap = jet(nPk);
    colorMap = colorMap(end:-1:1, :);
end

function CntStrCvr = adjustCovariate(CntStrCvr)
    CntStrCvr(CntStrCvr < -180) = CntStrCvr(CntStrCvr < -180) + 360;
    CntStrCvr(CntStrCvr > 180) = CntStrCvr(CntStrCvr > 180) - 360;
end


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
end
