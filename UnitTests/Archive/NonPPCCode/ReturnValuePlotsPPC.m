function ReturnValuePlotsPPC(RV,PltLayout,SvFld)
% Plots box plots and empirical CDFs of return values from Mrg and HT PPC analysis
% Inputs:
% - RV = data structure containing return values
%    * RV.X = [nRtrPrd x nSim x 3] matrix of return values where cols =
%    Hs,Srg|Hs,Srg
% - PltLayout = 0 or 1 indicating how to display results:
%     * 1 --> 3 subplots (Hs,Srg|Hs,Srg) multi-boxplots for return period 
%     * 0 --> 6 subplots (top row = Hs, bottom row = Srg|Hs and Srg combined, return periods in cols
% - SvFld = folder in which will save plots
% Output:
% Figure 1 = boxplot return values for 2 margins and Mrg2|Mrg1.
% Figure 2 = return value CDF 
%% 
if (PltLayout > 1) || (PltLayout < 0)
    error('PltLayout option must be 0 or 1')
end
if PltLayout==1 
    %% Box plots broken out by response, multi plots within each window for return period
    IndBox =  ones(RV.nSim,1)*RV.RtrPrd; %index for boxplots for diff return periods    
    
    figure(1);clf
    subplot(1,3,1)
    tRVHs =RV.X(:,:,1)';
    boxplot(tRVHs(:),IndBox(:))
    title(sprintf('%s (nSim = %d)',RV.Lbl{1},RV.nSim))
    xlabel('Return Period')
    ylabel('Hs')
    
    subplot(1,3,2)
    tRVSrgHs =RV.X(:,:,2)';
    tRVSrg =RV.X(:,:,3)';
    SrgLim = [min([tRVSrgHs(:);tRVSrg(:)])-0.1,max([tRVSrgHs(:);tRVSrg(:)])+0.1];
    boxplot(tRVSrgHs(:),IndBox(:))
    title(RV.Lbl{2})
    xlabel('Return Period')
    ylabel('SrgMdn')
    ylim(SrgLim)
    
    subplot(1,3,3)
    boxplot(tRVSrg(:),IndBox(:))
    title(RV.Lbl{3})
    xlabel('Return Period')
    ylabel('SrgMdn')
    ylim(SrgLim)
else
    %% Box plots broken out into 2 rows: Hs for different RtnPrds and Srg,Srg|Hs on same plot

    %prep data
    BoxWdt = 0.1;%box plot width 
    SrgLbl = repmat(RV.Lbl(3:-1:2),RV.nSim,1); %index for boxplots for diff return periods
    SrgLblI =  repmat(1:2,RV.nSim,1);
    tRVHs =RV.X(:,:,1)';
    tRVSrgHs =RV.X(:,:,2)';
    tRVSrg =RV.X(:,:,3)';
    HsLim = [min(tRVHs(:))-1,max(tRVHs(:))+1];
    SrgLim = [min([tRVSrgHs(:);tRVSrg(:)])-0.1,max([tRVSrgHs(:);tRVSrg(:)])+0.1];
    XLimHs = [(1-BoxWdt*2.5),(1+BoxWdt*2.5)];
    XLimSrg = [(1-BoxWdt*5),(2+BoxWdt*5)];
    %subplot layout
    if RV.nRtrPrd == 1
        nPltRw = 1;
        nPltCol = 2;
    else
        nPltRw = 2;
        nPltCol = RV.nRtrPrd;
    end
    figure(1);clf   
    for iRtrPrd = 1:RV.nRtrPrd
        % Hs subplots -----------------------------------------
        subplot(nPltRw,nPltCol,iRtrPrd)
        if 0 %Boxplot style (matlab or ours)
            boxplot(tRVHs(:,iRtrPrd)) %matlab
        else
            BoxPlt(tRVHs(:,iRtrPrd),ones(RV.nSim,1),BoxWdt); %our own boxplot function
            xlim(XLimHs)
        end
        if iRtrPrd ==1
            title(sprintf('%d yr (nSim = %d)',RV.RtrPrd(iRtrPrd),RV.nSim))
            ylabel('Hs')
        else
            title(sprintf('%d yr',RV.RtrPrd(iRtrPrd)))
        end
        set(gca,'xtick',[])
        ylim(HsLim)
        
%         %Make subplots closer together (top row)
%         if RV.nRtrPrd > 1
%             pos = get(gca, 'Position');  %[x y width height]
%             pos(1)=pos(1)-0.05; %shift down and left
%             pos(2)=pos(2)-0.05;
%             pos(3) = 1.2*pos(3); %inflate
%             pos(4) = 1.2*pos(4);
%             set(gca, 'Position', pos);
%         end
        
        % Srg subplots -------------------------------------
        subplot(nPltRw,nPltCol,iRtrPrd+RV.nRtrPrd)
        tRVSrgCmb =[tRVSrg(:); tRVSrgHs(:)]; %combine surge results
        
        if 0 %Boxplot style (matlab or ours)
            boxplot([tRVSrg(:,iRtrPrd); tRVSrgHs(:,iRtrPrd)],SrgLbl(:)); %matlab
        else
            BoxPlt([tRVSrg(:,iRtrPrd); tRVSrgHs(:,iRtrPrd)],SrgLblI(:),BoxWdt); %our own boxplot function
            set(gca,'XTickLabel',SrgLbl(1,:))      
            xlim(XLimSrg)
        end
        if iRtrPrd ==1
            ylabel('Surge')
        end
        ylim(SrgLim)
        
        
%         %Make subplots closer together (bottom row)
%         if RV.nRtrPrd > 1
%             pos = get(gca, 'Position');  %[x y width height]
%             pos(1)=pos(1)-0.05; %shift down and left
%             pos(2)=pos(2)-0.05;
%             pos(3) = 1.1*pos(3); %inflate
%             pos(4) = 1.1*pos(4);
%             set(gca, 'Position', pos);
%         end
    end    %loop over return periods
    
end %plot style
savePics(fullfile(SvFld,'Stg5_ReturnValBoxplots'),10,6,'jpg')

%% Return Value CDF Plot
%X values -----
tAsc = RV.X(:,:,2:3);
nX = 1000; %number of x values at which cdf is calculated
AscRspX = linspace(min(tAsc(:))-0.1,max(tAsc(:))+0.1,nX)'; %srg
tMain =  RV.X(:,:,1);
MainRspX = linspace(min(tMain(:)),max(tMain(:)),nX)'; %Hs
RV.CdfX = [MainRspX,AscRspX,AscRspX];%X vals Hs, Srg|Hs, Srg

%Prb vales -----
RV.CdfPrb = NaN(RV.nRtrPrd,nX,3);

%Empirical CDF Y2|Y1 large
for iRP = 1:RV.nRtrPrd %loop over return period
    tRvx = RV.X(iRP,:,:);
    for iX = 1:nX %loop over x
        I= bsxfun(@lt,squeeze(tRvx),RV.CdfX(iX,:));
        RV.CdfPrb(iRP,iX,:)=(1/RV.nSim).*sum(I,1);
    end
end
RVLgnd = cellstr(num2str(RV.RtrPrd', '%-dYr'));
tC = get(gca,'ColorOrder');
PltCol = tC(1:3,:);
%Figure -------
figure(2);clf;
subplot(2,1,1) %Hs
for iRP = 1:RV.nRtrPrd
    plot(RV.CdfX(:,1),RV.CdfPrb(iRP,:,1),'linewidth',2,'Color',PltCol(iRP,:))
    hold on;
end
ylim([0.01,0.99])
ylabel('Cumulative Probability')
xlabel(RV.Lbl(1))
legend(RVLgnd,'location','best');

subplot(2,1,2) %Srg
for iRP = 1:RV.nRtrPrd
    plot(RV.CdfX(:,2),RV.CdfPrb(iRP,:,2),'-.','linewidth',2,'Color',PltCol(iRP,:))
    hold on;
    plot(RV.CdfX(:,3),RV.CdfPrb(iRP,:,3),'-','linewidth',2,'Color',PltCol(iRP,:))
end
hold off
ylim([0.01,0.99])
ylabel('Cumulative Probability')
xlabel(RV.Lbl(3))
legend(cellstr(RV.Lbl(2:3)),'location','best');
savePics(fullfile(SvFld,'Stg5_ReturnValCdf'),10,6,'jpg')


end