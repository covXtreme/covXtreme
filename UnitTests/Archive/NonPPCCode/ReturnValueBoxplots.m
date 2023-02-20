function ReturnValuePlotsPPC(RV,PltLayout,SvFld)
% Plots box plots and empirical CDFs of return values from Mrg and HT PPC analysis
% Inputs:
% - RV = data structure containing return values
%    * RV.X = [nRtrPrd x nSim x 3] matrix of return values where cols =
%    Hs,Srg|Hs,Srg
% - PltLayout = 0 or 1 indicating how to display results:
%     * 1 --> 3 subplots (Hs,Srg|Hs,Srg) multi-boxplots for return period 
%     * 0 --> 2 subplots (top row = Hs, bottom row = Srg|Hs and Srg combined, return periods in cols
% - SvFld = folder in which will save plots
% Output:
% Figure 1 = boxplot return values for 2 margins and Mrg2|Mrg1.
% Figure 2 = return value CDF 
%% 
if (PltLayout > 1) || (PltLayout < 0)
    error('PltLayout option must be 0 or 1')
end
%% Box plots broken out into 2 rows: Hs for different RtnPrds and Srg,Srg|Hs on same plot

%prep data
BoxWdt = 0.1;%box plot width
SrgLbl = repmat(RV.Lbl(3:-1:2),RV.nSim,1); %index for boxplots for diff return periods
SrgLblI =  repmat(1:2,RV.nSim,1);
tRVHs =RV.X(:,:,1)';
tRVSrgHs =RV.X(:,:,2)';
tRVSrg =RV.X(:,:,3)';
%Samta
HsLim = [min(tRVHs(:))-1,max(tRVHs(:))+1];
%SrgLim = [min([tRVSrgHs(:);tRVSrg(:)])-0.1,max([tRVSrgHs(:);tRVSrg(:)])+0.1];
%XLimHs = [(1-BoxWdt*2.5),(1+BoxWdt*2.5)];
XLimSrg = [(1-BoxWdt*5),(2+BoxWdt*5)];
%subplot layout
if RV.nRtrPrd == 1
    nPltRw = 1;
    nPltCol = 2;
else
    nPltRw = 2;
    nPltCol = RV.nRtrPrd;
end
%figure(1);clf
for iRtrPrd = 1:RV.nRtrPrd
    % Srg subplots -------------------------------------
    tRVSrgCmb =[tRVSrg(:); tRVSrgHs(:)]; %combine surge results
    
    if 0 %Boxplot style (matlab or ours)
        boxplot([tRVSrg(:,iRtrPrd); tRVSrgHs(:,iRtrPrd)],SrgLbl(:)); %matlab
    else
        BoxPlt([tRVSrg(:,iRtrPrd); tRVSrgHs(:,iRtrPrd)],SrgLblI(:),BoxWdt); %our own boxplot function
        %set(gca,'XTickLabel',SrgLbl(1,:))
       
        xlim(XLimSrg)
    end
    if iRtrPrd ==1
        %ylabel('Surge')
    end
    %ylim(SrgLim)
    
end    %loop over return periods

savePics(fullfile(SvFld,'Stg5_ReturnValBoxplots'),10,6,'jpg')

end