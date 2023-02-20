function RVDgnPPC(Mrg,RV,SavFld);
%function RVDgnPPC(Mrg,RV,SavFld);
%
%S.Sam, P.Jonathan 20170605
%Modified from E.Ross code
%
%1. Plot marginal and conditional return values from PPC analysis
%2. Box-whisker plots

%% 1. Plot marginal and conditional return values from PPC analysis
%Mrg has Mrg(1) (condititioning) and Mrg(2) (conditioned)
%Mrg, HT and RV created by stages 1-4 of PPC
%Fld is string to relative folder location

clf;

nRtrPrd=RV.nRtrPrd;
for iRP=1:nRtrPrd;
    
    subplot(iRP,2,1+2*(nRtrPrd-1)); hold on;
    tPrb=Mrg(1).RVPrb(end,:)';
    plot(Mrg(1).RVX,tPrb,'k');
    tPrb=(1:size(RV.X,2))/(size(RV.X,2)+1);
    plot(sort(permute(RV.X(iRP,:,1),[2 3 1])),tPrb,'k--');
    box on; grid on;
    ylim([0.01,0.99]);
    ylabel 'Probability';
    xlabel(Mrg(1).Lbl);
    title(sprintf('%g-year',RV.RtrPrd(iRP)));
    legend('MrgFromMrg','MrgFromHT','location','SE');
    
    subplot(iRP,2,2+2*(nRtrPrd-1)); hold on;
    tPrb=Mrg(2).RVPrb(end,:)';
    plot(Mrg(2).RVX,tPrb,'k');
    tPrb=(1:size(RV.X,2))/(size(RV.X,2)+1);
    plot(sort(permute(RV.X(iRP,:,3),[2 3 1])),tPrb,'k--');
    plot(sort(permute(RV.X(iRP,:,2),[2 3 1])),tPrb,'r');
    box on; grid on;
    ylim([0.01,0.99]);
    ylabel 'Probability'
    xlabel(Mrg(2).Lbl);
    legend('MrgFromMrg','MrgFromHT','Cnd','location','SE');
    
end;

savePics(fullfile(SavFld,'Stg5_RtrValCdfDgn'),10,6,'jpg')

%% 2. Box plots broken out into 2 rows: Hs for different RtnPrds and Srg,Srg|Hs on same plot

%prep data
BoxWdt = 0.1;%box plot width
SrgLbl = repmat(RV.Lbl(3:-1:2),RV.nSim,1); %index for boxplots for diff return periods
tRVSrgHs =RV.X(:,:,2)';
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
    %
    %Use estimate of marginal Hs return value distribution from numerical
    %integration to create box plot
    [~,tUnq]=unique(Mrg(1).RVPrb(end,:)'); %unique values of probability
    t1=interp1(Mrg(1).RVPrb(end,tUnq)',Mrg(1).RVX(tUnq),[0.025,0.25,0.5,0.75,0.975]');
    %
    NewBoxPlt(t1,2); %our own boxplot function
    %
    if iRtrPrd ==1
        title(sprintf('%d yr (nSim = %d)',RV.RtrPrd(iRtrPrd),RV.nSim))
        ylabel('Hs')
    else
        title(sprintf('%d yr',RV.RtrPrd(iRtrPrd)))
    end
    set(gca,'xtick',[])
    %
    %Set sensible ylim
    ylim([t1(1)-range(t1)*0.05,t1(end)+range(t1)*0.05])
    
    % Srg subplots -------------------------------------
    subplot(nPltRw,nPltCol,iRtrPrd+RV.nRtrPrd)
    %
    %Use estimate of marginal Srg return value distribution from numerical
    %integration to create box plot
    [~,tUnq]=unique(Mrg(2).RVPrb(end,:)'); %unique values of probability
    t1=interp1(Mrg(2).RVPrb(end,tUnq)',Mrg(2).RVX(tUnq),[0.025,0.25,0.5,0.75,0.975]');
    %
    %Add empirical estimates for conditional surge distribution from
    %simulation performed
    t2=quantile(tRVSrgHs(:,iRtrPrd),[0.025,0.25,0.5,0.75,0.975],1);
    %
    NewBoxPlt([t1 t2]); %our own boxplot function
    set(gca,'XTickLabel',SrgLbl(1,:))
    %
    if iRtrPrd ==1
        ylabel('Surge')
    end
    %
    %Set sensible ylim
    t3=[min([t1;t2]);max([t1;t2])];
    ylim([t3(1)-range(t3)*0.05,t3(end)+range(t3)*0.05])
    
end    %loop over return periods

savePics(fullfile(SavFld,'Stg5_ReturnValBoxplotsNew'),10,6,'jpg')

return;