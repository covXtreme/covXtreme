clear; clf;
addpath('..\Code');

set(0,'defaultFigureColor',[1,1,1]);
%set(0, 'DefaultAxesFontName', fontname)
set(0, 'DefaultAxesBox', 'on');
set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on');
%% Stage4_FitH&T
%%%%%%%%%%%%%%%%%%%%%%%%%
FM=cellstr(ls('Output/MM*.mat'));
for iM=1:numel(FM) %load all marginal models
    load(sprintf('Output/MM%g',iM),'MM')  %Load margin 2
    if iM==1
        Mrg=MM;
    else
        Mrg=cat(1,Mrg,MM);
    end
    clear MM;
end
load('Output/HT','HT')

% TODO current setup is the only combination that does not work 
% RVMethod='retValX';
% A=[1,1,2,2];
% nA=numel(unique(A));1
% inputs
nRls=1000; % number of realisations
A=[1:4];% vector of bins
nA=numel(unique(A)); % number of bins

% selection of the different return value methods 
%RVMethod='randX';
RVMethod='retValX';
HT.SmpLclRsdOn=true;
% user defined X
switch RVMethod
    case 'userX' 
        RVX=4*ones(nA,HT.nBoot,HT.nRtr);
        Sml=SimulateMC(HT,Mrg,'userX',nRls,A,nA,RVX);                
    case 'retValX' % X is generated from the return value model
        Sml=SimulateMC(HT,Mrg,'retValX',nRls,A,nA);
                
        %% Check Marginal distributions
        for iMrg=1:2
            figure(2+iMrg);
            clf;
            iRtr=2;
            
            for iA=1:nA
                subplot(nA,1,iA)
                if nA==Mrg(iMrg).Bn.nBin
                    plot(Mrg(iMrg).RVX,Mrg(iMrg).RVPrb(iA,:,iRtr),'k-')
                elseif nA==1 %omni
                    plot(Mrg(iMrg).RVX,Mrg(iMrg).RVPrb(end,:,iRtr),'k-')
                end
                hold on
                plot(sort(Sml.Org(iA,:,iMrg,iRtr)),linspace(0,1,nRls),'r--')
                xlim([min(Sml.Org(iA,:,iMrg,iRtr)),max(Sml.Org(iA,:,iMrg,iRtr))])
                xlabel('Prob');
                ylabel(Mrg(1).RspLbl);
            end
        end         
    case 'randX'
        Sml=SimulateMC(HT,Mrg,'randX',nRls,A,nA);
                
         %% Check Marginal distributions
        for iMrg=1:2
            figure(2+iMrg);
            clf;            
            P=CDF(Mrg(iMrg),shiftdim(Mrg(iMrg).RVX,-2)); %mean over bootstraps
            for iA=1:nA
                subplot(nA,1,iA)
                if nA==Mrg(iMrg).Bn.nBin
                    plot(Mrg(iMrg).RVX,log10(1-squeeze(mean(P(iA,:,:),2))),'k-')
                elseif nA==1 %omni
                    plot(Mrg(iMrg).RVX,log10(1-squeeze(mean(prod(P(iA,:,:),1),2))),'k-')
                end
                hold on
                plot(sort(Sml.Org(iA,:,iMrg)),log10(1-linspace(0,1,nRls)),'r--')
                xlim([min(Sml.Org(iA,:,iMrg)),max(Sml.Org(iA,:,iMrg))])
                ylabel('log10(1-Prob)');
                xlabel(Mrg(1).RspLbl);
              
            end
            legend('CDF','MC')
            if iMrg==1
                sgtitle('Main')
            else
                sgtitle('Associated')
            end
        end 
end

% color scrambling
% figure showing scatter plots of original margins
iRtr=1;

figure(1)
clf;
iPlt=0;
for iA=1:nA
    iPlt=iPlt+1;
    subplot(2,nA,iPlt)
    %plot(Sml.Org(iA,:,1,iRtr),Sml.Org(iA,:,2,iRtr),'.k');
    scatter(Sml.Org(iA,:,1,iRtr),Sml.Org(iA,:,2,iRtr),14,Sml.A(iA,:,iRtr),'filled');
    caxis([1,HT.nBin])
    xlabel(Mrg(1).RspLbl);
    ylabel(Mrg(2).RspLbl);
    title(sprintf('Bin %d - original margins',iA));
end
iPlt=nA;

for iA=1:nA
    iPlt=iPlt+1;
    subplot(2,nA,iPlt)
    scatter(Sml.StnMrg(iA,:,1,iRtr),Sml.StnMrg(iA,:,2,iRtr),14,Sml.A(iA,:,iRtr),'filled');
    caxis([1,HT.nBin])
    xlabel(Mrg(1).RspLbl);
    ylabel(Mrg(2).RspLbl);
    title(sprintf('Bin %d - standard margins',iA));
end

figure(2)
clf;
for iMrg=1:2
    subplot(2,1,iMrg)
    scatter(Sml.X(:),reshape(Sml.Org(:,:,2,:),[],1),20,Sml.A(:),'filled')
    xlabel(Mrg(iMrg).CvrLbl);
    ylabel(Mrg(iMrg).RspLbl);
end
title(sprintf('X is the %s method', RVMethod));