%% Check data look like Shiraz's
clear;clc;
Lct='Clipper';
AnlType='Mxm';
Datpath=sprintf('R:/Emma.Ross/Extremes/CE_Surge/Analysis/%s/%s/Output',Lct,AnlType);
load(fullfile(Datpath,'Data.mat'))
load(fullfile(Datpath,'Bin.mat'))
BinAlc=BinAllocation(Dat.X,DrcEdg,Dat.Y,false); % Allocate Data to Bins and Plot
nBin=numel(unique(BinAlc));
Bins = unique(BinAlc);
figure(1);clf;
for iBin=1:nBin
    subplot(2,2,iBin)
    plot(Dat.Y(BinAlc==Bins(iBin),1),Dat.Y(BinAlc==Bins(iBin),2),'k.')
    ylabel(AnlType)
    %mod(iBin+1,nBin)
    xlabel(sprintf('Bin %d: [%d,-]',Bins(iBin),DrcEdg(Bins(iBin))))
end

figure(2);clf;
subplot(2,1,1)
plot(Dat.X,Dat.Y(:,1),'k.')
xlabel('Direction')
ylabel(Dat.Lbl{1})
subplot(2,1,2)
plot(Dat.X,Dat.Y(:,2),'k.')
xlabel('Direction')
ylabel(Dat.Lbl{2})



%% Boxplots to compare against Shiraz's analysis
% (stationary H&T)
% First nBin columns of HT.Prm are alphas for each bin. Remaining cols are 
%... Beta, Mu, Sigma (which don't vary by bin). 
HTpath=sprintf('R:/Emma.Ross/Extremes/CE_Surge/Analysis/%s/%s/Output/HT.mat',Lct,AnlType);
load(HTpath)

%alpha parameter
if HT.NonStat
    nAlpha = HT.nBin;
else
   nAlpha = 1; 
   iBin=1;
end
alphaHat=HT.Prm(1:nAlpha,:); %rows = alphaHat for each bin, cols = bootstrap
betaHat=HT.Prm(nAlpha+1,:);
muHat=HT.Prm(nAlpha+2,:);
sigmaHat=HT.Prm(nAlpha+3,:);

figure(1);clf;
boxplot(alphaHat(iBin,:))
title(Lct)
xlabel('\alpha')
ylabel(AnlType)
ylim([0,1])

figure(3);
boxplot(betaHat)
title(Lct)
xlabel('\beta')
ylabel(AnlType)
figure(4);
boxplot(muHat)
title(Lct)
xlabel('\mu')
ylabel(AnlType)
figure(5);
boxplot(sigmaHat)
title(Lct)
xlabel('\sigma')
ylabel(AnlType)