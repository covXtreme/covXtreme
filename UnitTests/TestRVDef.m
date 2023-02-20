clc; clear;


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
%load Heffernan and Tawn model
load('Output/HT','HT');

%% Definition 1 
%(current approach choose x value e.g. 10k Mrg RV look at distribution of
%Y| X=RV 10k
nG=40;
Edg=linspace(1,8,nG+1);
Md=(Edg(2:end)+Edg(1:end-1))./2;

qLv=[0.05,0.5,0.95];
QY=HT.ConditionalQuantile(Md',ones(HT.nBin,1),1,qLv,false);


%% N Yr RV simualtion 
Mrg(1).RtrPrd=logspace(0,8,9);
Mrg(2).RtrPrd=Mrg(1).RtrPrd;
Mrg(1).nRtr=numel(Mrg(1).RtrPrd);
nRls=10000;
HT2=HT.ConditionalReturnValue(Mrg,nRls);
RV1=HT2.RV;

CndRV{1}.RVX=quantile(squeeze(RV1.X(1,:,:)),exp(-1));
CndRV{1}.RVY=quantile(squeeze(RV1.Y(1,1,:,:)),0.5);

CndRV{3}.RVX=quantile(squeeze(RV1.X(1,:,:)),0.5);
CndRV{3}.RVY=quantile(squeeze(RV1.Y(1,1,:,:)),0.5);

CndRV{7}.RVX=mean(squeeze(RV1.X(1,:,:)));
CndRV{7}.RVY=mean(squeeze(RV1.Y(1,1,:,:)));
%% define RV on each bootstrap then get condiiotnal median for each then summarise
%get PSSRV return level --> X per bootstrap
nRls=1000;
nX=1000;
X=linspace(0,40,nX)';
P=Mrg(1).CDF(X); %nX x nB
% take quantile of R

tHT=HT;
tHT.RVMth=3;
CndRV{2}.RVX=NaN(Mrg(1).nRtr,Mrg(1).nBoot);
CndRV{2}.RVY=NaN(Mrg(1).nRtr,Mrg(1).nBoot);
CndRV{4}.RVX=NaN(Mrg(1).nRtr,Mrg(1).nBoot);
CndRV{4}.RVY=NaN(Mrg(1).nRtr,Mrg(1).nBoot);
for iB=1:Mrg(1).nBoot
    %% T -year
    for iRtr=1:Mrg(1).nRtr        
        alp=exp(-1);
        R=@(x)alp-exp(-Mrg(1).Rat(iB).*Mrg(1).RtrPrd(iRtr).*(1-Mrg(1).CDF(x,1,iB))); %nX x nB
        x0=Mrg(1).Thr(iB);
        CndRV{2}.SimRVX(iRtr,iB)=fzero(R,x0);
%         if     CndRV{2}.SimRVX(iRtr,iB)<3
%            'here'; 
%         end
    end
    tHT.RVFix=CndRV{2}.SimRVX(:,iB)';
    I=ones(nRls,1).*iB;
    tHT=tHT.ConditionalReturnValue(Mrg,nRls,I);
    
    CndRV{2}.SimRVY(:,iB)=quantile(squeeze(tHT.RV.Y(1,1,:,:)),0.5,1);
    
    %% Annual
    for iRtr=1:Mrg(1).nRtr
        alp=max(1-1./Mrg(1).RtrPrd(iRtr),exp(-1));
        R=@(x)log(alp)+Mrg(1).Rat(iB).*(1-Mrg(1).CDF(x,1,iB)); %nX x nB
        x0=Mrg(1).Thr(iB);
        CndRV{4}.SimRVX(iRtr,iB)=fzero(R,x0);
    end
    tHT.RVFix=CndRV{4}.SimRVX(:,iB)';
    I=ones(nRls,1).*iB;
    tHT=tHT.ConditionalReturnValue(Mrg,nRls,I);
    
    CndRV{4}.SimRVY(:,iB)=quantile(squeeze(tHT.RV.Y(1,1,:,:)),0.5,1);
    
end
CndRV{2}.RVX=mean(CndRV{2}.SimRVX,2);
CndRV{2}.RVY=mean(CndRV{2}.SimRVY,2);
CndRV{4}.RVX=mean(CndRV{4}.SimRVX,2);
CndRV{4}.RVY=mean(CndRV{4}.SimRVY,2);
CndRV{5}.RVX=median(CndRV{2}.SimRVX,2);
CndRV{5}.RVY=median(CndRV{2}.SimRVY,2);
CndRV{6}.RVX=median(CndRV{4}.SimRVX,2);
CndRV{6}.RVY=median(CndRV{4}.SimRVY,2);
%%

iRtr=2;
figure(1);
clf;
hold on
if 1
load('Output/DataLarge');
plot(Dat.Y(:,1),Dat.Y(:,2),'.','color',[1,1,1]*0.5,'handlevisibility','off')
end
plot(Mrg(1).Y,Mrg(2).Y,'k.','markersize',10,'handlevisibility','off')
hold on

C=lines(7);
plot(CndRV{1}.RVX,CndRV{1}.RVY,'-','color','r','linewidth',2)
plot(CndRV{3}.RVX,CndRV{3}.RVY,'-','color','b','linewidth',2)
%plot(CndRV{2}.RVX(end,:),CndRV{2}.RVY(end,:),'b.')
plot(CndRV{2}.RVX,CndRV{2}.RVY,'-','color','g','linewidth',2)
plot(CndRV{4}.RVX,CndRV{4}.RVY,'--','color','c','linewidth',2)
plot(CndRV{5}.RVX,CndRV{5}.RVY,':','color','m','linewidth',2)
plot(CndRV{6}.RVX,CndRV{6}.RVY,'-.','color','y','linewidth',2)
plot(CndRV{7}.RVX,CndRV{7}.RVY,':','color','b','linewidth',2)
hold on
plot(Md,QY,'k--','linewidth',1.5)
grid on
%set(gca,'xscale','log','yscale','log')
legend({'QuantileMCX 0.37and Y 0.5',
'MedianMC',
'MeanQuantile T_Year',
'MeanQuantile 1-year',
'MedianQuantile T_Year',
'MedianQuantile 1-year',
'MeanMC',
'QuantileMCYgXRV'},'location','SouthEast')

%%
if 1
    %iRtr=4;
    figure(2);
    clf;
    for iRtr=1:numel(Mrg(1).RtrPrd)
        subplot(3,3,iRtr);
        hold on
        if 0
            load('Output/DataLarge');
            plot(Dat.Y(:,1),Dat.Y(:,2),'.','color',[1,1,1]*0.5,'handlevisibility','off')
        end
        plot(Mrg(1).Y,Mrg(2).Y,'k.','markersize',10,'handlevisibility','off')
        hold on
        plot(RV1.X(1,:,iRtr)',squeeze(RV1.Y(1,1,:,iRtr)),'r.')
        U=[0.25,0.25];
        nG=40;
        [ptsX,ptsY]=meshgrid(linspace(0,30,nG),linspace(0,30,nG));
        pts=[ptsX(:),ptsY(:)];
        f=ksdensity([RV1.X(1,:,iRtr)',squeeze(RV1.Y(1,1,:,iRtr))],pts,'Bandwidth',U);
        contour(ptsX,ptsY,reshape(f,[nG,nG]),'color','r');
        
        f=ksdensity([CndRV{2}.SimRVX(iRtr,:)',CndRV{2}.SimRVY(iRtr,:)'],pts,'Bandwidth',U);
        contour(ptsX,ptsY,reshape(f,[nG,nG]),'color','b');
        hold on
        plot(CndRV{2}.SimRVX(iRtr,:)',CndRV{2}.SimRVY(iRtr,:)','b.')
        plot(Md,QY,'g-','linewidth',1.5)
        grid on
        axis tight
        title(sprintf('RtrPrd %g',Mrg(1).RtrPrd(iRtr)))
        xlim([0,10])
        ylim([0,10])
    end
end

%set(gca,'xscale','log','yscale','log')
%legend('QuantileMCX 0.37and Y 0.5','QuantileMC 0.5 and Y 0.5','MeanQuantile','QuantileMCYgXRV','location','NorthWest')