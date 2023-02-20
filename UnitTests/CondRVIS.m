clc; clear; close all;
load('Output\HT')
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

iRtr=1;

%% FUll Brute Force Monte Carlo

nRls=2e5;
Mrg(1).RtrPrd=100;

BruteOn=false;

if BruteOn
    fprintf('Brute Force MC\n');
    NRep=400;
    MCSml.X=NaN(NRep,1);
    MCSml.Y=NaN(NRep,1);    
    
    for iRep=1:NRep
        fprintf('#')
        if rem(iRep,20)==0
            fprintf('\n')
        end
        Cnt=poissrnd(sum(Mrg(1).Rat(:,1)).*Mrg(1).RtrPrd(iRtr));
        MC=Simulate(HT,Mrg,Cnt);
        [MCSml.X(iRep),I]=max(MC.Org(:,1));
        MCSml.Y(iRep)=MC.Org(I,2);
    end
end
%% Conditional RV current code


HT=HT.ConditionalReturnValue(Mrg,nRls);

%%

Sml=SimulateIS(HT,Mrg,nRls);
C=[hsv(HT.nBin);[0,0,0]];
%LT=Mrg(1).RtrPrd.*Mrg(1).Rat
PGP=Mrg(1).CDF(Sml.Org(:,1)',(1:HT.nBin),Sml.I); %nRls x nBin
logfGP=Mrg(1).LogPDF(Sml.Org(:,1)',(1:HT.nBin),Sml.I); %nRls x nBin

LT=Mrg(1).Rat(:,Sml.I).*shiftdim(Mrg(1).RtrPrd(iRtr),-1);
RP=exp(-LT.*(1-PGP)); %return value CDF nBin x nRls x nRtr
fxRV=exp(-LT.*(1-PGP)+log(LT)+logfGP);  %return value density in single bin  nBin x nRls x nRtr

RPOmni=prod(RP,1); %1 x nRtr Omni CDF
fxRVOmni=sum(fxRV./RP,1).*RPOmni;% nRls x nRtr 
fxRVOmni(prod(RP,1)==0)=0;
fxRVOmni=fxRVOmni';
fxRV2=NaN(nRls,1);
for iBn=1:3
    A=Sml.A==iBn;
    %fxRV2(A)=fxRV(iBn,A,iRtr);   
    fxRV2(A)=fxRV(iBn,A,iRtr);   
end
wIS=(fxRV./RP).*(RPOmni./fxRVOmni');% nRls x nRtr
wIS=wIS./sum(wIS,1);
wIS(isnan(wIS))=0;
wIS(isinf(wIS))=0;

wIS2=NaN(nRls,1);
for iBn=1:3
    A=Sml.A==iBn;
    %fxRV2(A)=fxRV(iBn,A,iRtr);   
    wIS2(A)=wIS(iBn,A,iRtr);   
end
%fxRVOmni(prod(RP,1)==0)=0;
%fxRVOmni=shiftdim(fxRVOmni,1);
% fxRVOmni=fxRV2./(sum(fxRV(:,:,iRtr),1))';
% fxRVOmni(isnan(fxRVOmni))=0;

%Sml.logfxRV=-LT.*(1-PGP)+log(LT)+logfGP;  %log return value density  nBin x nRls x nRtr
nX=300;

% 
% Sml.logfxRV2=NaN(nRls,1);
% Sml.logfGP2=NaN(nRls,1);
% for iBn=1:3
%     A=Sml.A==iBn;
% %    Sml.logfxRV2(A)=Sml.logfxRV(iBn,A,iRtr);
%     Sml.logfGP2(A)=logfGP(iBn,A);
% end
% 

%%

for j=1:2 %looping X and Y 
    figure(j);
    clf;
    X=linspace(min(Sml.Org(:,j)),max(Sml.Org(:,j))*1.2,nX)';
    for iBn=1:4 %3 bins + omni
    
        if j==1  %X
            P=WeightedCDF(X,HT.RV.X(iBn,:,iRtr)');
        else %Y|X
            P=WeightedCDF(X,squeeze(HT.RV.Y(iBn,1,:,iRtr)));
        end
        %% by Bin
        if iBn<4 %bin 
            A=Sml.A==iBn;
            if j==1 %X
                 
                 %r(x)  = return value density
                CndRV=WeightedCDF(X,Sml.Org(A,j),fxRV2(A));
            else %Y|X
                %f(y|x)= exp(sum(Sml.logfog(:,2)));
                %r(x)  = return value density in bin A
                CndRV=WeightedCDF(X,Sml.Org(A,j),exp(Sml.logfog(A,j)).*fxRV2(A));
            end
            
        else
            %% omni
            %f(y|x)= exp(sum(Sml.logfog(:,2)));
            %r(x)
            if j==1% X
                CndRV=WeightedCDF(X,Sml.Org(:,j),Sml.W.*fxRVOmni);
            else %Y| X
                fog=exp(Sml.logf(:,j)-sum(Sml.logg,2));
                fog(isnan(fog))=0;
                
                % f(y | x, theta) f(theta | x) f(x)
                CndRV=WeightedCDF(X,Sml.Org(:,j),fog.*wIS2.*fxRVOmni);
            end
        end
        
        subplot(2,1,1)
        plot(X,P,'color',C(iBn,:));
        hold on
        plot(X,CndRV,'--','color',C(iBn,:),'handlevisibility','off');               
        if iBn==4
        legend([Mrg(1).Bn.BinLbl,'Omni'],'location','northwest')
        end
        if BruteOn
            plot(sort(MCSml.Y),linspace(0,1,NRep),'-m','linewidth',2,'handlevisibility','off')
        end
        if j==1
            title('Marginal Return Value CDF X')
        else
            title('Conditional Return Value CDF Y|X')
        end
        
        subplot(2,1,2)
        plot(X(2:end),diff(P),'color',C(iBn,:));
        hold on
        hold on
        plot(X(2:end),diff(CndRV),'--','color',C(iBn,:));
        
        grid on
        
        %xlim([0,1])
        if j==1
            title('Marginal Return Value PDF X')
        else
            title('Conditional Return Value PDF Y|X')
        end

    end
end
figure(1);
subplot(2,1,1)
xlim([30,50])
subplot(2,1,2)
xlim([30,50])
%%
% figure(6);
% clf;
nG=60;
Edg=linspace(min(HT.RV.X(iBn,:,iRtr)),max(HT.RV.X(iBn,:,iRtr)),nG+1);
Md=(Edg(2:end)+Edg(1:end-1))/2;
% Cnt=NaN(nG,3);
% for iBn=1:4
% %     subplot(2,1,1)
% %     plot(sort(HT.RV.X(iBn,:,iRtr)),linspace(0,1,nRls),'color',C(iBn,:),'linewidth',2);
% %   hold on
% %     subplot(2,1,2)
%     Cnt(:,iBn)=histcounts(HT.RV.X(iBn,:,iRtr),Edg);
% %     histogram(HT.RV.X(iBn,:,iRtr),Edg,'EdgeColor',C(iBn,:),'displaystyle','stairs','linewidth',2);
% %     hold on  
% end
% %title('Marginal Return Value X');
% EmpRat=Cnt./sum(Cnt,2);


A=discretize(HT.RV.X(4,:,iRtr)',Edg);
Cnt2=NaN(nG,3);
for i=1:3   
    IMx=HT.RV.X(1:3,:)==HT.RV.X(4,:);
    Cnt2(:,i)=accumarray(A,IMx(i,:)',[nG,1],@sum,0);
end
EmpRat2=Cnt2./sum(Cnt2,2);
%%
figure(3);
clf
for j=1
   % subplot(1,2,j)
    for i=1:3
        J=Sml.A==i;
        plot(Sml.Org(J,j),wIS2(J),'.','color',C(i,:),'markersize',1,'handlevisibility','off');
        hold on
       % plot(Md,EmpRat(:,i),'-','color',C(i,:)*0.8,'linewidth',4);
        plot(Md,EmpRat2(:,i),'--','color',C(i,:)*0.8,'linewidth',1);
    end
    sgtitle('Return Probability it came from bin')
    
    ylabel('Weight')
    if j==1
        xlabel('X')
    else
        xlabel('Y')
    end
    legend(Mrg(1).Bn.BinLbl)
    xlim([0,200])
end
axis tight
%%
figure(4);
clf;
I=Sml.Org(:,1)>20;
subplot(2,1,1)
scatter(Sml.Org(I,1),Sml.Org(I,2),3,exp(sum(Sml.logfog(I,:),2)),'filled')
title('Joint Density')

subplot(2,1,2)
scatter(Sml.Org(I,1),Sml.Org(I,2),3,exp(Sml.logfog(I,2)),'filled')
title('Conditional Density Y | X')
caxis([0,100])