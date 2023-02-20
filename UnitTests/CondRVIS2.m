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

%% FUll Monte Carlo
fprintf('Brute Force MC\n');
nRls=1e5;

MC=Simulate(HT,Mrg,nRls);

%% IS

IS=SimulateIS(HT,Mrg,nRls);

%% Conditional RV current code
figure(1);
clf;
nG=30;
Qnt=[0.025,0.05,0.25,exp(-1),0.5,0.75,0.95,0.975];
nQ=numel(Qnt);
QCol=jet(nQ);
for iBn=1:Mrg(1).Bn.nBin+1
    
    %% Monte Carlo Version
    if iBn<=Mrg(1).Bn.nBin
        tXMC=MC.Org(MC.A==iBn,1);
        tYMC=MC.Org(MC.A==iBn,2);
    else
        tXMC=MC.Org(:,1);
        tYMC=MC.Org(:,2);
    end
    Edg=linspace(min(tXMC),max(tXMC)*1.2,nG+1);
    Md=(Edg(1:end-1)+Edg(2:end))/2;
    A=discretize(tXMC,Edg);
    
    QMC=accumarray(A,tYMC,[nG,1],@(x){quantile(x,Qnt)},{NaN(1,nQ)});
    QMC=cell2mat(QMC);
    
    %% IS version old
    
    if iBn<=Mrg(1).Bn.nBin
        AIS=IS.A==iBn;
        tfogBn=exp(sum(IS.logfog(AIS,:),2)); %nRlsBin x 1  %exp(log(f(x))+log(f(y|x)))
    else
        AIS=true(numel(IS.A),1);
        %     logfog=reshape(IS.logfog,[nRls,Mrg(1).Bn.nBin,2]);
        %    W=reshape(IS.W,[nRls,Mrg(1).Bn.nBin]);
        tfogBn=exp(sum(IS.logfog(AIS,:),2)).*IS.W(AIS); %nRlsBin x 1
        
        %exp(log(f(x))+log(f(y|x)))*W
        
        % sum(w_i.*f_i(x))
        
    end
    tXIS=IS.Org(AIS,1); %nRlsBin x 1
    tYIS=IS.Org(AIS,2);   %nRlsBin x 1
    
    A=discretize(tXIS,Edg); %nRlsBin x 1
    
    QISOld=NaN(size(QMC));
    for iG=1:nG
        [tYSrt,srtI]=sort(tYIS(A==iG));
        tfogGrBn=tfogBn(A==iG);
        tfogGrBn=tfogGrBn(srtI);
        tP=cumsum(tfogGrBn)./sum(tfogGrBn);
        
        if any(tP)
            [sP,J]=unique(tP);
            if numel(sP)>1
                QISOld(iG,:)=interp1(sP,tYSrt(J),Qnt);
            end
        end
    end
    
    %% IS version new
    if 1
        if iBn<=Mrg(1).Bn.nBin
            AIS=IS.A==iBn;
            tfogBn=exp(sum(IS.logfog(AIS,:),2)); %nRlsBin x 1  %exp(log(f(x))+log(f(y|x)))
        else
            AIS=true(nRls,1);
            logfog=reshape(IS.logfog,[nRls,Mrg(1).Bn.nBin,2]);
            W=reshape(IS.W,[nRls,Mrg(1).Bn.nBin]);
            tfogBn=sum(exp(sum(logfog,3)).*W,2); %nRlsBin x 1
            
            %exp(log(f(x))+log(f(y|x)))*W
            
            % sum(w_i.*f_i(x))
            
        end
        tXIS=IS.Org(AIS,1); %nRlsBin x 1
        tYIS=IS.Org(AIS,2);   %nRlsBin x 1
        % tfogBn=exp(sum(IS.logfog(AIS,:),2)).*IS.W(AIS); %nRlsBin x 1
        
        A=discretize(tXIS,Edg); %nRlsBin x 1
        
        QIS=NaN(size(QMC));
        for iG=1:nG
            [tYSrt,srtI]=sort(tYIS(A==iG));
            tfogGrBn=tfogBn(A==iG);
            tfogGrBn=tfogGrBn(srtI);
            tP=cumsum(tfogGrBn)./sum(tfogGrBn);
            
            if any(tP)
                [sP,J]=unique(tP);
                if numel(sP)>1
                    QIS(iG,:)=interp1(sP,tYSrt(J),Qnt);
                end
            end
        end
    end
    
    %% Plot
    if  iBn<=Mrg(1).Bn.nBin
        figure(1);
        nPlt=ceil(sqrt(Mrg(1).Bn.nBin));
        subplot(nPlt,nPlt,iBn)
    else
        figure(2);
        clf;
    end
    
    plot(tXMC,tYMC,'k.','markersize',3)
    hold on
    
    for iQ=1:nQ
        plot(Md, QMC(:,iQ),'-','color',QCol(iQ,:),'linewidth',2)
        hold on
         plot(Md, QIS(:,iQ),'--','color',QCol(iQ,:),'linewidth',2)
        hold on
        plot(Md, QISOld(:,iQ),':','color',QCol(iQ,:),'linewidth',2)
        
    end
    axis tight
end