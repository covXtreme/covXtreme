%% Stage6_StormTrajectorySimulation

%% Set-up
RtrPrd=100;
nSml=10;

%% Using marginal model for dominant variate only
iDmn=1; %Always iDmn=1 for dominant variate

load(sprintf('Output/MM%g',iDmn),'MM');

% Simulate event set
Sml=cell(nSml,1);
for iS=1:nSml;
    nStr = poissrnd(MM.nDat*RtrPrd/MM.Yrs);
    Sml{iS}=sample_MC(MM,nStr); %This is input to storm matching
end;

% Diagnostic plot
load('Output/Data','Dat');
clf; hold on;
for iS=1:nSml; 
    tSml=sample_MC(MM,MM.nDat);
    plot(sort(Dat.Y(:,1)),sort(tSml.Org),'k.');
end;
xlabel('Sorted originial sample');
ylabel('Sorted simulated sample');

%% Using HT model for all variables

% Load HT model
load('Output/HT','HT');
smlMethod="randX";

% Load and combine marginal models
FM=cellstr(ls('Output/MM*.mat'));
for iM=1:numel(FM) %load all marginal models
    load(sprintf('Output/MM%g',iM),'MM');
    if iM==1
        Mrg=MM;
    else
        Mrg=cat(1,Mrg,MM);
    end
    clear MM;
end

% Simulate event set
SmlHT=cell(nSml,1);
fprintf(1,'HT simulation: ');
for iS=1:nSml;
    fprintf(1,'%g ',iS); if rem(iS,10)==0; fprintf(1,'\n'); end;
    nStr = poissrnd(Mrg(1).nDat*RtrPrd/Mrg(1).Yrs);
    SmlHT{iS}=SimulateMC(HT,Mrg,smlMethod,nStr); %This is input to storm matching
end;

% Diagnostic plot
load('Output/Data','Dat');
tSml=cell(nSml,1);
fprintf(1,'HT simulation: ');
for iS=1:nSml; 
    fprintf(1,'%g ',iS); if rem(iS,10)==0; fprintf(1,'\n'); end;
    tSml{iS}=SimulateMC(HT,Mrg,smlMethod,Mrg(1).nDat); 
end;
clf;
for iDmn=1:HT.nDmn;
    subplot(1,HT.nDmn,iDmn); hold on;
    for iS=1:nSml;
        plot(sort(Dat.Y(:,iDmn)),sort(tSml{iS}.Org(1,:,iDmn)),'k.');
    end;
    xlabel('Sorted originial sample');
    ylabel('Sorted simulated sample');
    title(Dat.RspLbl(iDmn));
end; %iDmn

%% Naive matching

load('Output/Data','Dat');

load('Output/Bin','Bn');

% Loop over storm peaks in Sml

% For each storm peak, identify the nearest nNgh storms in Dat.StrTrj
% This will use 
% - only dominant variate for "Kevin's" method, 
% - all variates for HT method
% Only match with a storm which has its storm peak in the same covariate bin

% Adjust matched storm trajectory based on physical constraints (do later)

% Save the full matched storm trajectory
% - bin index
% - nAsc+1 values of dominant and associated variates


%% Create storm trajectory bin allocations
fprintf(1,'Calculating bin allocations for historical trajectories\n');
for iS=1:Bn.n;

    tCvr=zeros(Bn.n,1);
    tStrTrjCvr=Dat.StrTrj.Cvr{iS,:};
    tCvr(1:size(tStrTrjCvr,1))=tStrTrjCvr;
    iVrb = 0;
    tBn=BinAllocation(Bn,tCvr, iVrb);
    Dat.StrTrj.A{iS,:}=tBn.A(1:size(tStrTrjCvr,1));

end;
save('Output/Data','Dat');

%% Storm matching to allocate historical trajectories to simulated storms
% Match using the dominant variate only, and in the covariate bin
% corresponding to the dominant variate only.
nNgh=10; %Number of near neighbour historical storm trajectories to consider
fprintf(1,'Storm matching to allocate historical trajectories to simulated storms:\n');
iE=0;
for iS=1:nSml;

    fprintf(1,'%g ',iS); if rem(iS,10)==0; fprintf(1,'\n'); end;

    for iR=1:Sml{iS}.nRls;

        iE=iE+1;

        tBn=Sml{iS}.A(iR); %bin for realisation iR from simulation iS

        tMB=find(Bn.A==tBn); %all historical matched with correct bin allocation

        [~,tD2]=sort((Dat.Y(tMB,1)-Sml{iS}.Org(iR)).^2);

        tMV=tMB(tD2(1:nNgh));

        tMtc=tMV(randi(nNgh,1));

        Sml{iS}.StrTrj.RA{iR,:}=Dat.StrTrj.RA{tMtc,:};
        Sml{iS}.StrTrj.Cvr{iR,:}=Dat.StrTrj.Cvr{tMtc,:};
        Sml{iS}.StrTrj.A{iR,:}=Dat.StrTrj.A{tMtc,:};
        Dgn(iE,:)=[max(Dat.StrTrj.RA{tMtc,1}) Sml{iS}.Org(iR) Sml{iS}.A(iR)];

    end; %iR
end; %iS

%% Diagnostic plot for storm matching
fprintf(1,'Diagnostic plot for storm matching\n');
clf; 
for iB=1:Bn.nBin;
    subplot(ceil(sqrt(Bn.nBin)),ceil(sqrt(Bn.nBin)),iB); hold on;
    plot(Dgn(Dgn(:,3)==iB,1),Dgn(Dgn(:,3)==iB,2),'k.');
    tMxm=max(Dgn(:,1));
    plot([1 tMxm],[1 tMxm],'r');
end;

% Next steps
%
% 1. Tidy up code
% 1.1 Tidy up diagnostic plot (borrow from covXtreme)
% 2. Rescale matches trajectories to correct "storm peak" values at the
% dominant variable.
% 3. Estimate CDF for storm peak variables and for sea state variables(Omni
% & covariate bin).
% 4. can we write a faster version avoiding structures??





