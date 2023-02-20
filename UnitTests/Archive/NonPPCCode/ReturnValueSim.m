function RV = ReturnValueSim(HT,Mrg,nSim,SvFld)
%Function to calculate return values (given marginal and HT fits in Stg1 to
%Stg4
%Inputs: 
% - HT = [1x1] instance of HeffernanTawn class
% - Mrg = [2x1] instance of MarginalModel class
% - RtrPRd = [nRtrPrd x 1] vector of return periods in years
% - nSim = scalar number of simulation repetitions
% Output:
%  - RV = data structure containing return values
%    o RV.X = [nRtrPrd x nSim x 3] matrix of return values where cols =
%    Hs,Srg|Hs,Srg
%% 
RV.nSim = nSim; %no. repetitions of simulation
RV.RtrPrd = Mrg(1).RtrPrd; %return period(s) in years
RV.nRtrPrd = numel(RV.RtrPrd); %number of return periods considered
nObsPerYr=numel(Mrg(1).Y)/Mrg(1).Yrs; %ave no. obs per year
RV.nStr = int64(RV.RtrPrd.*nObsPerYr); %no. of storms required for return period (HT.Simulate requires integer input)
RV.X = NaN(RV.nRtrPrd,RV.nSim,4); %col1 = Hs, col2 = Srg|Hs, col3 = Srg, col4=Drc
RV.Lbl = {Mrg(1).Lbl,sprintf('%s | %s',Mrg(2).Lbl,Mrg(1).Lbl),sprintf(Mrg(2).Lbl,'%s')}; 

nSct=size(Mrg(1).Thr,1);

%% Set up parfor for simulation
nRtrPrd = numel(RV.RtrPrd);
nPar = RV.nRtrPrd*RV.nSim;
tNmbRls = RV.nStr;
%initialise Hs, Srg|Hs and Srg matrices
tHs = cell(nPar,1);
tSrgHs = cell(nPar,1);
tSrgHsSct = cell(nPar,nSct);
tSrg = cell(nPar,1);
tDrc = cell(nPar,1);
ParfInd = NaN(nPar,2);
tSimNmbr = repmat(1:RV.nSim,nRtrPrd,1)';
tRtrPrdNmbr = repmat(1:nRtrPrd,RV.nSim,1)';
ParfInd(:,1) = tSimNmbr(:); %simulation index
ParfInd(:,2) = tRtrPrdNmbr(:); %return period index
%parfor iPar = 1:nPar %loop over sim repetitions
for iPar = 1:nPar    
    nStr = tNmbRls(ParfInd(iPar,2));
    Sml=Simulate(HT,Mrg,nStr); %Simulate RtnPrd-worth of peaks
    [tHs{iPar},IMaxHs] = max(Sml.Org(:,1)); %largest Hs
    tSrgHs{iPar} = Sml.Org(IMaxHs,2);  % Srg associated with (given) largest Hs
    tSrg{iPar} = max(Sml.Org(:,2)); %largest Srg
    tDrc{iPar} = Sml.A(IMaxHs);
    for iS=1:nSct;
        tOrg=Sml.Org(Sml.A==iS,:);
        [~,IMaxHs] = max(tOrg(:,1)); %largest Hs
        tSrgHsSct{iPar,iS} = tOrg(IMaxHs,2);
    end;
end
%reshape results of parfor
RV.X(:,:,1) = reshape(cell2mat(tHs),RV.nRtrPrd,RV.nSim);
RV.X(:,:,2) = reshape(cell2mat(tSrgHs),RV.nRtrPrd,RV.nSim);
RV.X(:,:,3) = reshape(cell2mat(tSrg),RV.nRtrPrd,RV.nSim);
RV.X(:,:,4) = reshape(cell2mat(tDrc),RV.nRtrPrd,RV.nSim);
for iS=1:nSct;
    RV.Sct(:,:,iS)=reshape(cell2mat(tSrgHsSct(:,iS)),RV.nRtrPrd,RV.nSim);
end;

save(fullfile(SvFld,'RV'),'RV')

end