% Copyright © [2023] Shell Global Solutions International B.V. All Rights Reserved.
% SPDX-License-Identifier: Apache-2.0

function Dat=PeakPicking(Rsp,Cvr,Asc,IsPrd,NEP,RspLbl,CvrLbl,IsStrTrj)
%% INPUTS
% Rsp     [n x 1] vector of main response 
% Cvr     [n x nCvr] matrix of covariates (e.g. direction, season)
% Asc     [n x nAsc] matrix of other associated variables 
% IsPrd   [nCvr x 1] (boolean) vector flag for if covariate is periodic
% NEP     [1 x 1] quantile level used to find threshold
% RspLbl  [nAsc+1 x 1] cell arrray of labels for the response (first) and the associated  variables
% CvrLbl  [nCvr x 1]  cell arrray of labels for the covariates
%% OUTPUTS
% Dat Peak picked data set
% Dat.Y  [nPk x (1+nAsc)] vector of data with main response first
% Dat.X  [nPk x 1] vector of covariate (direction) data 

if nargin<8
    IsStrTrj=0;
end
validateattributes(Rsp, {'numeric'},{'vector'},'PeakPick','Y',1);
n=numel(Rsp);  %number of observations
validateattributes(Cvr, {'numeric'},{'nrows',n},'PeakPick','Thet',2);
nCvr=size(Cvr,2);
validateattributes(Asc, {'numeric'},{'nrows',n},'PeakPick','Asc',3);
validateattributes(IsPrd, {'numeric','logical'},{'numel',nCvr,'integer'},'PeakPick','IsPrd',4);
nAsc=size(Asc,2); %number of associated variables
validateattributes(NEP, {'numeric'},{'scalar','>=',0,'<=',1},'PeakPick','NEP',5);
validateattributes(RspLbl, {'cell'},{'numel',nAsc+1},'PeakPick','RspLbl',6);
validateattributes(CvrLbl, {'cell'},{'numel',nCvr},'PeakPick','CvrLbl',7);
validateattributes(IsStrTrj, {'numeric'},{'scalar','binary'},'PeakPick','IsStrTrj',8);

%remove any nans in orginal data.
I=isnan(Rsp) | any(isnan(Cvr),2) | any(isnan(Asc),2);
Rsp=Rsp(~I);
Cvr=Cvr(~I,:);
Asc=Asc(~I,:);

%% Exceedence threshold
Thr=quantile(Rsp,NEP); %Find threshold corresponding to non-exceedance probability NEP
Obs=(1:n)'; %observation numbers for original data
IExc=Rsp>Thr; %index for exceedences
ObsExc=Obs(IExc); %observation numbers for exceedence data

%% find storm periods
ObsDiff=diff(ObsExc); %difference between observation numbers of exceedances
GapInd=ObsDiff>1;  %location of gaps bigger than 1 observation
GapInd=[GapInd;true];  %add final observation as a gap
Up=circshift(GapInd,1);  %upcrossings
Dw=GapInd; %downcrossings
Prd=[ObsExc(Up),ObsExc(Dw)];  %storm periods: col1=start of storm; col2=end of storm

%% turn periods into index on exceedences
Del=0.1*1;
tPrd=[Prd(:,1)-Del,Prd(:,2)+Del]; %expand periods slightly to aviod boundary issues in histc.
[~,I]=histc(ObsExc,[-Inf;reshape(tPrd',[],1)]);  %Find indices of data belonging to each storm period (add inf to make sure periods are even)
I(mod(I,2)==1)=0;  %remove non-exceedance periods between storms (will have odd indices)
Ind=I/2; %Put storm indices on 1:1:nStorm scale instead of 1:2:nStorm

nExc=max(Ind); %number of exceedences
%% Get exceedences
RspExc=Rsp(IExc);
CvrExc=Cvr(IExc,:);
AscExc=Asc(IExc,:);

%% find storm peak maximum index
maxInd=accumarray(Ind,RspExc,[],@findmax,NaN); %index of maxima within storm
SSCnt=accumarray(Ind,Ind,[],@numel,NaN); %sea state count per storm
cmSSCnt=[0;cumsum(SSCnt(1:end-1))];%cumulative sea state count per storm
maxIndOrg=maxInd+cmSSCnt; %index of maxima within exceedences

%% Populate output vector
Dat.Y=NaN(nExc,1+nAsc); %Initialise empty response matrix
Dat.Y(:,1)=RspExc(maxIndOrg); %maxima within storm
Dat.Y(:,2:end)=AscExc(maxIndOrg,:); %value of associated response @ maxima
Dat.X=CvrExc(maxIndOrg,:);  %value of covariate (direction) @ maxima
Dat.RspLbl=RspLbl; %response label
Dat.CvrLbl=CvrLbl; %covariate label
Dat.IsPrd=IsPrd;  %periodic covariate flag
Dat.Prd=Prd; %storm periods: col1=start of storm; col2=end of storm

if IsStrTrj==1
    Dat.StrTrj.RA=cell(nExc,1+nAsc); %Initialise empty cell array for storm trajectories
    Dat.StrTrj.Cvr=cell(nExc, nCvr); %Initialise empty cell array for covariate trajectories.
    %% For each storm peak, get its trajectory.
    for iExc = 1:nExc
        % Find the start and end of each storm
        startIdx = Dat.Prd(iExc, 1); % start of the storm
        endIdx = Dat.Prd(iExc, 2);   % end of the storm

        % Extract the storm trajectory within the storm period
        stormRsp = Rsp(startIdx:endIdx);  % Response trajectory
        stormAsc = Asc(startIdx:endIdx, :);  % Associated variables trajectory
        stormCvr = Cvr(startIdx:endIdx, :); % Covariates

        % Store the storm trajectory
        Dat.StrTrj.RA{iExc, 1} = stormRsp; % Store response trajectory
        for iAsc = 1:nAsc
            Dat.StrTrj.RA{iExc, iAsc+1} = stormAsc(:, iAsc);  % Store associated variables trajectory
            Dat.StrTrj.Cvr{iExc, iAsc} = stormCvr(:, iAsc);
        end %iAsc
    end %iExc
end


nDmn=size(Dat.Y,2);
%% Plotting
%% creation of the figure directory
if ~exist('Figures','dir')
   mkdir('Figures') 
end
%% marginal plot
figure(1);
clf;
c=0;
for i=1:nDmn
    for iC=1:nCvr
        c=c+1;
        subplot(nDmn,nCvr,c)        
        if i==1
            plot(Cvr(:,iC),Rsp,'.','color',[1,1,1].*0.7);
        else
            plot(Cvr(:,iC),Asc(:,i-1),'.','color',[1,1,1].*0.7);
        end
        hold on
        grid on
        plot( Dat.X(:,iC),Dat.Y(:,i),'k.');        
        xlabel(CvrLbl{iC})        
        axis tight
        if IsPrd(iC)
            set(gca,'xtick',0:45:360,'xlim',[0,360])
        end
        ylabel(RspLbl{i})
    end %iC
end %i
savePics('Figures/Stg1_Data_Margins')

%% joint plot
if nDmn>1
    figure(2);
    clf;   
    nAsc=size(Asc,2);
    for iA=1:nAsc
        subplot(1,nAsc,iA);
        plot(Rsp,Asc(:,iA),'.','color',[1,1,1].*0.7);
        hold on
        grid on
        plot( Dat.Y(:,1),Dat.Y(:,iA+1),'k.');
        xlabel(RspLbl{1})
        ylabel(RspLbl{iA+1})
    end
end
savePics('Figures/Stg1_Data_Joint')

end

function I=findmax(Y)
%find location of maximum
[~,I]=max(Y);
end