function Sml=SimulateMC(HT,Mrg,smlMethod,nRls,A,nA,RVX)
%Compute conditional quantiles of Y|X for given X values.

%Inputs
% HT Heffernan and Tawn model - fit upfront
% Mrg marginal model object, which includes the precalculated values of X
% smlMethod - method choice
% randX - X is generated at random from the marginal model
% retValX - X is generated from the return value model
% userX - user provided X, which is provided as a user input
% in that case RVX [nA X nBt x nRtr]
% A vector of bins to compute; (default computes omni over all bins)
% nA number of bins (default max(A));
% RVX - [nA x nBt x nRtr] return value calculations of X
%
% Outputs
% Sml = Populated structure of Y|X
% Y: [nA x nRls x nRtr] cond quantiles of Y given X
% StnMrg: [nA x nRls x nRtr] cond quantiles of Y given X on std
% Unif: [nA x nRls x nRtr] cond quantiles of Y given X on uniform margins
% X: [nA x nRls x nRtr] X
% X_Stn: [nA x nRls x nRtr] X std margins
% I: [nRls x 1] bootstrap index

% validate
if nargin<=2
    smlMethod="randX";
else
    smlMethod=validateString(smlMethod,["randX","retValX","userX"]);
end

% optional arguments for nRls, A and nA
if nargin<=3
    nRls=1000; %default number of MC draws
end

%Input A and Na super bins definitions
if nargin<=4
    % A=ones(HT.nBin,1); omni situation
    A=(1:HT.nBin)';% per original bin
end

%Input 5 number of super bins
if nargin<=5
    nA=max(A);
end

% populating conditional values of X
if nargin<=6 && strcmp(smlMethod,"UserX")
    error('User defined return values need to be provided in the case where smlMethod is UserX.')
end

% creating empty Sml objects for original, uniform and standard margins
Sml.A=NaN(nA,nRls);
Sml.Org=NaN(nA,nRls,HT.nDmn,HT.nRtr);
Sml.StnMrg=NaN(nA,nRls,HT.nDmn,HT.nRtr);
Sml.Unf=NaN(nA,nRls,HT.nDmn,HT.nRtr);

% get existing X should be resized to in terms of A instead nBin
Sml.I=randi(HT.nBoot,nRls,1);%decide which bootstrap sample to use over all bootstraps HT.nBoot
% simulate bin and sector allocation
if Mrg(1).Bn.nBin>1
    for iA=1:nA
        %Simulate covariate with right rate
        Rat=Mrg(1).Rat(A==iA,Sml.I);  %get rate of all observations
        RatCdf=cumsum(Rat)/sum(Rat); %get rate cdf
        Sml.A(iA,:)=sum(bsxfun(@gt,rand(1,nRls),RatCdf),1)'+1;%bin allocation
    end %iA
    Sml.X=reshape(SampleCovariateFromBin(Mrg(1).Bn,Sml.A(:)),size(Sml.A));
else
    Sml.A=ones(1,nRls);
    Sml.X=rand(1,nRls)*360; %todo think about non periodic case here.
end

switch smlMethod
    %sample uniform value
    case "randX"
        Sml.Unf(:,:,1)=rand(nA,nRls);
        Sml.Org(:,:,1)=Mrg(1).INV(Sml.Unf(:,:,1),Sml.I);
    case "retValX"
        U=NaN(nA,nRls,HT.nRtr); %U should be in the range [ P0, 1] where P0 is the non occurence rate.
        for iRtr=1:HT.nRtr
            for iA=1:nA
                rho=Mrg(1).Rat(A==iA,Sml.I); %annual rate of occurence
                LT=rho*Mrg(1).RtrPrd(iRtr); %poisson Rate
                UX=rand(1,nRls); %U should be in the range [ P0, 1] where P0 is the non occurence rate.
                U(iA,:,iRtr)=1+log(UX)./(LT);
            end %iA
        end %iRtr
        %adjust for return period  (inv of exp(-L*T*(1-C));
        U(U<HT.NEP(Sml.I)')=NaN;
        Sml.Unf(:,:,1,:)=permute(U,[1,2,4,3]);
        Sml.Org(:,:,1,:)=Mrg(1).INV(U,Sml.I);
    case "userX"
        %         for iBt=1:HT.nBoot
        %             Sml.Org(:,Sml.I==iBt,1,:)=RVX(:,iBt,:);
        %         end %iBt
        Sml.Org(:,:,1,:)=permute(RVX(:,Sml.I,:),[1,2,4,3]);
        Sml.Unf(:,:,1,:)=Mrg(1).CDF(Sml.Org(:,:,1,:));
end
%% simulate on standard scale
Sml.StnMrg(:,:,1,:)=INV_Standard(Mrg(1),Sml.Unf(:,:,1,:));
Sml.StnMrg(:,:,2:end,:)=HT.SampleYgX_Stn(Sml.StnMrg(:,:,1,:),Sml.I,Sml.A);
%% Transform to uniform
Sml.Unf=CDF_Standard(Mrg(1),Sml.StnMrg);
%% Transform back to original margin
for iDmn=2:HT.nDmn %loop over dimension
    Sml.Org(:,:,iDmn,:)=Mrg(iDmn).INV(Sml.Unf(:,:,iDmn,:),Sml.I,Sml.A);
end
% Populated Sml object
end %SimulateMC