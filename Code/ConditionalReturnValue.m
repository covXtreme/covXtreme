% Copyright © [2023] Shell Global Solutions International B.V. All Rights Reserved.
%
% SPDX-License-Identifier: Apache-2.0
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%      http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

function HT=ConditionalReturnValue(HT,Mrg,nRls,I)
%compute conditional return value Y | X using monte carlo simulation use return period as defined in marginal model.
%% INPUTS
% HT precomputed Heffernan and Tawn model 
% Mrg precomputed marginal model 
% nRls number of realisations in the Monte Carlo sample 
% I index for the bootstrap and number of samples 
%% OUTPUTS
% HT with populated simulatiob object 
if nargin<=2
    nRls=1000; %default number of MC draws
end

HT.nRtr=numel(Mrg(1).RtrPrd);
%preallocate
nAsc=HT.nDmn-1;
tRV.X_Stn=NaN(HT.nBin,nRls,HT.nRtr);
tRV.Y_Stn=NaN(HT.nBin,nAsc,nRls,HT.nRtr);
tRV.Y=NaN(HT.nBin,nAsc,nRls,HT.nRtr);

%% function - Get X part of the function
switch HT.RVMth
    case 1 %Random X from T-yr Maxium
        tRV.X=NaN(HT.nBin,nRls,HT.nRtr);
        I=randi(HT.nBoot,nRls,1);  %bootstrap samples to use
    case 2 %Fixed value from cuser input
        %nBin x nRtrPrd Fixed return value definition (only used in RVMth=2).
        tRV.X=permute(HT.RVFix,[1,3,2]);%gives you a [nBin x 1 x nRtrPrd]
        tRV.X=repmat(tRV.X,1,nRls,1);%need to check whether this is necessary
        %comment to see if it is necessary to include
        I=randi(HT.nBoot,nRls,1);  %bootstrap samples to use
        
    case 3 %Fixed bootstrap and user X value
        tRV.X=permute(HT.RVFix,[1,3,2]);%gives you a [nBin x 1 x nRtrPrd]
        tRV.X=repmat(tRV.X,1,nRls,1);%need to check whether this is necessary

        
end
%% function - Get residual components for the Monte Carlo analysis
tRV.I=I; %store bootstraps for use in threshold plot
%draw random resdiuals

if ~HT.SmpLclRsdOn   %if too many bins (few obs per bin), sample residuals globally
    Z=cell2mat(cellfun(@(x)x(randi(numel(x),HT.nBin,1))',HT.Rsd(I,:),'uniformoutput',false))'; 
    Z=permute(Z,[1,3,2]);
    
else
    %when your bins are big enough (decent no. of obs per bin), sample residuals locally from bin
    Z=NaN(nRls,nAsc,HT.nBin);
    for iBin=1:HT.nBin
        tRsd=cellfun(@(x,y)y(x==iBin,:),HT.RsdInd(I,:),HT.Rsd(I,:),'uniformoutput',false);
        J=cellfun(@length,tRsd)>0;
        Z(J,:,iBin)=cell2mat(cellfun(@(x)x(randi(size(x,1)),:),tRsd(J,:),'uniformoutput',false));
    end %iBin
    Z=permute(Z,[3,2,1]);
end
%% Return period obtaining for each of the methods
for iRtr=1:HT.nRtr %loop over return periods
    
    %% Sample from return value distribution within each bin
    rho=Mrg(1).Rat(:,I); %annual rate of occurence
    LT=rho*Mrg(1).RtrPrd(iRtr); %poisson Rate
    
    switch HT.RVMth
        case 1
            UX=rand(HT.nBin,nRls); %U should be in the range [ P0, 1] where P0 is the non occurence rate.
            P=1+log(UX)./(LT);    %adjust for return period  (inv of exp(-L*T*(1-C));
            P(bsxfun(@lt,P,HT.NEP(I)'))=NaN;  %this is the non-exceedence rate on standard scale           
            tRV.X(:,:,iRtr)=Mrg(1).INV(P,I); %RVX values on original scale
        case {2,3}
            P=Mrg(1).CDF(tRV.X(:,1,iRtr),(1:HT.nBin)',I); %transform uniform to standard margins using CDF
    end
    tRV.X_Stn(:,:,iRtr)=INV_Standard(Mrg(1),P); %transform form uniform to standard margins using CDF
    %compute Y using conditional model given X
    tX=permute(tRV.X_Stn(:,:,iRtr),[1,3,2]);
    tRV.Y_Stn(:,:,:,iRtr)=bsxfun(@times,HT.Alp(:,:,I),tX) + bsxfun(@power,tX,HT.Bet(:,:,I)).*bsxfun(@plus,HT.Mu(:,:,I),bsxfun(@times,HT.Sig(:,:,I),Z));
    % perform for the associated case
    for iAsc=1:nAsc
        %transform Y from standard to uniform margins using CDF
        UY=permute(Mrg(iAsc+1).CDF_Standard(tRV.Y_Stn(:,iAsc,:,iRtr)),[1,3,2]); %nBin x nAsc x nRls x nRtr
        %Transform Y to back original margins
        tRV.Y(:,iAsc,:,iRtr)=permute(Mrg(iAsc+1).INV(UY,I),[1,3,2]);
    end %iAsc
end %iRtr

%% Get Omni Value (covariate free) return value and its associated conditions
if HT.nBin > 1
    
    XOmni_Stn=max(tRV.X_Stn,[],1);
    [XOmni,J]=max(tRV.X,[],1);    %J is index of location of max
    % set up empty data frames for omnidirectional outcome
    YOmni=NaN(1,nAsc,nRls,HT.nRtr);
    YOmni_Stn=NaN(1,nAsc,nRls,HT.nRtr);
    for iRtr=1:HT.nRtr
        I=sub2ind([HT.nBin,nRls],J(:,:,iRtr),(1:nRls)); %need composite index of realisation and max location
        %original margins
        tYOmni=reshape(permute(tRV.Y(:,:,:,iRtr),[1,3,2]),HT.nBin*nRls,nAsc);
        YOmni(1,:,:,iRtr)=permute(tYOmni(I,:),[3,2,1]);
        %standard margins
        tYOmni_Stn=reshape(permute(tRV.Y_Stn(:,:,:,iRtr),[1,3,2]),HT.nBin*nRls,nAsc);
        YOmni_Stn(1,:,:,iRtr)=permute(tYOmni_Stn(I,:),[3,2,1]);
    end %iRtr
    
    tRV.Y=cat(1,tRV.Y,YOmni);
    tRV.Y_Stn=cat(1,tRV.Y_Stn,YOmni_Stn);
    tRV.X=cat(1,tRV.X,XOmni);
    tRV.X_Stn=cat(1,tRV.X_Stn,XOmni_Stn);   
end

%% Store Return Value simulation
HT.RV=tRV;
HT.RV.nRls=nRls;

end %ConditionalReturnValue
