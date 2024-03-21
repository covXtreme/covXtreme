% Copyright © [2023] Shell Global Solutions International B.V. All Rights Reserved.
% SPDX-License-Identifier: Apache-2.0

function F=rndasymlgs(alp,theta,nD,n)
%generate random values from nD dimension asymmetric logistic function with
%Frechet margins
%
%Algorithm 2.1 Alec Stephenson, "Simulating Multivariate Extreme Value
%Distributions of Logistic Type", Extremes 6, 49-59, 2003
%
%% INPUTS
%alp   : dependency parameter
%theta : weighting parameter
%nD    : is number of dimesions
%n     : the number of realisations
%% OUTPUTS
%F random valeus

%find all sets
nB=factorial(nD)+1;
B=NaN(nB,nD);
%find all possible sets for numbers 1:nD
C=0;
for iD=1:nD
    t1=nchoosek((1:nD)',iD);
    for jD=1:size(t1,1)
        B(jD+C,1:iD)=t1(jD,:);
    end
    C=C+size(t1,1);
end %iD
%-------------------------------------------------------------------------
%generate symmetric values
Z=cell(nB,1);
for b=1:nB
    nD=sum(~isnan(B(b,:))); %number of dimensions
    if nD==1 %1d is just unit frechet
        Z{b}=1./(-(log(rand(1,n))));
    else
        Z{b}=rndsymlgs(alp(b),nD,n);
    end        
end %b
%-------------------------------------------------------------------------
%find max over sets for assymetric values
F=NaN(nD,n);
for iD=1:nD
    I=find(any(B==iD,2));    
    tI=find(B(I(1),:)==iD);
    F(iD,:)=theta{I(1)}(tI).*Z{I(1)}(tI,:);
    for iB=2:length(I)
        tI=find(B(I(iB),:)==iD);
        F(iD,:)=max(F(iD,:),theta{I(iB)}(tI).*Z{I(iB)}(tI,:));
    end %iB
end %iD
