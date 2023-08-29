% Copyright 2023 covXtreme
%
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

function F=rndsymlgs(alp,nD,n)
%generate random values from nD dimension symmetric logistics function with
%Frechet margins
%
%Algorithm 1.1 Alec Stephenson, Simulating Multivariate Extreme Value
%Distributions of Logistic Type, Extremes 6, 49-59, 2003
%% INPUTS
%alp 1 x 1 or n x 1 vector of the dependency parameter
%nD  is number of dimensions
%n is the number of realisations
%% OUTPUTS
%F nD x n random values on frechet margins

if unique(alp)>1
    n=numel(alp);
    nalp=n;
else
    nalp=1;
end

%step 1)  get T from unit exponentials
W=exprnd(1,[nD,n]);
T=W./sum(W);

%find mixture probabilities from recurrence relationship
p=NaN(nD,nD,nalp);
for iA=1:nalp
    p(1,1,iA)=1;
    for iD=2:nD
        p(iD,1,iA)=gamma(iD-alp(iA))./(gamma(iD).*gamma(1-alp(iA)));
        for jD=2:(iD-1)
            t1=(iD-1-alp(iA)*jD).*p(iD-1,jD)+alp(iA).*(jD-1).*p(iD-1,jD-1);
            p(iD,jD,iA)=t1./(iD-1);
        end
        p(iD,iD,iA)=alp(iA).^(iD-1);
    end
end
 
P=cumsum(shiftdim(p(end,:,:),1)); %[nD x nalp]  %cumulative probability

%step 2) Find k
U=rand(1,n);
k=sum(U>P,1)+1;
%step 3) get Z from gamma(k,1);
Z=gamrnd(k,1,[1,n]);
%step 4) find F on frechet margins
F=1./(Z.*(T.^(alp(1:nalp))));

