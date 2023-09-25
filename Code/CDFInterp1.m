% Copyright © [2023] Shell Global Solutions International B.V. All Rights Reserved.
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

function Pi=CDFInterp1(Z,Zi,P)
%"next" neighbour interpolation for a CDF
%% INPUTS
%Z should be sorted [n x 1]
%P (optional) should be an [n x 1] linspace vector of probabilities.
%Zi should be [p x q] set of place to interpolate onto
%% OUTPUTS
% Pi p x q probabilities on [0,1] corresponing to Zi

n=numel(Z);

if n==0
    Pi=ones(size(Zi));
elseif n==1
    Pi=double(Zi>=Z);
else  
    if nargin<=2 %cumulants evenly spaced                                          
        P=(1:n)'./n;
    end
    [~,I]=histc(Zi,[Z;Inf]);
    I(I==n+1)=n;
    Pi=zeros(size(Zi));
    Pi(I>0)=P(I(I>0));
    %end
end
