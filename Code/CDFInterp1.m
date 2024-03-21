% Copyright © [2023] Shell Global Solutions International B.V. All Rights Reserved.
% SPDX-License-Identifier: Apache-2.0

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
