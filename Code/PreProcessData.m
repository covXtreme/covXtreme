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

function [Rsp,Cvr,Asc] = PreProcessData(RspOrg,CvrOrg,AscOrg, VslHd)
   %Function to pre-process case study data: remove "0" responses due to 
   %...heading measurement error giving 90deg
   %% OUTPUTS
   % Rsp - main response
   % Cvr - covariate
   % Asc - associated response
   KpInd = ~any(VslHd==90 | AscOrg== 0,2);
   Asc = AscOrg(KpInd,:); %main response
   Rsp = RspOrg(KpInd); %covariate (direction)
   Cvr = CvrOrg(KpInd,:); %associated variable(s)

end
