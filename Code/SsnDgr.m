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

function SsnDgr=SsnDgr(DatJln)
%For a given vector of Julian dates, finds the seasonal degree
%All years have seasonal degrees in [0,360) (regardless of leap years etc)
%
%% INPUTS
%  DatJln [n x 1] Julian dates for storm peak data
%% OUTPUTS
%  SsnDgr [n x 1] Seasonal degree (in [0,360)) corresponding to Julian dates
% dealing with different versions of MATLAB
if verLessThan('matlab','8.4')
    tYr=datevec(DatJln); %datevector (has current year in first column)
    tYr=tYr(:,1); %extract year only
    tDay=DatJln-datenum(tYr,1,1)+1; %Day of year (in [0,360))
    tNmb=datenum(tYr,12,31)-datenum(tYr,1,1)+1; %Number of days in year (365 or 366)
    SsnDgr=360*tDay./tNmb; %Season in [0,360) (regardless of number of days in year)    
else
    tD =datetime(DatJln,'convertfrom','datenum');
    yr=year(tD);
    ly=@(yr)~rem(yr,400)|rem(yr,100)&~rem(yr,4);
    L=ly(yr);
    nDay=365*ones(length(tD),1);
    nDay(L)=366;
    SsnDgr=day(tD,'dayofyear')./nDay*360;
end
return