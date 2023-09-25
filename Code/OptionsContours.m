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

classdef OptionsContours
    %Class of options for estimation of contours using PPC
    %
    %Four contour estimation methods are included
    % Exc = Exceedance contour
    % HTDns = Heffernan and Tawn Density contour
    % Hus = Huseby contour
    % HusCln = Huseby Cleaned contour
    
    properties %(with default values)
        nGrd=50;        %Number of grid points for each of x and y when rectangular gridding needed
        nPnt=200;       %Number of points on the contour
        SmtWdtC=5;      %Smoothing width for the huseby contour C function (see Huseby et al 2015)
        BndWdtScl=0.02; %Band width scale for HTDns contour     
        Mth={};         %Cell array of set of contour Methods to use
        nSml=1e6;       %Number to simulate: number of importance samples to use
    end %properties
    
    methods
        % Check valid input for nSml
        function obj=set.nSml(obj,nSml)            
            validateattributes(nSml,{'numeric'},{'integer','positive','scalar'});
            obj.nSml=nSml;            
        end %set.nSml
        
        % Check valid input for nGrd
        function obj=set.nGrd(obj,nGrd)
            validateattributes(nGrd,{'numeric'},{'integer','positive','scalar'});
            obj.nGrd=nGrd;            
        end %set.nGrd
        
        % Check valid input for nPnt
        function obj=set.nPnt(obj,nPnt)            
            validateattributes(nPnt,{'numeric'},{'integer','positive','scalar'});
            obj.nPnt=nPnt;            
        end %set.nPnt

        % Check valid input for SmtWdtC
        function obj=set.SmtWdtC(obj,SmtWdtC)            
            validateattributes(SmtWdtC,{'numeric'},{'integer','positive','scalar'});
            obj.SmtWdtC=SmtWdtC;            
        end %set.SmtWdtC
        
        % Check valid input for BndWdtScl
        function obj=set.BndWdtScl(obj,BndWdtScl)            
            validateattributes(BndWdtScl,{'numeric'},{'positive','scalar'});
            obj.BndWdtScl=BndWdtScl;            
        end %set.BndWdtScl
        
        % Check valid input for Mth
        function obj=set.Mth(obj,Mth)            
            Mth=string(Mth); %handles Mth specified as text string not cell
            nMth=length(Mth);
            obj.Mth=cell(nMth,1);
            for i=1:nMth
                tmp=validatestring(Mth{i},{'Exc','HTDns','Hus','HusOld'});
                obj.Mth{i}=tmp; %tmp contains corrected names (e.g. case errors)
            end
            obj.Mth=unique(obj.Mth); %eliminate duplicates if present            
        end %set.Mth
        
    end %methods
    
end %OptionsContours