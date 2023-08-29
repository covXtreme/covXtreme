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

function savePics(Nme,width,height,typ)
% Saving pictures function
%% INPUTS 
% Nme - string filename
% width - saved width of figure in inches
% height - saved height of figure in inches
% typ - filetype 
if nargin==1
    width = 11;     % Width in inches
    height = 9*2/3;    % Height in inches
end

if nargin<=3
   typ='png'; 
end

papersize = [width,height];

if ~strcmp(get(gcf,'WindowStyle'),'docked')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
end

set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
set(gcf, 'PaperSize',papersize);
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
fontname='Helvetica';
set(findall(gcf,'type','text'),'fontname',fontname)
set(findall(gcf,'Type','axes'),'fontname',fontname)

switch typ    
    case 'pdf'
        print(Nme,'-dpdf','-r300'); %prints to pdf
    case 'jpg'
        print(Nme,'-djpeg','-r150'); %prints to jpg
    case 'png'
        print(Nme,'-dpng','-r150'); %prints to png
end
