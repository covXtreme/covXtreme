function savePics(Nme,width,height,typ)

if nargin==1
width = 13;     % Width in inches
height = 13*2/3;    % Height in inches
end
if nargin<=3
   typ='pdf'; 
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


%print(Nme,'-djpeg','-r150');
switch typ    
    case 'pdf'
        print(Nme,'-dpdf','-r300'); %prints to pdf
    case 'jpg'
        print(Nme,'-djpeg','-r150'); %prints to jpg
    case 'png'
        print(Nme,'-dpng','-r150'); %prints to png
end
% print(Nme,'-depsc2','-r600');  

