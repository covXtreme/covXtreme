function h=NewBoxPlt(Qnt,EffNmb);
%function h=NewBoxPlt(Qnt,EffNmb);
%
%S.Sam, P.Jonathan 20170605
%Modified from E.Ross code
%
%Uses 5 quantiles of distribution to set box-whisker limits
%(Previous code read in sample and worked out quantiles)
%
%Inputs:
% Qnt 5 x nBox quantiles to use (5 quantiles for each of nBox box plots) 

hold on;

nBox=size(Qnt,2);
Pos=(1:nBox)';
BoxWdt=0.25;

tMd=NaN(nBox,1);
for iB=1:nBox;

    tCmb=Qnt(:,iB);
    
    tMd(iB)=tCmb(3);
    
    tXBx=Pos(iB);    
    plot([tXBx,tXBx],[tCmb(1),tCmb(5)],'-','color','k','linewidth',3)
    
    tX=[tXBx-BoxWdt,tXBx+BoxWdt,tXBx+BoxWdt,tXBx-BoxWdt];
    tY=[tCmb(2),tCmb(2),tCmb(4),tCmb(4)];
    h=patch(tX,tY,'k','EdgeColor','k') ;
    
end
plot(Pos,tMd,'.','color','k','MarkerFaceColor','k','markersize',12)
plot(Pos,tMd,'.','color','w','MarkerFaceColor','w','markersize',10)
plot(Pos,tMd,'.','color','k','markersize',5)

if size(Pos,1)==1;
    set(gca,'xlim',[min(Pos)-EffNmb*0.5 max(Pos)+EffNmb*0.5]);
else;
    set(gca,'xlim',[min(Pos)-0.5 max(Pos)+0.5]);
end;

if strcmp(version,'8.3.0.532 (R2014a)');
    set(gca,'xtick',round(sort(Pos)))
else;
    set(gca,'xtick',round(sort(Pos),1))
end;