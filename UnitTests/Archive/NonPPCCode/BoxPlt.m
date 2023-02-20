function h=BoxPlt(Z,I,Scl,Pos,C,OutlineOn);
%Inputs:
% Z = data in long vector, data for different boxes concatenated  
% I = index on box data belongs to
% Scl = width of boxes relative to positions
% Pos = position of box plots on x axis
% C = colour, [nBx x 3] matrix
% OutlineOn = flag for drawing black outline around coloured boxes
  
%clf;
Cat=unique(I);
nCat=length(Cat);
if nargin<4
    Pos=1:nCat;
end
if nargin<5
    C='k';
end
if nargin < 6
   OutlineOn=0; 
end
if nCat > 1
    BoxWdt=Scl*(Pos(nCat)-Pos(1));
else
    BoxWdt=Scl;
end


count=0;
tMd=NaN(nCat,1);
for iBx=1:nCat;
    if (size(C,2)== 3) && (size(C,1)>1) %if entered different colour for each boxplot
        tC=C(iBx,:);
    else
        tC=C; 
    end
        
    tD=Z(I==Cat(iBx));
    count=count+1;
    
    tCmb=quantile(tD,[0.025,0.25,0.5,0.75,0.975],1);
    
    tXBx=Pos(iBx);
    plot([tXBx,tXBx],[tCmb(1),tCmb(5)],'-','color',tC,'linewidth',3)
    hold on
    tX=[tXBx-BoxWdt,tXBx+BoxWdt,tXBx+BoxWdt,tXBx-BoxWdt];
    tY=[tCmb(2),tCmb(2),tCmb(4),tCmb(4)];
    
    if OutlineOn  %black outline for boxes
        plot([tXBx-0.75,tXBx+0.75],[tCmb(1),tCmb(1)],'-','color','k','linewidth',1)
        plot([tXBx-0.75,tXBx+0.75],[tCmb(5),tCmb(5)],'-','color','k','linewidth',1)
        h=patch(tX,tY,tC,'EdgeColor','k') ;
    else
        h=patch(tX,tY,tC,'EdgeColor',tC) ;  
    end
    
    tMd(iBx)=tCmb(3);
end
plot(Pos,tMd,'.','color','w','MarkerFaceColor','w','markersize',10)
plot(Pos,tMd,'.','color','k','markersize',5)
if strcmp(version,'8.3.0.532 (R2014a)');
    set(gca,'xtick',round(sort(Pos)))
else;
    set(gca,'xtick',round(sort(Pos),1))
end;