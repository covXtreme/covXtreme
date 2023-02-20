function XNew=MovingMeanPeriodic(XOld,SmtWdt);

if nargin==0;
    
    XOld=sin(0.1*(1:100)');
    
end;

nX=size(XOld,1);
tX=[XOld;XOld;XOld];

tX=movmean(tX,SmtWdt);

XNew=tX(nX+1:2*nX);

return;