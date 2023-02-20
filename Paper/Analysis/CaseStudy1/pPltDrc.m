function pPltDrc(Drc,X);

if nargin==0;
    X=rand(100,2);
    Drc=rand(100,1)*360;
end;

Bin=(22.5:45:337.5)';

nX=size(X,1);
A=nan(nX,1);

for iO=1:8;
    if iO<8;
        t=Drc>Bin(iO) & Drc<=Bin(iO+1);
    else;
        t=Drc>Bin(iO) | Drc<=Bin(1);
    end;
    A(t)=iO;
end;

Lct=[2;3;6;9;8;7;4;1];
for iO=1:8;
    subplot(3,3,Lct(iO));
    plot(X(A==iO,1),X(A==iO,2),'ko');
    pAxsLmt;
end;

return;
   