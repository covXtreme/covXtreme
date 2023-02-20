function pPltDrcSct(Drc,X,RspLbl,CvrLbl);

if nargin==0;
    clf;
    X=rand(100,2);
    Drc=rand(100,1)*360;
    CvrLbl={'Covariate'};
    RspLbl={'Response1';'Response2'};
end;

Bin=(22.5:45:337.5)';

nX=size(X,1);
A=nan(nX,1);

BinLbl=cell(8,1);
for iO=1:8;
    if iO<8;
        t=Drc>Bin(iO) & Drc<=Bin(iO+1);
        BinLbl{iO}=sprintf('(%g,%g]',Bin(iO),Bin(iO+1));
    else;
        t=Drc>Bin(iO) | Drc<=Bin(1);
        BinLbl{iO}=sprintf('%s in (%g,360] or (0,%g]',CvrLbl{1},Bin(iO),Bin(1));
    end;
    A(t)=iO;
end;

Lct=[3;6;9;8;7;4;1;2];
for iO=1:8;
    subplot(3,3,Lct(iO));
    plot(X(A==iO,1),X(A==iO,2),'color',pClr(iO),'marker','o','linestyle','none');
    if Lct(iO)==7;
        xlabel(RspLbl{1});
        ylabel(RspLbl{2});
    end;
    title(BinLbl{iO});
    box on; grid on;
end;

subplot(3,3,5); hold on;
for iO=1:8;
    plot(X(A==iO,1),X(A==iO,2),'color',pClr(iO),'marker','o','linestyle','none');
end;
box on; grid on;

return;
   