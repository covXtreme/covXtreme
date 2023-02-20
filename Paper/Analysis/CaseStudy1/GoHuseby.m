clear all;
clc;
clf;

load PhlCntXY Cnt;

for iR=1:3;
    for iD=1:9;
        subplot(3,9,(iR-1)*9+iD); hold on;
        X=Cnt.XY{1}(:,:,iD,iR);
        plot(X(:,1),X(:,2),'k.-');
        X=CleanHuseby(X);
        plot(X(:,1),X(:,2),'r.-');
        drawnow;
    end;
end;
pDatStm('x: different directions, y: different return periods');
pGI('CleanHuseby',2);

