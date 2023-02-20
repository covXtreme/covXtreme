
clc; clear; 

n=1e5;

Mu=[0,0];
rho=0.5;
Sgm=[1,rho;rho,1];

LcPt=[3,0]; %lock point

%% Empircal density contour
Z=mvnrnd(Mu,Sgm,n);
clf;
subplot(2,1,1)
plot(Z(:,1),Z(:,2),'.k')
nG=10;
Bn=linspace(-4,4,nG);
[N,XEDGES,YEDGES] = histcounts2(Z(:,1),Z(:,2));

%%find bin for Lock point
A1=discretize(LcPt(1),XEDGES);
A2=discretize(LcPt(2),YEDGES);

XG=(XEDGES(1:end-1)+XEDGES(2:end))/2;
YG=(YEDGES(1:end-1)+YEDGES(2:end))/2;

[X,Y]=ndgrid(XG,YG);


hold on 
plot(LcPt(1),LcPt(2),'m.','markersize',20)
hold on 
contour(X,Y,N,N(A1,A2)*[1,1],'linewidth',2,'color','r')

xlim([-4,4])
ylim([-4,4])
%% IS density contour

nSmp=1e4;
U=(lhsdesign(nSmp,2)-0.5)*8;

f=mvnpdf(U,Mu,Sgm);

subplot(2,1,2)
scatter(U(:,1),U(:,2),40,f,'filled')
colorbar

nG=40;
BinEdg=linspace(-4,4,nG+1);
Grd=(BinEdg(1:end-1)+BinEdg(2:end))/2;

[X,Y]=ndgrid(Grd,Grd);

%%find bin for Lock point
A1=discretize(LcPt(1),BinEdg);
A2=discretize(LcPt(2),BinEdg);

Ind1=discretize(U(:,1),BinEdg);
Ind2=discretize(U(:,2),BinEdg);

A=sub2ind([max(Ind1),max(Ind2)],Ind1,Ind2);

fBin=accumarray(A,f,[],@mean,NaN);
fBin=reshape(fBin,[nG,nG]);


hold on 
contour(X,Y,fBin,fBin(A1,A2)*[1,1],'linewidth',2,'color','r')
hold on 
plot(LcPt(1),LcPt(2),'m.','markersize',20)
xlim([-4,4])
ylim([-4,4])