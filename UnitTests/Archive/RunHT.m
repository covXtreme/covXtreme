clear ;clc; close all;
%% simulate some data from assymertic logistic model
nD=2;  %number of dimensions
n=10000;  %number of data
alp=[NaN,NaN,0.5];  %three parameters in nD=2;
theta={0,0.5,[1,0.5]};
F=rndasymlgs(alp,theta,nD,n)';  %simulated on unit frechet scale

%% transform to gumbel scale (this would be slightly different if we were fitting the marginal model).
G=log(F);

%% Fit heffernan and Tawn model
Thr=4;  %threshold
HT=HeffernanTawn(G(:,1),G(:,2),Thr);

%% simulate under HT model
nRls=10000;
[SmlY,SmlX]=simulate(HT,nRls);

%% X range over which to draw contour
CndX=7;  %X value to draw contour from  (could be Return Value)
nPon=100;  %how many points for which to draw contour
XRng=linspace(Thr,CndX,nPon); %range in X for contour

%% Constant Exceedence
% contours of constant
% 1) define extreme X value and find conitional Y (could be median or mode)
% 2) 
CndI=SmlX>CndX;  %exceedences of CndX
CndY=median(SmlY(CndI));  %median(Y | X > CndX)
Iup=CndI & SmlY>CndY;    %exceedence points above median Y  (for top part of contour)
Idw=CndI & SmlY<=CndY;   %exceedence points above median Y  (for bottom part of contour)
nUp=sum(Iup);   %number of exceedences above median Y
nDw=sum(Idw);   %number of exceedences above median Y

I=bsxfun(@gt,SmlX,XRng);  
nAll=sum(I);  %total observations for each XRng value

%find contour above and below which preserves the number of exceedences
YUp=NaN(nPon,1);
YDw=NaN(nPon,1);
for iP=1:nPon  %loop over contour points
    YDw(iP)= quantile(SmlY(I(:,iP)),nDw/nAll(iP));   
    YUp(iP)= quantile(SmlY(I(:,iP)),(nAll(iP)-nUp)/nAll(iP));
end

%% Constant density
%find kernel density of HT residuals to define shape density
[f,x]=ksdensity(HT.Rsd);  %only works for unimodal
[maxdns,J]=max(f);  %find density at the mode
% Density in X is Gumbel 
Gmlfmax=exp(-(CndX +exp(-CndX)));   %X density at the max 
GmlfRng=exp(-(XRng +exp(-XRng)));   %X density over the range
c=Gmlfmax*maxdns;  %c is  joint density level to preserve.
ZDns=Gmlfmax*maxdns./GmlfRng;  %density of Rsd which preserves joint density c
ZDw=interp1(f(1:J),x(1:J),ZDns);  %interpolate kernel onto needed densities
ZUp=interp1(f(J:end),x(J:end),ZDns); 
%Transform HT.Rsd distribution onto response
ZDwHT=HT.Prm(1).*XRng +XRng.^HT.Prm(2).*(HT.Prm(3) + HT.Prm(4).*ZDw);
ZUpHT=HT.Prm(1).*XRng +XRng.^HT.Prm(2).*(HT.Prm(3) + HT.Prm(4).*ZUp);

%% plots
clf;
hold on
plot(SmlX,SmlY,'.','color',[1,1,1]*0.5);
plot(G(:,1),G(:,2),'k.')
% contour of constant exceedence
plot(XRng,YUp,'r-','linewidth',1.5);
plot(XRng,YDw,'r-','linewidth',1.5,'handlevisibility','off');
% contour of constant density
plot(XRng,ZUpHT,'g-','linewidth',1.5);
plot(XRng,ZDwHT,'g-','linewidth',1.5,'handlevisibility','off');
%threshold
plot(Thr*[1,1],ylim,'b--','linewidth',2);
legend('HT simulation','Data','Contour Cnst Exc','Contour Cnst Dns','HT Threshold')

