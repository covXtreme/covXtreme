%% Olso Tangent Method
%Speed up tried: solving optimisation problem to find optimal tangent angle
%and radius. This didn't work --> unstable

clf; hold on; 
clear;
%% Non-exceedance for quantile
Tau=[0.90]';
CrrLng=1; %correlation length (set to something big for convex)
Cnt=[-1 1]; %sensible centre point (not critical)

%% Set up interesting sample
n1=700*2;
n2=200*2;
n3=100*2;
n=n1+n2+n3;
Sigma1=2*[1 0.95;0.95 1];
X1=randn(n1,2)*sqrtm(Sigma1);
Sigma2=[0.5 0;0 2];
X2=randn(n2,2)*sqrtm(Sigma2);
Sigma3=0.5*[1 0;0 1];
X3=randn(n3,2)*sqrtm(Sigma3);
X=[X1+ones(n1,1)*[0 1];X2+ones(n2,1)*[-2 2];X3+ones(n3,1)*[-1 1]];

%% Preallocate
dBeta=0.1;
Beta=(0:dBeta:2*pi+dBeta)';nB=size(Beta,1);

nTau=size(Tau,1);
C=nan(nB,1); %contour length from Cnt

Mth='Kernel';

opts=optimset('display','off','largescale','off');
tic
for iB=1:nB
    
    T=[Cnt(1)+R.*cos(Beta), Cnt(2)+R.*sin(Beta)];
    [Ang,Dst]=cart2pol(X(:,1)-T(:,1)',X(:,2),T(:,1)');

    
    
    C(iB)=fminunc(@(R)objfun(R,Tau,X,Beta(iB),Cnt,Mth,CrrLng),R0,opts);           
end;
toc       

clf; hold on;
plot(X(:,1),X(:,2),'k.');
CX=Cnt(1)+C.*cos(Beta);
CY=Cnt(2)+C.*sin(Beta);
CI=convhull(CX,CY);
plot(CX,CY,'r-','linewidth',2);
plot(CX(CI),CY(CI),'b--','linewidth',1);
%pAxsLmt; pDfl; pDatStm;
%pGI('LocalTangentContour',2);