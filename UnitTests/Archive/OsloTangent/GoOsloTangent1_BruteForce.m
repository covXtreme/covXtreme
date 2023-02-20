%% Original Oslo Tangent Code
% Lots of loops, search for solution by brute force

clf; hold on; 
clear;
i=sqrt(-1);

%Non-exceedance for quantile
Tau=[0.90]';nTau=size(Tau,1);

%Set up interesting sample
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

%Parameters for fit
CrrLng=100; %correlation length (set to something big for convex)
Cnt=[-1 1]; %sensible centre point (not critical)
dBeta=0.1;
dTheta=0.1;
dR=0.1;
MxmR=10;%this needs to scale with the problem

%Set up search to find local tangent
Beta=(0:dBeta:2*pi+dBeta)';nB=size(Beta,1);
Theta=(0:dTheta:2*pi+dTheta)';nT=size(Theta,1);
R=(dR:dR:MxmR)';nR=size(R,1); 

C=nan(nB,1); %contour length from Cnt
A=nan(nB,1); %tangent angle on contour
tic
for iB=1:nB
    bt=Beta(iB);
    OptPrb=nan(nR,1);
    OptTh=nan(nR,1);
    for iR=1:nR
        T=[Cnt(1)+R(iR)*cos(Beta(iB)) Cnt(2)+R(iR)*sin(Beta(iB))];
        dX=[X(:,1)-T(1) X(:,2)-T(2)]; 
        Dst=sqrt(dX(:,1).^2+dX(:,2).^2);
        Ang=angle(dX(:,1)+i*dX(:,2));
        Prb=nan(nT,1);
        for iT=1:nT
            th=Theta(iT);
            Prj=Dst.*cos(th-Ang); %projection onto the normal to the tangent
            Wgh=normpdf(Dst.*sin(th-Ang),0,CrrLng); %weighted projection onto the tangent
            Prb(iT)=sum(Wgh(Prj<=0))/sum(Wgh); %weighted probability "weighted quantile"
        end;
        
        t1=Theta(Prb==max(Prb));t1=t1(1);
        OptTh(iR)=t1;
        OptPrb(iR)=max(Prb); %tangent is the one which has max non-exceedence prb
        
    end;
    t=abs(OptPrb-Tau);
    t2=find(t==min(t));t2=t2(1);
    C(iB)=R(t2);
    A(iB)=OptTh(t2);
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