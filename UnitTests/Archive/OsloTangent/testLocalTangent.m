

CntPnt=[-1,1];
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
KrnWdt=0.5; %correlation length (set to something big for convex)
Cnt=[-1 1]; %sensible centre point (not critical)
dBeta=0.1;
dTheta=0.1;
dR=0.1;

Prb=0.9;
nGrd=30;
[XVal,YVal] = LocalTangent(X,Prb,nGrd,CntPnt,KrnWdt)

% 
%        
clf; hold on;
plot(X(:,1),X(:,2),'k.');
plot(XVal,YVal,'r-','linewidth',2);
CI=convhull(XVal,YVal);
plot(XVal(CI),YVal(CI),'b--','linewidth',1);
%pAxsLmt; pDfl; pDatStm;
%pGI('LocalTangentContour',2);
