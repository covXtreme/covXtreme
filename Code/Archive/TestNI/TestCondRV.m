
clc; clf; clear;

n=3e3;

muZ=[0,0];
rho=0.9;  %correlation between Y,Z
muY=[0,0];


%% Monte Carlo
Z=randn(n,2)+muZ;  %independent in direction
[mZ,J]=max(Z,[],2);

%Conditional simulation
EY=muY+rho*(Z-muZ);
VarY=1-rho.^2;
StDevY=sqrt(VarY);
Y=randn(n,2)*StDevY+EY;
cY=Y(sub2ind([n,2],(1:n)',J)); %conditional
mY=max(Y,[],2); %omni max


%% NI

X=linspace(-4,10,100)';
dZ=diff(X(1:2));

clear PZ PY
%% marginal Z
PZ(:,1)=normcdf(X,muZ(1),1);
PZ(:,2)=normcdf(X,muZ(2),1);
PZ(:,3)=prod(PZ,2);

%% conditional Y | Z

%FW=normcdf(X-muY-rho.*(shiftdim(X,-2)-muZ),0,StDevY);
FW=normcdf(X-rho.*(shiftdim(X,-2)-muZ),muY,StDevY);

%bin pdf
fZ=normpdf(X,muZ,1);
fZ=permute(fZ,[3,2,1]);
%fZ=fZ./sum(fZ,2);

%omni pdf
fOZ=[normpdf(X,muZ(1),1)./normcdf(X,muZ(1),1).*PZ(:,3),normpdf(X,muZ(2),1)./normcdf(X,muZ(2),1).*PZ(:,3)];
fOZ=permute(fOZ,[3,2,1]);
%fOZ=fOZ./sum(fOZ,3);

%conditional Y
PY(:,1:2)=normcdf(X,muY,1);
% %PY(:,1:2)=PY(:,1:2)./max(PY(:,1:2));
%FW=normcdf(X-muY-rho.*(ZRng-muZ),0,sqrt(VarY));;
PY(:,3)=prod(PY,2);


%PY(:,4)=sum(sum(FW.*fZ.*fOZ,2),3);
PY(:,4)=sum(sum(FW.*fOZ,2),3)*dZ;
PY(:,4)=PY(:,4);
%PY(:,3)=prod(PY,2);

%% Plot
clf;
subplot(1,2,1)
plot(Z(:,1),Y(:,1),'.','markersize',10)
hold on
plot(Z(:,2),Y(:,2),'.','markersize',10)
hold on
plot(mZ,cY,'ko','markersize',5)
xlabel('Z')
ylabel('Y|Z')
legend('1','2','Mx','location','northwest')
hold on 
tx=xlim;
ty=ylim;

subplot(2,2,2)
plot(sort(Z),linspace(0,1,n),'r-','linewidth',2)
hold on
plot(X,PZ(:,1:2),'c--','linewidth',2)
plot(sort(mZ),linspace(0,1,n),'m-','linewidth',2)
plot(X,PZ(:,3),'b--','linewidth',2)
legend('MC 1','MC 2','NI 1', 'NI 2','MC Omni','NI Omni','location','northwest')
xlim(tx)
title('P(Z)')

subplot(2,2,4)
plot(sort(Y),linspace(0,1,n),'r-','linewidth',2)
hold on
plot(X,PY(:,1:2),'c--','linewidth',2)
plot(sort(cY),linspace(0,1,n),'g-','linewidth',2)
plot(X,PY(:,4),'k--','linewidth',2)
% plot(sort(mY),linspace(0,1,n),'m-','linewidth',2)
% plot(X,PY(:,3),'b--','linewidth',2)
legend('MC 1','MC 2','NI 1', 'NI 2','MC Cond','NI Cond','location','northwest')
xlim(ty)
title('P(Z|Y)')
