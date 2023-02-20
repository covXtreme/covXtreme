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


%% Trun data into polar coordinates
[A,R]=cart2pol(X(:,1)-Cnt(1),X(:,2)-Cnt(2));
[A,I]=sort(A); R=R(I);

nA=20; %number of angles
AngBin=linspace(-pi,pi,nA)'; %angle bins

R0=quantile(R,Tau)*ones(nA,1);  %inital radii
G=normpdf(R.*sin(AngBin'-A),0,CrrLng);



clf;
subplot(1,2,1)
plot(A,R,'k.');
hold on
plot(A(I),R(I),'g.');
axis tight
hold on;
plot(AngBin,R0,'r--')
title('Polar Coordinates')

subplot(2,2,2)
plot(X(:,1),X(:,2),'k.');
hold on
plot(X(I,1),X(I,2),'g.');
axis tight
hold on;
CX=Cnt(1)+R0.*cos(AngBin);
CY=Cnt(2)+R0.*sin(AngBin);
plot(CX,CY,'r--')
title('Cartesian Coordinates')
axis equal

subplot(2,2,4)
plot(AngBin,1-sum(IExc)./n,'k-')
title('Exceedence Probability')