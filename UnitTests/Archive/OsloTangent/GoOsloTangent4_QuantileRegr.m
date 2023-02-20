%% Oslo Tangent Method
% David trying to turn it into a Quntile regresion problem


clf; hold on; 
clear;
%% Non-exceedance for quantile
Tau=[0.8]';
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

nA=40; %number of angles
AngBin=linspace(-pi,pi,nA)'; %angle bins

R0=quantile(R,Tau)*ones(nA,1);  %inital radii

%% Gaussian Kernel Basis. Quantile Regression
fprintf('Quantile Regression')
tic
G=NaN(n,nA);
for iA=1:nA
    G(:,iA)=normpdf(A,AngBin(iA),CrrLng)+normpdf(A-2*pi,AngBin(iA),CrrLng)+normpdf(A+2*pi,AngBin(iA),CrrLng);
end
beta = ipqr (Tau, R, G);

CR=G*beta;

[CX,CY]=pol2cart(A,CR);
CX=CX+Cnt(1);
CY=CY+Cnt(2);
toc

%% Oslo Tangent
fprintf('Oslo Tangent\n')   
tic
G2=NaN(n,nA);
for iA=1:nA
    ThetmAng=AngBin(iA)-A;     
    G2(:,iA)=normpdf(bsxfun(@times,R,sin(ThetmAng)./cos(ThetmAng)),0,CrrLng); %weighted projection onto the tangent
end
beta2 = ipqr (Tau, R, G2);

CR2=G2*beta;

[CX2,CY2]=pol2cart(A,CR2);
CX2=CX2+Cnt(1);
CY2=CY2+Cnt(2);
toc


 
clf;
subplot(1,2,1)
plot(A,R,'k.');
axis tight
hold on;
plot(AngBin,R0,'r--')
title('Polar Coordinates')

subplot(1,2,2)
plot(X(:,1),X(:,2),'k.');
axis tight
hold on;
%CX=Cnt(1)+R0.*cos(AngBin);
%CY=Cnt(2)+R0.*sin(AngBin);
plot(CX,CY,'r--')
hold on
plot(CX2,CY2,'g--')
title('Cartesian Coordinates')
axis equal
% 
% subplot(2,2,4)
% plot(AngBin,1-sum(IExc)./n,'k-')
% title('Exceedence Probability')