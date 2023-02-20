%% Oslo Tangent Method
% Speed up tried: vectorisation (didn't improve speed much)

clf; hold on; 
clear;
%% Non-exceedance for quantile
Tau=[0.7]';
CrrLng=0.5; %correlation length (set to something big for convex)
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

%% Parameters for fit

dBeta=0.1;
dTheta=0.1;
dR=0.1;
MxmR=10;%this needs to scale with the problem

%% Set up search to find local tangent
Beta=(0:dBeta:2*pi+dBeta)';nB=size(Beta,1);

nT=200;
Theta=linspace(-pi/4,pi/4,nT)';
R=(dR:dR:MxmR)';nR=size(R,1); 

%% Preallocate
nTau=size(Tau,1);
C=nan(nB,1); %contour length from Cnt

figure(1);
clf;
subplot(2,1,1)
plot([R(1),R(end)],Tau*[1,1],'k-','linewidth',2)
hold on
Col=hsv(nB);

tic
for iB=1:nB
    DmT=Beta(iB)-Theta;
    
    %% Turn data into polar coordinates
    T=shiftdim([Cnt(1)+bsxfun(@times,R,cos(Beta(iB))) Cnt(2)+bsxfun(@times,R,sin(Beta(iB)))],-1);
    [Ang,Dst]=cart2pol(bsxfun(@minus,X(:,1),T(:,:,1)),bsxfun(@minus,X(:,2),T(:,:,2)));
    
    %% weighted projection onto the normal to the tangent
    ThetmAng=bsxfun(@minus,shiftdim(DmT,-2),Ang);
    Prj=bsxfun(@times,Dst,cos(ThetmAng)); 
    Wgh=normpdf(bsxfun(@times,Dst,sin(ThetmAng)),0,CrrLng); %weighted projection onto the tangent
    
    %% weighted probability "weighted quantile"
    Prb=sum(bsxfun(@times,Wgh,Prj<=0))./sum(Wgh); %sum over data; 
    [OptPrb,I]=max(Prb,[],3); %location of min over theta
    
    
    if 0
    clf
    imagesc(DmT*180/pi,R,squeeze(Prb))
    axis xy;
    ylabel('Radius')
    xlabel('Angle')
%     clf;
    subplot(2,1,1)
    plot(R,OptPrb,'-','color',Col(iB,:));
    hold on
    subplot(2,1,2)
    plot(R,OptTheta,'-','color',Col(iB,:));
    end
            
    t=abs(OptPrb-Tau); %find closest match to desired Tau
        
    [~,J]=min(t);          
    C(iB)=R(J);               
end;
toc       

figure(2);
clf; hold on;
plot(X(:,1),X(:,2),'k.');
CX=Cnt(1)+C.*cos(Beta);
CY=Cnt(2)+C.*sin(Beta);
CI=convhull(CX,CY);
plot(CX,CY,'r-','linewidth',2);
plot(CX(CI),CY(CI),'b--','linewidth',1);
%pAxsLmt; pDfl; pDatStm;
%pGI('LocalTangentContour',2);