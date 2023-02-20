%% Development code: adaptation of "Oslo Tangent" method which allows local non-convexity
% Phil/Emma/David ,16.03.17
% Speed-ups attempted:
%%
% 1. Reduce the range of tangent angles we try by only looking at angles
% close to angle at neighbouring point. --> Local non-convexity means can
% get big turns around a corner in the contour - restricting angles tried
% doesn't work
%
% 2. Reduce range of radii searched. Again: can't base it on radius close
% to neighbouring optimal radius; can see big step changes. 
% Therefore tried a binary segmentation search for best radius (due to
% monotonicity of exceedance probability along radius, given that have
% chosen the tangent angle well). This is what is implemented here. 

% NOTE: code still deemed too slow to include in PPC. Can look @ Arne's
% code to see where he has saved computation time. 
%%

tic;
clf; hold on;
clear;
i=sqrt(-1);

%Non-exceedance for quantile
Tau=0.90;

%% Data
n1=300*2;
n2=100*2;
n3=50*2;
n=n1+n2+n3;
Sigma1=2*[1 0.95;0.95 1];
X1=randn(n1,2)*sqrtm(Sigma1);
Sigma2=[0.5 0;0 2];
X2=randn(n2,2)*sqrtm(Sigma2);
Sigma3=0.5*[1 0;0 1];
X3=randn(n3,2)*sqrtm(Sigma3);
X=[X1+ones(n1,1)*[0 1];X2+ones(n2,1)*[-2 2];X3+ones(n3,1)*[-1 1]];

%Parameters for fit
CrrLng=0.5; %correlation length (set to something big for convex)
Cnt=[-1.5 0.5]; %sensible centre point (not critical)
dBeta=0.01;
dTheta=0.01;
MxmR=10;%this needs to scale with the problem

%Set up search to find local tangent
Beta=(0:dBeta:2*pi)';nB=size(Beta,1);
Theta=(0:dTheta:2*pi)';nT=size(Theta,1);

C=nan(nB,1); %contour length from Cnt
A=nan(nB,1); %tangent angle on contour
for iB=1:nB; 
    
    Bet=Beta(iB);
    
    R=[0;MxmR/2;MxmR];
    Prb=[0;NaN;1];
    OptTht=[NaN;NaN;NaN];
    
    for iH=1:10;
        
        T=[Cnt(1)+R(2)*cos(Beta(iB)) Cnt(2)+R(2)*sin(Beta(iB))];
        dX=[X(:,1)-T(1) X(:,2)-T(2)];
        Dst=sqrt(dX(:,1).^2+dX(:,2).^2);
        Ang=angle(dX(:,1)+i*dX(:,2));
        LclPrb=nan(nT,1);
        
        for iT=1:nT;
            
            Tht=Theta(iT);
            Go=0;
            if iB==1;
                Go=1;
            else
                dTht=angle(cos(Tht-A(iB-1))+i*sin(Tht-A(iB-1)));
                if abs(dTht)<10*dTheta;
                    Go=1;
                end;
            end;
            if Go==1;
                Spr=10;
                Ltr=Dst.*sin(Tht-Ang);
                tInd=abs(Ltr)<Spr*CrrLng;
                tLtr=Ltr(tInd);
                tDst=Dst(tInd);
                tAng=Ang(tInd);
                tPrj=tDst.*cos(Tht-tAng); %projection onto the normal to the tangent
                if Spr<5;
                    tWgh=normpdf(tDst.*sin(Tht-tAng),0,CrrLng); %weighted projection onto the tangent
                else;
                    tWgh=tDst>0;
                end;
                LclPrb(iT)=sum(tWgh(tPrj<=0))/sum(tWgh); %weighted probability "weighted quantile"
            end;
        end;
        
        t1=Theta(LclPrb==nanmax(LclPrb));
        OptTht(2)=t1(1);
        
        tPrb=Prb;
        tR=R;
        tOptTht=OptTht;
        if Tau>max(LclPrb);
            Prb(1)=max(LclPrb);
            Prb(3)=tPrb(3);
            R(1)=tR(2);
            R(3)=tR(3);
            OptTht(1)=tOptTht(2);
            OptTht(3)=tOptTht(3);
        else;
            Prb(1)=tPrb(1);
            Prb(3)=max(LclPrb);
            R(1)=tR(1);
            R(3)=tR(2);
            OptTht(1)=tOptTht(1);
            OptTht(3)=tOptTht(2);
        end;
        
        R(2)=(R(1)+R(3))/2;
        
    end;
    
    C(iB)=R(1);
    A(iB)=OptTht(1);
    
end;

if 1;
    HlfWdt=5;
    Inc=0.2;
    %sC=pSmt([C(end-HlfWdt+1:end);C;C(1:HlfWdt)],HlfWdt,'MEAN','C');
    %sC=sC(HlfWdt+1:end-HlfWdt);
    sC=C;
    clf; hold on;
    CX=Cnt(1)+sC.*cos(Beta);
    CY=Cnt(2)+sC.*sin(Beta);
    plot(CX,CY,'r-','linewidth',2);
    DX1=CX-Inc.*cos(A-pi/2);
    DX2=CX+Inc.*cos(A-pi/2);
    DY1=CY-Inc.*sin(A-pi/2);
    DY2=CY+Inc.*sin(A-pi/2);
    for iB=1:nB;
        plot([DX1(iB);DX2(iB)],[DY1(iB);DY2(iB)],'g');
    end;
    plot(CX,CY,'r-','linewidth',2);
    plot(X(:,1),X(:,2),'k.');
end;
pAxsLmt; pDfl; pDatStm;
pGI('LocalTangentContour',2);


toc;