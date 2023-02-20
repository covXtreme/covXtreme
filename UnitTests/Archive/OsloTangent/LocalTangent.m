function [XVal,YVal] = LocalTangent(X,Tau,nGrd,CntPnt,KrnWdt)
%INPUTS:
% - X data n x 2
% - Prb 1 x 1 non-exceedance for quantile
% - nGrd 1 x 1 number of angles beta to draw contour at
% - KrnWdt=100; %kernel width (set to something big for convex)
% - CntPnt=[-1 1]; %sensible centre point of data (not critical)
% OUTPUTS:
% - [Cnt.XVal,Cnt.YVal]

%Parameters for fit
MxmR=10;   %TODO: automate this! needs to scale with the problem, rng x^2 + rngY^2 / 2
dR=MxmR/100;

%Set up search to find local tangent
Beta=linspace(0,2*pi,nGrd)';nB=size(Beta,1);
nT=nGrd;
Theta=linspace(-pi/2,pi/2,nT)'; %only look in half plane around Beta to avoid big jumps in tangent location
R=(0:dR:MxmR)';

C=nan(nB,1); %contour length from CntPnt

for iB=1:nB %loop over angle around centre point
    BtmT=Beta(iB)-Theta; 
    T=shiftdim([CntPnt(1)+R.*cos(Beta(iB)) CntPnt(2)+R.*sin(Beta(iB))],-1);
    [Ang,Dst]=cart2pol(bsxfun(@minus,X(:,1),T(:,:,1)),bsxfun(@minus,X(:,2),T(:,:,2)));
    
    ThetamAng=bsxfun(@minus,shiftdim(BtmT,-2),Ang);
    Prj=bsxfun(@times,Dst,cos(ThetamAng)); %projection onto the normal to the tangent
    Wgh=normpdf(bsxfun(@times,Dst,sin(ThetamAng)),0,KrnWdt); %kernel-weighted projection onto the tangent: give higher weight to data closer to the line
    
    Prb=sum(bsxfun(@times,Wgh,Prj<=0))./sum(Wgh);
    t=abs(max(Prb,[],3)-Tau); 
    [~,J]=min(t);
    C(iB)=R(J);
     
end % end loop over angle around center point

%Output conour
XVal=CntPnt(1)+C.*cos(Beta);
YVal=CntPnt(2)+C.*sin(Beta);


