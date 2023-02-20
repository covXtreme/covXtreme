function fval=objfun(R,Tau,X,Beta,Cnt,Mth,CrrLng)
%R    1 x 1 radius  (objective parameter)
%Tau  1 x 1
%Ang  nDat x
%Beta value 1 x 1
%Cnt 2 x1 center point
%Mth   Method 'Convex','Kernel'
%CrrLng  1 x 1 correlation length  (for Gaussian Kernel)


%Output
%fval  function value


%% Now we do the line search over theta
nT=50;
Theta=linspace(0,2*pi,nT);
ThetmAng=bsxfun(@minus,Theta,Ang);
Prj=bsxfun(@times,Dst,cos(ThetmAng)); %projection onto the normal to the tangent
switch Mth %smoothing method
    case 'Convex' %none (gives convex curve)
        Prb=sum(Prj<=0)./size(Prj,1);
    case 'Kernel' %guassian kernel
        Wgh=normpdf(bsxfun(@times,Dst,sin(ThetmAng)),0,CrrLng); %weighted projection onto the tangent
        Prb=sum(Wgh.*(Prj<=0))./sum(Wgh); %weighted probability "weighted quantile"    
end
OptPrb=max(Prb,[],2); %location of min over theta

fval=abs(OptPrb-Tau); %find closest match to desired Tau

