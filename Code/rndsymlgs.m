function F=rndsymlgs(alp,nD,n)
%generate random values from nD dimension symmetric logistics function with
%Frechet margins
%
%Algorithm 1.1 Alec Stephenson, Simulating Multivariate Extreme Value
%Distributions of Logistic Type, Extremes 6, 49-59, 2003
%INPUT
%alp 1 x 1 or n x 1 vector of the dependency parameter
%nD  is number of dimensions
%n is the number of realisations
%OUTPUT
%F nD x n random values on frechet margins

if numel(alp)>1
    n=numel(alp);
    nalp=n;
else
    nalp=1;
end

%step 1)  get T from unit exponentials
W=exprnd(1,[nD,n]);
T=W./sum(W);

%find mixture probabilities from recurrence relationship
p=NaN(nD,nD,nalp);
for iA=1:nalp
    p(1,1,iA)=1;
    for iD=2:nD
        p(iD,1,iA)=gamma(iD-alp(iA))./(gamma(iD).*gamma(1-alp(iA)));
        for jD=2:(iD-1)
            t1=(iD-1-alp(iA)*jD).*p(iD-1,jD)+alp(iA).*(jD-1).*p(iD-1,jD-1);
            p(iD,jD,iA)=t1./(iD-1);
        end
        p(iD,iD,iA)=alp(iA).^(iD-1);
    end
end
 
P=cumsum(shiftdim(p(end,:,:),1)); %nD x nalp  %cumulative probability

%step2) Find k
U=rand(1,n);
k=sum(U>P,1)+1;
%step3) get Z from gamma(k,1);
Z=gamrnd(k,1,[1,n]);
%step4) find F on frechet margins
F=1./(Z.*(T.^(alp)));

