function PLOGL=gplikenonstationary(PARAMS,Dat,BinAllocation,L)
%PLOGL=gplikenonstationary(PARAMS,Dat,BinAllocation,L)
%Piecewise constant gplike
%INPUTS:
% - PARAMS = [(constant) xi, sig_1,..., sig_nBins]
% - Dat = (nDat x 1) vector GP data (=exceedences - thr)
% - BinAllocation = (nDat x 1) vector of bin indices
% - L = penalty imposing smoothness in sigmas as move between bins

%OUTPUT:
% Scalar Penalised neg. log likelihood

nBins=max(BinAllocation);
if size(PARAMS,1)~=(nBins+1)
   error('Input dimension mismatch: dmn input1 (#parameters+1) =?= max input3 (#bins)') 
end    
if size(Dat,1) ~= size(BinAllocation,1)
    error('Input dimension mismatch: dmn Dat =?= dmn BinAllocation')
end

if nargin<=3
   L=0;  %unpenalised 
end
Xi=PARAMS(1);
Sig=exp(PARAMS(2:end)); 

NLOGL=0;
for iBin=1:nBins %Additive like over bins
    I=BinAllocation==iBin;
    if any(I)
        NLOGL=NLOGL+gplike([Xi,Sig(iBin)],Dat(I));
    end
end

PLOGL=NLOGL+L.*sum((Sig-mean(Sig)).^2);




