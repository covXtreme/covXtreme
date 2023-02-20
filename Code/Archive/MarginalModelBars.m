function Mdl=MarginalModel(AnalysisType,BasisType,Dat,DrcEdg,NEP)
%INPUTS
%BasisType: Constant,Spline,PPC
%Dat:  Dat.X  nObs x nCvr
%      Dat.Y  nObs x 1
%NEP: Quantile Regression non exceedence probability

%% Input checks
validateattributes(NEP,{'numeric'},{'<=',1,'>=',0});

%TODO check data
%TODO default inputs??


%% Threshold estimation
if NEP==0 %no need for threshold
    n=length(Dat.Y);
    Mdl.QR.Prm.PrdDat=zeros(n,1);
    Mdl.IExc=true(n,1);
else
    switch BasisType
        case 'Constant'
            Bss=Basis.Constant(AnalysisType,Dat.X);
        case 'Spline'
            Opts=Options;  %empty structure
            Bss=Basis.SplinenD(AnalysisType,Dat.X,Opts.Bss);
        case 'PPC'
            Bss=Basis.PPC(AnalysisType,Dat.X, DrcEdg );
        otherwise
            error('BasisType not recognised')
    end
        
    %% Quantile Regression
    Mdl.QR=Distribution.QuantileRegression(Bss,Dat.Y,NEP);
    Mdl.QR.Response='QuantileRegression';
    Mdl.QR.AnalysisType=AnalysisType;
    Mdl.QR.nDmn=length(AnalysisType);
    Mdl.QR=SimpleFit(Mdl.QR);
    
    %% Get Exceedences
    Mdl.IExc=Dat.Y>Mdl.QR.Prm.PrdDat;       
end

%%Need exceedence basis 
 switch BasisType
        case 'Constant'
            BssExc=Basis.Constant(AnalysisType,Dat.X(Mdl.IExc));
        case 'Spline'
            Opts=Options;  %empty structure
            BssExc=Basis.SplinenD(AnalysisType,Dat.X(Mdl.IExc),Opts.Bss);
        case 'PPC'
            BssExc=Basis.PPC(AnalysisType,Dat.X(Mdl.IExc), DrcEdg );
 end
    
%% Generalised Pareto
Mdl.GP=Distribution.GeneralisedPareto(BssExc,Dat.Y(Mdl.IExc)-Mdl.QR.Prm.PrdDat(Mdl.IExc));
Mdl.GP.BckAlg='Independent';
Mdl.GP.Response='GeneralisedPareto';
Mdl.GP.AnalysisType=AnalysisType;
Mdl.GP.nDmn=length(AnalysisType);
Mdl.GP=Fit(Mdl.GP);

%% Possion
%when we do poisson penalise rate not bin count!! because all bins have
%different volume
Mdl.Cnt=bincount(Mdl.GP.Prm(1).Bss);
Mdl.Pss=Distribution.Poisson(BssExc,Mdl.Cnt);
Mdl.Pss.Response='Poisson';
Mdl.Pss.AnalysisType=AnalysisType;
Mdl.Pss.nDmn=length(AnalysisType);
Mdl.Pss=SimpleFit(Mdl.Pss);
