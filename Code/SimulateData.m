function Dat=SimulateData(Mdl)
%Simulate model data from a directional model

%% Check input parameter consistency
if Mdl.nBin ~= size(Mdl.DrcEdg,1) || Mdl.nBin ~= size(Mdl.Rat,1)
   error('Dimension of DrcEdg and/or Rate inconsistent with nBins')
end
if (Mdl.nDmn > 1) %nDmn =2
    if Mdl.nBin ~= size(Mdl.MM(2).Scl,1)
        error('Dimension of Scl inconsistent with nBins in dimension 2')
    end
    if size(Mdl.MM(2).Shp,1) > 1
        error('GP shape assumed contant, Mdl.MM(2).Shp should be a scalar')
    end
else %nDmn=1
    if Mdl.nBin ~= size(Mdl.MM(1).Scl,1)
        error('Dimension of Scl inconsistent with nBins in dimension 1')
    end
    if size(Mdl.MM(1).Shp,1) > 1
        error('GP shape assumed contant, Mdl.MM(1).Shp should be a scalar')
    end
end

if numel(Mdl.Jnt.Alp) > 1
    if size(Mdl.Jnt.Alp,1) ~= Mdl.nBin
        error('Alpha should be [1 x 1] or [Mdl.nBin x 1]')
    end
end


%Simulate covariate (directions) with specified rate (Mdl.Rat)
RatCdf=cumsum(Mdl.Rat)./sum(Mdl.Rat); %get rate cdf 
A=discretize(rand(Mdl.nDat,1),[0;RatCdf]);  %simulate unif values on [0,1] and allocate to drc bins 1:nBins by probabilities specified by RatCdf

%% Simulate Joint tail behaviour on Uniform Margins
if Mdl.nDmn>1
    switch Mdl.Jnt.Mth
        case 'ASL' %Assymetric Logistic on Frechet margins
            %check input parameters
            if (max(Mdl.Jnt.Alp)>1) || (min(Mdl.Jnt.Alp)<0)
               error('Please provide a dependency parameter Alpha in [0,1] for ASL dependency model') 
            end
            if any(Mdl.Jnt.Theta>1) || any(Mdl.Jnt.Theta <0)
               error('Please provide weighting parameters Theta in [0,1] for ASL dependency model') 
            end
            %reshape for input to nD rnsasymlgs function
            Alp=[NaN,NaN,Mdl.Jnt.Alp];  %dependency parameter(s) - 3 parameters in nD=2;
            Th={Mdl.Jnt.Theta(1),Mdl.Jnt.Theta(2),[1-Mdl.Jnt.Theta(1),1-Mdl.Jnt.Theta(2)]}; %weighting parameter
            F=rndasymlgs(Alp,Th,Mdl.nDmn,Mdl.nDat)';  %simulated on frechet scale
            U=exp(-(1./F));% transform to uniform scale
        case 'LGS' %Logistic on Frechet margins
            if (max(Mdl.Jnt.Alp)>1) || (min(Mdl.Jnt.Alp)<0)
               error('Please provide a dependency parameter Alpha in [0,1] for LGS dependency model') 
            end
            F=rndsymlgs(Mdl.Jnt.Alp(A),Mdl.nDmn,Mdl.nDat)';
            U=exp(-(1./F));   % transform to uniform scale
        case 'MVN' %Multivariate normal on Normal margins
            if max(Mdl.Jnt.Rho)>1 || min(Mdl.Jnt.Rho) <-1
                error('Please provide a dependency parameter Rho in [0,1] for MVN dependency model')
            elseif any(Mdl.Jnt.Rho< 0) || any(Mdl.Jnt.Rho>= 1)
                error('Please provide a dependency parameter Rho in [0,1] for MVN dependency model - later conditional H&T model is set up only for positive dependence. To model negative dependence, simulate with Rho in [0,1] and multiply by -1') 
            end
            nRho=size(Mdl.Jnt.Rho,2);
            if nRho ==1
                Sig=[1,Mdl.Jnt.Rho;Mdl.Jnt.Rho,1];  %rho
                Z=mvnrnd(zeros(Mdl.nDmn,1),Sig,Mdl.nDat);
            end
            if nRho > 1
                tabRho=tabulate(A);
                Z=NaN(Mdl.nDat,Mdl.nDmn);
				NumberRho=tabRho(:,2)';
                NoRho=size(NumberRho,2);
                for iRho=1:NoRho
                    SigChange=[1,Mdl.Jnt.Rho(iRho);Mdl.Jnt.Rho(iRho),1];  %rho
                    Z(A==iRho,:)=mvnrnd(zeros(Mdl.nDmn,1),SigChange,NumberRho(iRho));
                end    
            end
            U=normcdf(Z); % transform to uniform scale
        otherwise
            error('Joint tail method not recognised, should be one of ASL, LGS or MVN');
    end        
else %if 1D don't need to use joint model
    U=rand(Mdl.nDat,1);
end

%% Transform data to Generalised Pareto margins

%Initialise empty response(s) and covariate
Dat.Y=NaN(Mdl.nDat,Mdl.nDmn);  %response
Dat.X=NaN(Mdl.nDat,1);  %covariates


Edg=[Mdl.DrcEdg(end)-360;Mdl.DrcEdg]; %bin edges (add final bin for wrapping)
BinSze=Edg(2:end)-Edg(1:end-1);  %bin width

%Simulate covariate (directions) with specified rate (Mdl.Rat) - which has been determined for the directional model 
Dat.X=mod(rand(Mdl.nDat,1).*BinSze(A)+Edg(A),360);  %For each storm assigned to a given bin, simulate a direction/covariate value within that bin  

for i=1:Mdl.nDmn   %loop over dimension                     
    %Obtain response Y: transform joint uniform data to GP margins with
    %scale and threshold of covariate bin (associated with direction X)
    Dat.Y(:,i)=gpinv(U(:,i),Mdl.MM(i).Shp,Mdl.MM(i).Scl(A),Mdl.MM(i).Thr(A));
end

Dat.IsPrd=true;
Dat.RspLbl={'Response','Associated'};
Dat.CvrLbl={'Direction'};

%% Plot 
if ~exist('Figures','dir')
   mkdir('Figures') 
end

%% marginal
figure(1);
clf;
for i=1:Mdl.nDmn
    subplot(Mdl.nDmn,1,i)
    plot( Dat.X,Dat.Y(:,i),'k.');
    xlabel(Dat.CvrLbl)
    set(gca,'xtick',0:45:360,'xlim',[0,360])
    ylabel(Dat.RspLbl{i})
end
savePics('Figures/Stg1_Data_Simulated_Margins')

%% joint
%TODO generalised to nD
if Mdl.nDmn>1
    figure(2);
    clf;   
    plot( Dat.Y(:,1),Dat.Y(:,2),'k.');   
    xlabel(Dat.RspLbl{1})
    ylabel(Dat.RspLbl{2})
    title('Original Margins')
 
end
savePics('Figures/Stg1_Data_Simulated_Joint')

