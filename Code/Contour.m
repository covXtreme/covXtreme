% Copyright © [2023] Shell Global Solutions International B.V. All Rights Reserved.
% SPDX-License-Identifier: Apache-2.0

classdef Contour
    % Environmental contour class
    
    properties
        XRng %[nPnt x (nBin+1) x Cnt.nLvln], conditioned values for contour
        XY; %[nMth x 1] cell array
        %In case of Exc and Hus
        %XY{iMth} is [nPnt x 2 x (nBin+1)x Cnt.nLvl x Cnt.nAsc] defining contour lines
        % in case of HTDNs
        %XY{iMth} is [(nBin+1)x Cnt.nLvl x Cnt.nAsc] cell with sub-elements [2 x nPnt] defining contour
        %in this case nPnt varies for each contour bin, return period and associated variable.
        Mth;  %[nMth x 1], (cell array) of contour methods used
        MthLabel; %Method label for plotting different contour methods
        nMth; %1 x 1, number of contouring methods
        nBin; %1 x 1, number of covariate bins
        nLvl; %1 x 1, number of contour levels chosen
        nAsc; %1 x 1 number of associated variables
        Sml; %structure importance sampled simulation under the model
        PltOn=false; %flag to switch on explanatory diagrams
        LvlOrg; %[(nBin+1) x nLvl x 2], contour level on orginal scale of conditioned variable
        nGrd; %Number of grid locations in each of x and y when x-y domain needs discretisation
        nPnt; %[1 x 1], how many points from which to draw contour
        SmtWdtC; %Smoothing width for the huseby contour C function
        BndWdtScl; %Band width scale for HTDns contour
        nSml; %Number of Simulations
    end
    
    methods  
        
        function Cnt=Contour(HT,Mrg,OptCnt)
            % Contour CONSTRUCTOR function
            %
            %% INPUTS
            % HT Heffernan and Tawn class object
            % Mrg Marginal model class object
            % OptCnt Options for Contour estimation class object
            %
            %% OUTPUTS
            % Cnt Contours class object
            
            %% Check inputs
            if nargin==0
                return;
            end
            if ~isa(Mrg,'MarginalModel')
                error('Mrg should be a nDmn x 1 Marginal Model');
            end
            
            if ~isa(HT,'HeffernanTawn')
                error('HT should be a HeffernanTawn class');
            end
            
            if ~isa(OptCnt,'OptionsContours')
                error('OptCnt should be an OptionsContours class');
            end
            
            %% Populate properties of Cnt from OptCnt and elsewhere
            % From OptCnt
            Cnt.nSml=OptCnt.nSml;
            Cnt.nGrd=OptCnt.nGrd;
            Cnt.nPnt=OptCnt.nPnt;
            Cnt.SmtWdtC=OptCnt.SmtWdtC;
            Cnt.BndWdtScl=OptCnt.BndWdtScl;
            Cnt.Mth=OptCnt.Mth;
            Cnt.nMth=numel(Cnt.Mth);
            % From elsewhere
            Cnt.nLvl=numel(Mrg(1).RtrPrd);
            Cnt.MthLabel=cell(Cnt.nMth,1);
            Cnt.nBin=Mrg(1).Bn.nBin;
            Cnt.nAsc=HT.nDmn-1;
            
            %% Compute conditional return value
            fprintf('Simulating under Heffernan and Tawn model using importance sampling\n');
            nSmlBin=ceil(Cnt.nSml./HT.nBin);
            Cnt.Sml=HT.SimulateIS(Mrg,nSmlBin);
            
        end %Contour constructor
        
        function Cnt = makeContours(Cnt,Mrg,HT,A,nA)
            % Estimate contours
            %
            %% INPUTS
            % Cnt Contours class containing settings/options for contours
            % Mrg Marginal model class output
            % HT  Heffernan and Tawn class output
            % A   bin Allocation vector
            % nA  nNumber of possible bin Allocations
            %
            %% OUTPUTS
            % Cnt Contours class object
            
            if nargin<=3
                A=Cnt.Sml.A; %use original bins
                nA=Cnt.nBin;
            end
            
            fprintf(1,'Computing Contour Curves\n')
            if isempty(Cnt.LvlOrg)
                Cnt=Cnt.GetLockPoint(HT);
            end
            % bins + omni
            if nA>1
                nBinsplsOmni = nA+1;
            else
                nBinsplsOmni = 1;
            end
            
            %% Get Contour levels from Marginal and Conditional analysis            
            Cnt.XRng=NaN(Cnt.nPnt,nBinsplsOmni, Cnt.nLvl);  %X range for contours
            
            %% Compute lock point for each bin            
            %Initialise XVal
            for iBin=1:nBinsplsOmni
                for iLvl=1:Cnt.nLvl
                    Cnt.XRng(:,iBin,iLvl)=linspace(min(Mrg(1).Y),Cnt.LvlOrg(iBin,iLvl,1),Cnt.nPnt);
                end
            end
            %Initialise YVal
            Cnt.XY = cell(Cnt.nMth,1);
  
            %% Constant Exceedence contour (invariant to margins);
            if any(strcmp(Cnt.Mth,'Exc'))
                fprintf(1,'Exceedence Contour\n')
                tic
                Cnt=Exceedence(Cnt,A,nA);
                toc
            end
            
            %% Huseby Contour Method  (Original Scale);
            if any(strcmp(Cnt.Mth,'HusOld'))
                fprintf(1,'Huseby Contour\n')
                tic
                Cnt=HusebyContour(Cnt,A,nA);
                toc
            end
            
            %% Huseby Cleaned Contour Method  (Original Scale);
            if any(strcmp(Cnt.Mth,'Hus'))
                fprintf(1,'Huseby Contour Cleaned\n')
                tic
                Cnt=HusebyContourCleaned(Cnt,A,nA);
                toc
            end
            
            %% Constant Hefferenan and Tawn Density contour  (must be on StandardMargins margins);
            if any(strcmp(Cnt.Mth,'HTDns'))
                fprintf(1,'Heffernan and Tawn Density contour\n')
                tic
                Cnt=HTDensity(Cnt,A,nA);
                toc            
            end
            
        end %makeCounters
        
        function Cnt=GetLockPoint(Cnt,HT)
            % Estimate lock point location
            %% INPUTS 
            % Cnt Contours class containing settings/options for contours
            % HT  Heffernan and Tawn class output
            %% OUTPUTS
            % Cnt Contours class populated with the lock point
            
            sOrig=HT.GetSummary('quantile',0.5);
            
            %Contours just computed at the bin
            Cnt.LvlOrg= permute(sOrig,[1,3,2]);  %range in X for contour on original scale
        end %GetLockPoint
        
        function Cnt=Exceedence(Cnt,A,nA)
            % Constant Exceedence probability contour
            %% INPUTS
            % Cnt Contours class containing settings/options for contours
            % A   bin Allocation vector
            % nA  nNumber of possible bin Allocations
            
            if nargin<=1
                A=Cnt.Sml.A;
                nA=Cnt.nBin;
            end
            
            IMth=strcmp(Cnt.Mth,'Exc');
            Cnt.MthLabel{IMth}='JntExc';

            % # bins + omni
            if nA>1
                nBinsplsOmni = nA+1;
            else
                nBinsplsOmni = 1;
            end
            
            % initialise
            Cnt.XY{IMth}=NaN(2.*Cnt.nPnt,2,nBinsplsOmni,Cnt.nLvl,Cnt.nAsc);
            
            % Loop over contour levels
            for iQ=1:Cnt.nLvl %contour levels
                fprintf('#');
                for iB=1:nBinsplsOmni
                    % pre-allocate
                    YDw = nan(Cnt.nPnt,Cnt.nAsc);
                    YUp = nan(Cnt.nPnt,Cnt.nAsc);
                    % loop over associated variables
                    for iAsc = 1:Cnt.nAsc
                        %lock point not defined for this bin/return level combo
                        if any(isnan(Cnt.LvlOrg(iB,iQ,:)))
                            continue
                        end
                        if iB==nBinsplsOmni %omni case
                            SmlX = Cnt.Sml.Org(:,1);
                            logfogX = Cnt.Sml.logfog(:,1)+Cnt.Sml.W;
                            SmlY = Cnt.Sml.Org(:,iAsc+1);
                            logfogY = Cnt.Sml.logfog(:,iAsc+1);
                        else
                            %pick out simulations from each bin in turn
                            if ~any(A==iB)
                                continue
                            end
                            SmlX = Cnt.Sml.Org(A==iB,1);
                            logfogX = Cnt.Sml.logfog(A==iB,1)+Cnt.Sml.W(A==iB);
                            SmlY = Cnt.Sml.Org(A==iB,iAsc+1);
                            logfogY = Cnt.Sml.logfog(A==iB,iAsc+1);
                        end
                        % compute overall f(.)/g(.)
                        fog = exp(logfogX+logfogY);
                        
                        % compute P(X>x_{loc},Y>y_{loc})
                        I_gtLoc_Up = (SmlX>Cnt.LvlOrg(iB,iQ,1))&...
                            (SmlY>Cnt.LvlOrg(iB,iQ,iAsc+1));
                        I_gtLoc_Dwn = (SmlX>Cnt.LvlOrg(iB,iQ,1))&...
                            (SmlY<=Cnt.LvlOrg(iB,iQ,iAsc+1));
                        p_Loc_Up = sum(fog(I_gtLoc_Up))./sum(fog);
                        p_Loc_Dwn = sum(fog(I_gtLoc_Dwn))./sum(fog);
                        
                        % grid in main and associated
                        x_grd = Cnt.XRng(:,iB,iQ);
                        y_grd = linspace(min(SmlY),max((SmlY)),Cnt.nGrd)';
                        
                        % P(X>x)
                        P_x = sum((x_grd>SmlX').*fog',2)./sum(fog);
                        % P(Y>y|X>x)
                        P_ygx = nan(Cnt.nGrd,Cnt.nPnt);
                        for iG = 1:Cnt.nPnt
                            I_x = SmlX(:,1)>x_grd(iG);
                            if any(I_x)
                                P_ygx(:,iG) = Cnt.WeightedCDF(y_grd,SmlY(I_x),fog(I_x));
                            end
                            
                        end %iG
                        
                        % find y_Dw closest to this probability level
                        tJnt=(P_ygx).*(1-P_x'); %joint Prb(Y>y,X>x) for a sequence of values of x
                        YDw(:,iAsc)=Cnt.exceedanceSmooth(tJnt,y_grd,p_Loc_Dwn); % using linear interpolation (produces smoother output)
                        
                        % find y_Up closest to this probability level                        
                        tJnt=(1-P_ygx).*(1-P_x'); %joint Prb(Y<=y,X>x) for a sequence of values of x
                        YUp(:,iAsc)=Cnt.exceedanceSmooth(tJnt,y_grd,p_Loc_Up); % using linear interpolation (produces smoother output)
                        
                    end %iAsc
                    
                    %store the 2 parts of the contour (up and down) together
                    for iAsc=1:Cnt.nAsc
                        Cnt.XY{IMth}(:,2,iB,iQ,iAsc) = [permute(YDw(:,iAsc),[1,3,4,5,2]);flipud(permute(YUp(:,iAsc),[1,3,4,5,2]))]; % last dimention is Associated
                        Cnt.XY{IMth}(:,1,iB,iQ,iAsc) = [Cnt.XRng(:,iB,iQ);flipud(Cnt.XRng(:,iB,iQ))]; 
                    end %iAsc
                end %iB
            end %iQ
            fprintf('\n');
            
            
            %% diagram
            if Cnt.PltOn
                clf;
                
                rectangle('position',[Cnt.LvlOrg(iB,iQ,1),Cnt.LvlOrg(iB,iQ,2),20,20],'curvature',0,'edgecolor','none','facecolor','b');
                hold on;
                rectangle('position',[Cnt.LvlOrg(iB,iQ,1),0,20,Cnt.LvlOrg(iB,iQ,2)],'curvature',0,'edgecolor','none','facecolor','r')
                plot(Cnt.Sml.Org(:,1),Cnt.Sml.Org(:,2),'k.','markersize',5)
                hold on
                plot(Cnt.XY{1}(1:Cnt.nPnt,1,end,iQ,1),Cnt.XY{1}(1:Cnt.nPnt,2,end,iQ,1),'r-','linewidth',2)
                plot(Cnt.XY{1}((Cnt.nPnt+1):end,1,end,iQ,1),Cnt.XY{1}((Cnt.nPnt+1):end,2,end,iQ,1),'b-','linewidth',2)
                
                plot(Cnt.LvlOrg(iB,iQ,1),Cnt.LvlOrg(iB,iQ,2),'go','markersize',10,'linewidth',2)
                axis tight
                plot([Cnt.LvlOrg(iB,iQ,1),Cnt.LvlOrg(iB,iQ,1)*2],Cnt.LvlOrg(iB,iQ,2)*[1,1],'k--','linewidth',2)
                plot(Cnt.LvlOrg(iB,iQ,1)*[1,1],ylim,'k--','linewidth',2)
                
                xlim([min(Cnt.Sml.Org(:,1)),max(Cnt.Sml.Org(:,1))]);
                ylim([min(Cnt.Sml.Org(:,2)),max(Cnt.Sml.Org(:,2))]);
                
                savePics('ExcDiagram')
            end
            
        end %Exceedence
        
        function Cnt=HTDensity(Cnt,A,nA)
            % Contours of constant HT density
            %
            %% INPUTS
            % Cnt Contours class containing settings/options for contours
            % A   bin Allocation vector
            % nA  nNumber of possible bin Allocations

            if nargin<=1
                A=Cnt.Sml.A;
                nA=Cnt.nBin;
            end
            
            %% Parse inputs
            % method
            IMth=strcmp(Cnt.Mth,'HTDns');
            Cnt.MthLabel{IMth}='Heffernan Tawn density';
            % # bins + omni
            if nA>1
                nBinsplsOmni = nA+1;
            else
                nBinsplsOmni = 1;
            end
            
            % initialise
            Cnt.XY{IMth}=cell(nBinsplsOmni,Cnt.nLvl,Cnt.nAsc);
            
            % loop over associated variables
            for iAsc = 1:Cnt.nAsc
                % loop over bins
                for iBin = 1:nBinsplsOmni
                    fprintf('#')
                    Lck=shiftdim(Cnt.LvlOrg(iBin,:,[1,iAsc+1]),1);
                    
                    if all(isnan(Lck(:,1))) || all(isnan(Lck(:,2)))%bad case
                        continue
                    end
                    
                    %% use contour function to get iso-line.
                    %setup gridded to bin density into
                    I_notinf=~isinf(Cnt.Sml.Org(:,iAsc+1));
                    
                    
                    YLimFct = 1.3;
                    if iBin>nA %omni case!!
                        fog=exp(sum(Cnt.Sml.logfog(:,[1,iAsc+1]),2)).*Cnt.Sml.W;
                        Lmt=squeeze(Cnt.LvlOrg(end,end,[1,iAsc+1]));
                        Edgx=linspace(min(Cnt.Sml.Org(I_notinf,1)),Lmt(1)*YLimFct,Cnt.nGrd+1);
                        Edgy=linspace(min(Cnt.Sml.Org(I_notinf,iAsc+1)),max(Cnt.Sml.Org(I_notinf,iAsc+1))*YLimFct,Cnt.nGrd+1);
                    else
                        if ~any(A==iBin) % dealing with restricted domain cases with there is no data in a bin
                            continue
                        end
                        fog=exp(sum(Cnt.Sml.logfog(:,[1,iAsc+1]),2)).*Cnt.Sml.W;
                        Lmt=squeeze(Cnt.LvlOrg(iBin,end,[1,iAsc+1]));
                        Edgx=linspace(min(Cnt.Sml.Org(A==iBin & I_notinf,1)),Lmt(1)*YLimFct,Cnt.nGrd+1);
                        Edgy=linspace(min(Cnt.Sml.Org(A==iBin & I_notinf,iAsc+1)),max(Cnt.Sml.Org(A==iBin & I_notinf,iAsc+1))*YLimFct,Cnt.nGrd+1);
                        
                    end
                    
                    Grdx=(Edgx(1:end-1)+Edgx(2:end))/2;
                    Grdy=(Edgy(1:end-1)+Edgy(2:end))/2;
                    
                    [GX,GY]=ndgrid(Grdx,Grdy);
                    G=[GX(:),GY(:)];
                    
                    Ax=discretize(Cnt.Sml.Org(:,1),Edgx);  %returns indices of the bins (Grdx) that SmlOrg falls into
                    Ay=discretize(Cnt.Sml.Org(:,iAsc+1),Edgy);
                    
                    GrdInd=sub2ind([Cnt.nGrd Cnt.nGrd],Ax,Ay);
                    
                    %% find gridded density estimate for current bin
                    if iBin>nA %omni case!!
                        I=~isnan(GrdInd); %in current bin and not a nan
                    else
                        I=  (~isnan(GrdInd)) & (A==iBin); %in current bin and not a nan
                    end
                    %% which bin is lock point in?
                    
                    Lx=discretize(Lck(:,1),Edgx);
                    Ly=discretize(Lck(:,2),Edgy);
                    L=sub2ind([Cnt.nGrd Cnt.nGrd],Lx,Ly);
                    
                    BW=Cnt.BndWdtScl.*(range(Cnt.Sml.Org(I,[1,iAsc+1])));
                    f=reshape(ksdensity(Cnt.Sml.Org(I,[1,iAsc+1]),G,'Weights',fog(I),'Bandwidth',BW),[Cnt.nGrd Cnt.nGrd]);
                    %% find lock point density value for each return level
                    IL=find(~isnan(L));
                    
                    if isempty(IL)
                        continue;
                    end
                    Lvl=f(L(IL));
                    if numel(Lvl)==1
                        Lvl=Lvl.*[1,1];
                    end
                    % compute contour using low level function
                    C =contourc(Grdx,Grdy,f',Lvl);
                    
                    [I,J]=ismember(C(1,:),Lvl);   %identify different contour segments inside C by locating the separatprs 'Lvl' (see help file on contour...odd output style)
                    Ind=J(I); %which contour does each segment belong to
                    Count=cumsum(C(2,I)+1);  %cumulative sum of the number points in each contour segment
                    
                    c=0;
                    for i=1:numel(Ind) %loop over line segments
                        iLvl=IL(Ind(i)); %contour level;
                        tI=c+2:Count(i);
                        c=Count(i);
                        %assign line segement to right level;
                        if isempty(Cnt.XY{IMth}{iBin,iLvl,iAsc})
                            Cnt.XY{IMth}{iBin,iLvl,iAsc}=C(:,tI);
                        else
                            Cnt.XY{IMth}{iBin,iLvl,iAsc}=[Cnt.XY{IMth}{iBin,iLvl,iAsc},[NaN;NaN],C(:,tI)];
                        end
                    end %i
                end %iBin
            end%iAsc
            fprintf('\n');
            
            
        end %HTDensity
        
        function Cnt=HusebyContour(Cnt,A,nA)
            % Drawing Huseby contours from a bivariate sample
            %
            %% INPUTS
            % Cnt Contours class containing settings/options for contours
            % A   bin Allocation vector
            % nA  nNumber of possible bin Allocations
            
            if nargin<=1
                A=Cnt.Sml.A;
            end
            
            %% Parse inputs
            % method
            IMth=strcmp(Cnt.Mth,'HusOld');
            Cnt.MthLabel{IMth}='TangExc';
            % # bins + omni
            if nA>1
                nBinsplsOmni = nA+1;
            else
                nBinsplsOmni = 1;
            end
            % initialise
            Cnt.XY{IMth}=NaN(Cnt.nPnt,2,nBinsplsOmni,Cnt.nLvl,Cnt.nAsc);
            
            % set of angles for contour eval.
            angles = linspace(-pi,pi,Cnt.nPnt); %set of angles
            
            % loop over associated variables
            for iAsc = 1:Cnt.nAsc
                
                % loop over bins
                for iBin = 1:nBinsplsOmni
                    fprintf('#')
                    
                    J=[1,iAsc+1];
                    if (iBin==nBinsplsOmni)&&(nBinsplsOmni>1)
                        x=Cnt.Sml.Org(:,J);
                        fog=exp(sum(Cnt.Sml.logfog(:,J),2)).*Cnt.Sml.W;
                    else
                        if ~any(A==iBin) % check whether data exist in the bin
                            continue
                        else
                            tI=A==iBin;
                            x=Cnt.Sml.Org(tI,J);
                            fog=exp(sum(Cnt.Sml.logfog(tI,J),2)).*Cnt.Sml.W(tI);
                        end
                    end
                    
                    % E[X]
                    IGood=~isinf(x(:,2));
                    E_X = (sum(x(IGood,:).*fog(IGood),1)./sum(fog(IGood)))';
                    % Cov[X]
                    Cov_X = (x(IGood,:)-E_X')'*(fog(IGood).*(x(IGood,:)-E_X'))./sum(fog(IGood));
                    % chol(Cov[X])
                    cholCov_X = cholcov(Cov_X);                   
                    z = (x-E_X')/cholCov_X;
                    
                    % compute prob. associated with lock point
                    IGd=find(~isnan(Cnt.LvlOrg(iBin,:,1)) & Cnt.LvlOrg(iBin,:,1)>E_X(1));
                    
                    p_c = sum(bsxfun(@gt,x(:,1),Cnt.LvlOrg(iBin,IGd,1)).*fog,1)./sum(fog,1);
                                       
                    % set up knots for the IS density
                    Z_knt = linspace(min(z(:))-1e-2,max(z(:))+1e-2,20000)';
                    
                    
                    %% Calculate the function C_theta
                    C_theta=NaN(Cnt.nPnt, numel(p_c)); %preallocate
                    for iAng=1:length(angles)
                        % X_p = c*X
                        Z_p = z(:,1)*cos(angles(iAng)) + z(:,2)*sin(angles(iAng));
                        % Compute importance sampled density
                        P_Zp = 1-Cnt.WeightedCDF(Z_knt,Z_p,fog);
                        % find closest point to req. probability
                        [~,I_min] = min(abs(P_Zp-p_c),[],1);
                        % get C(\theta)
                        C_theta(iAng,:) = Z_knt(I_min);
                    end %iAng
                    C_theta=movmean(C_theta,5); %smooth Huseby to get rid of 'spikey corners'
                    
                    %% Calculate the intersection points
                    
                    for iRtr=1:numel(p_c) %loop over return periods
                        for iAng =1:Cnt.nPnt-1  %loop over points
                            % compute contour points for standardised
                            % variables
                            tXY = nan(1,2);
                            tXY(1) = (sin(angles(iAng+1))*C_theta(iAng, iRtr)-sin(angles(iAng))*C_theta(iAng+1,iRtr))./...
                                (sin(angles(iAng+1))*cos(angles(iAng))-sin(angles(iAng))*cos(angles(iAng+1)));
                            tXY(2) = (-cos(angles(iAng+1))*C_theta(iAng, iRtr)+cos(angles(iAng))*C_theta(iAng+1, iRtr))./...
                                (sin(angles(iAng+1))*cos(angles(iAng))-sin(angles(iAng))*cos(angles(iAng+1)));
                            % transform back to original scale
                            Cnt.XY{IMth}(iAng,:,iBin,IGd(iRtr),iAsc) = tXY*cholCov_X + E_X';

                        end %iAng
                    end %iRtr
                end %iBin loop over bins
                
            end %loop over associated variables
            fprintf('\n')
            
        end % HusebyContour
        
        function Cnt=HusebyContourCleaned(Cnt,A,nA)
            %Huseby cleaned contours
            %
            %% INPUTS
            % Cnt contour object
            % A   optional alternative bin allocation vector
            % nA  integer scalar : largest value bin indicator (might be larger than largest value in A, if there are empty bins)
            
            %Use existing covariate bins if non specified
            if nargin<=1
                A=Cnt.Sml.A;
            end
            
            %% Parse inputs
            % method
            IMth=strcmp(Cnt.Mth,'Hus');
            Cnt.MthLabel{IMth}='TangExc';
            
            % # bins + omni
            if nA>1
                nBinsplsOmni = nA+1;
            else
                nBinsplsOmni = 1;
            end
            
            % initialise
            Cnt.XY{IMth}=NaN(Cnt.nPnt,2,nBinsplsOmni,Cnt.nLvl,Cnt.nAsc);
            
            % set of angles for contour evaluation
            angles = linspace(-pi,pi,Cnt.nPnt); %set of angles
            
            % loop over associated variables
            for iAsc = 1:Cnt.nAsc

                % loop over bins
                for iBin = 1:nBinsplsOmni
                    fprintf('#')
                    
                    %% Access the appropriate importance sampling "fog" ratio
                    J=[1,iAsc+1];
                    if (iBin==nBinsplsOmni)&&(nBinsplsOmni>1)
                        x=Cnt.Sml.Org(:,J);
                        fog=exp(sum(Cnt.Sml.logfog(:,J),2)).*Cnt.Sml.W;
                    else
                        if ~any(A==iBin) % check whether data exist in the bin
                            continue
                        else
                            tI=A==iBin;
                            x=Cnt.Sml.Org(tI,J);
                            fog=exp(sum(Cnt.Sml.logfog(tI,J),2)).*Cnt.Sml.W(tI);
                        end
                    end
                    
                    %% Sphere the data to make contour estimation smoother
                    % E[X]
                    IGood=~isinf(x(:,2));
                    E_X = (sum(x(IGood,:).*fog(IGood),1)./sum(fog(IGood)))';
                    % Cov[X]
                    Cov_X = (x(IGood,:)-E_X')'*(fog(IGood).*(x(IGood,:)-E_X'))./sum(fog(IGood));
                    % chol(Cov[X])
                    cholCov_X = cholcov(Cov_X);
                    z = (x-E_X')/cholCov_X;
                                        
                    %% Compute probability associated with lock point
                    % Identify good data
                    IGd=find(~isnan(Cnt.LvlOrg(iBin,:,1)) & Cnt.LvlOrg(iBin,:,1)>E_X(1));
                    p_c = sum(bsxfun(@gt,x(:,1),Cnt.LvlOrg(iBin,IGd,1)).*fog,1)./sum(fog,1);
                                       
                    %% Locations for calculation of density from importance
                    %sampling
                    Z_knt = linspace(min(z(:))-1e-2,max(z(:))+1e-2,20000)';
                    
                    %% Calculate the Huseby et al. function C_theta
                    C_theta=NaN(Cnt.nPnt, numel(p_c));
                    for iAng=1:length(angles)
                        Z_p = z(:,1)*cos(angles(iAng)) + z(:,2)*sin(angles(iAng));
                        % Compute importance sampled density
                        P_Zp = 1-Cnt.WeightedCDF(Z_knt,Z_p,fog);
                        % find closest point to req. probability
                        [~,I_min] = min(abs(P_Zp-p_c),[],1);
                        % get C(\theta)
                        C_theta(iAng,:) = Z_knt(I_min);
                    end %iAng
                    
                    %% Apply a moving median smoother
                    C_theta=movmean(C_theta,Cnt.SmtWdtC); %smooth Huseby to get rid of 'spikey corners'
                                      
                    %% Calculate the contour
                    % See Huseby, A. B., Vanem, E., Natvig, B. (2015). 
                    % Alternative environmental contours for structural reliability analysis.
                    % Structural Safety 54 32–45
                    for iRtr=1:numel(p_c) %loop over return periods
                        
                       %% Loop over angles, finding intersection points
                        for iAng =1:Cnt.nPnt-1  %loop over points
                            % compute contour points for standardised
                            % variables
                            tXY = nan(1,2);
                            tXY(1) = (sin(angles(iAng+1))*C_theta(iAng, iRtr)-sin(angles(iAng))*C_theta(iAng+1,iRtr))./...
                                (sin(angles(iAng+1))*cos(angles(iAng))-sin(angles(iAng))*cos(angles(iAng+1)));
                            tXY(2) = (-cos(angles(iAng+1))*C_theta(iAng, iRtr)+cos(angles(iAng))*C_theta(iAng+1, iRtr))./...
                                (sin(angles(iAng+1))*cos(angles(iAng))-sin(angles(iAng))*cos(angles(iAng+1)));
                            % transform back to original scale
                            Cnt.XY{IMth}(iAng,:,iBin,IGd(iRtr),iAsc) = tXY*cholCov_X + E_X';
                        end %iAng
                        
                        %% Call to CleanHuseby
                        % The original Huseby contour contains uninteresting "bow ties"
                        % Dropping non-physical angles, padding array with
                        % following NaN
                        Cnt.XY{IMth}(:,:,iBin,IGd(iRtr),iAsc)=Cnt.CleanHuseby(Cnt.XY{IMth}(:,:,iBin,IGd(iRtr),iAsc));    

                    end %iRtr loop over return periods
                    
                end %iBin loop over bins
                
            end %iAsc loop over associated variables
            fprintf('\n');
            
        end % HusebyContourCleaned
        
        function Plot(Cnt,Mrg)
            
            %Use different line styles to distinguish return periods
            LinStl={'-';'--';'-.';':'};
            
            %% plot omni contour
            figure(1);
            clf;
            for iAsc=1:Cnt.nAsc
                subplot(1,Cnt.nAsc,iAsc);
                
                plot(Mrg(1).Y,Mrg(iAsc+1).Y,'k.','markersize',10,'handlevisibility','off')
                hold on
                grid on
                box on
                
                xlabel(sprintf('%s: Conditioning',Mrg(1).RspLbl))
                ylabel(sprintf('%s: Conditioned',Mrg(iAsc+1).RspLbl))
                C=lines(Cnt.nMth);
                
                for iQ=1:Cnt.nLvl %return period
                    if iQ<=4
                        tLS=LinStl{iQ};
                    else
                        tLS='-';
                    end
                    for iCnt=1:Cnt.nMth  %method
                        %General plotting code:
                        if strcmp(Cnt.Mth(iCnt),'HTDns')
                            if any(Cnt.XY{iCnt}{end,iQ,iAsc}(:))
                                plot(Cnt.XY{iCnt}{end,iQ,iAsc}(1,:),Cnt.XY{iCnt}{end,iQ,iAsc}(2,:),'linestyle',tLS,'color',C(iCnt,:),'linewidth',2);
                                hold on
                            end
                        else
                            plot(Cnt.XY{iCnt}(:,1,end,iQ,iAsc),Cnt.XY{iCnt}(:,2,end,iQ,iAsc),'linestyle',tLS,'color',C(iCnt,:),'linewidth',2);
                            hold on
                        end
                    end %iCnt
                    plot(Cnt.LvlOrg(end,:,1),Cnt.LvlOrg(end,:,iAsc+1),'go','markersize',8, 'LineWidth',2);
                end %iQ
                
                title(sprintf('%s|%s: Omni-covariate',Mrg(iAsc+1).RspLbl,Mrg(1).RspLbl))
                legend(Cnt.Mth,'location','best')
            end
            
            savePics('Figures/Stg5_Contour_1_Omni');
            
            %% plot contour by bin
            if Cnt.nBin >1  %only produce if nBin>1 (else this is identical to above
                
                for iAsc=1:Cnt.nAsc
                    figure(iAsc+1);
                    clf;
                    nPlt1=ceil(sqrt(Cnt.nBin)); %max size nPlt x nPlt
                    nPlt2=ceil(Cnt.nBin./nPlt1);
                    
                    
                    for iBin=1:Cnt.nBin
                        subplot(nPlt2,nPlt1,iBin)
                        hold on
                        grid on
                        box on
                        J=Mrg(1).Bn.A==iBin;
                        plot(Mrg(1).Y(J),Mrg(iAsc+1).Y(J),'k.','markersize',10,'handlevisibility','off')
                        xlabel(sprintf('%s: Conditioning',Mrg(1).RspLbl))
                        ylabel(sprintf('%s: Conditioned',Mrg(iAsc+1).RspLbl))
                        C=lines(Cnt.nMth);
                        
                        for iQ=1:Cnt.nLvl %return period
                            if iQ<=4
                                tLS=LinStl{iQ};
                            else
                                tLS='-';
                            end
                            for iCnt=1:Cnt.nMth %method
                                if strcmp(Cnt.Mth(iCnt),'HTDns')
                                    if any(Cnt.XY{iCnt}{iBin,iQ,iAsc}(:))
                                        plot(Cnt.XY{iCnt}{iBin,iQ,iAsc}(1,:),Cnt.XY{iCnt}{iBin,iQ,iAsc}(2,:),'linestyle',tLS,'color',C(iCnt,:),'linewidth',2);
                                        hold on
                                    end
                                else
                                    plot(Cnt.XY{iCnt}(:,1,iBin,iQ,iAsc),Cnt.XY{iCnt}(1:end,2,iBin,iQ,iAsc),'linestyle',tLS,'color',C(iCnt,:),'linewidth',2);
                                    hold on
                                end
                            end %iCnt
                            
                        end %IQ
                        plot(Cnt.LvlOrg(iBin,:,1),Cnt.LvlOrg(iBin,:,iAsc+1),'go','markersize',8, 'LineWidth',2);  %BUG: NaN in Cnt.LvlOrg(iBin,:,iAsc+1)
                        
                        if iBin==1
                            title(sprintf('Bin %s',Mrg(1).Bn.BinLbl{iBin}))
                            legend(Cnt.Mth,'location','best')
                        else
                            title(sprintf('Bin %s',Mrg(1).Bn.BinLbl{iBin}))
                        end
                        
                    end %iBin
                    savePics(sprintf('Figures/Stg5_Contour_2_Binned_%s',Mrg(iAsc+1).RspLbl));
                end %iAsc
            end
            
        end %plot
    end%methods
     
    methods(Static)
        function P=WeightedCDF(X,Y,W)
            %Efficiently Compute Weighted CDF at requested locations. Used in computing empirical CDF and importance
            %sampled CDF
            % Direct copy of the ReturnValue repository WeightedCDF.m function
            %
            % P(X<=x) = int (X<=Y).*W  over samples
            %
            %% INPUTS
            %X  nX x Sz(1) x Sz(2)  which points to compute CDF at
            %Y: nRls x Sz(1) x Sz(2)  observations to derive CDF for
            %W: nRls x Sz(1) x Sz(2)  weights to use in compute cdf. (importance weights)
            %q: [nQ x 1]  also provide quantile levels to compute
            %
            % Sz need to match for binary singleton expansion
            %% OUTPUTS
            % P = probabilities  nX x Sz(1) x Sz(2)
            % Yq = quantiles of data.
            %
            %Reference https://en.wikipedia.org/wiki/Importance_sampling           
            
            %% Input defaults
            %interpolate knots
            validateattributes(X,{'numeric'},{},'IS_CDF','X',1)
            
            if isvector(X)
                X=X(:); %make sure X is nX x 1;
            end
            nX=size(X,1);
            %observations
            validateattributes(Y,{'numeric'},{},'IS_CDF','Y',2)
            SzY=size(Y);
            nRls=SzY(1);
            SzX=size(X);
            if nargin>=3
                SzW=size(W);
                validateattributes(Y,{'numeric'},{'nrows',nRls},'IS_CDF','Y',2)
                W_On=true;
                Sz=NaN(1,max([numel(SzW),numel(SzY),numel(SzX)])-1);
            else
                W_On=false;
                Sz=NaN(1,max([numel(SzY),numel(SzX)])-1);
            end
            
            %% match dimension sizes
            
            Sz(1:(numel(SzY)-1))=SzY(2:end);
            for i=2:numel(SzX)
                Sz(i-1)=max(Sz(i-1),SzX(i));
            end %i
            if W_On
                SzW=size(W);
                for i=2:numel(SzW)
                    Sz(i-1)=max(Sz(i-1),SzW(i));
                end
            end
            SzOther=prod(Sz); %number of passive dimensions
            
            %repmat to match sizes
            X=reshape(X.*ones([1,Sz]),[],SzOther);
            Y=reshape(Y.*ones([1,Sz]),[],SzOther);
            if W_On
                W=reshape(W.*ones([1,Sz]),[],SzOther);
            end
            %% Sort data
            [~,I]=sort([X;Y],1); %place together edge vector with data and sort
            Ijoint=I+(nX+nRls).*(0:SzOther-1); %2D index
            IndEdg=ismember(I,(1:nX)'); %edge index where the edges lie in sorted matrix
            
            if W_On
                sW=[zeros(nX,SzOther);W]; %compute integral for all weights use zero weight for edges
            else
                sW=[zeros(nX,SzOther);ones(nRls,SzOther)];
            end
            cW=cumsum(sW(Ijoint),1); %sum weights to give integral at the data
            P=reshape(cW(IndEdg),nX,SzOther)./cW(end,:); %pull out integral at bin edges and renormalise
            P(:,cW(end,:)==0)=1; %deal with bad cases
            
            P=reshape(P,[nX,Sz]);
            
            P(P>1)=1;
            P(P<0)=0;
            
        end %WeightedCDF

        function XNew=CleanHuseby(XOld)
            % Clean up Huseby contours around "cusps"
            %
            % Algorithm works by searching for and eliminating points that cannot be
            % accessed from the contour median without crossing a line segment (between
            % adjacent points on the contour)
            %% INPUTS
            % XOld uncleaned Huseby contour
            %% OUTPUTS
            % Xnew cleaned Huseby contour
            
            % Test mode
            if nargin==0
                Om=(0:0.1:2*pi)';
                nOm=size(Om,1);
                XOld=[cos(Om)+0.1*randn(nOm,1) sin(Om)+0.1*randn(nOm,1)];
            end
            
            % Drop missing values
            X=XOld(sum(~isnan(XOld),2)==2,:);
            nX=size(X,1);
            
            % Kep is an indicator to good points (to keep)
            Kep=ones(nX,1);
            nKep=sum(Kep);
            for iX=1:nX
                
                % Drop a single point
                Kep=Contour.CleanHusebyDropOne(X,Kep);
                if sum(Kep)==nKep-1
                    nKep=nKep-1;
                else
                    % Stop if no points for dropping
                    break;
                end
                
            end %iX
            
            % Return the cleaned contour, padded with NaN if necessary
            XNew=nan(size(XOld,1),2);
            XNew(1:sum(Kep==1),:)=X(Kep==1,:);
            
        end %CleanHuseby
        
        function Kep=CleanHusebyDropOne(X,Kep)
            % function Kep=CleanHusebyDropOne(X,Kep);
            %
            % Subfunction for CleanHuseby
            %
            % The algorithm drops a single point from a Huseby contour corresponding to a "bow-tie" point
            %
            % The point to drop is identified iteratively, by searching around the contour from start to end
            % Once a point to drop is identified, the algorithm terminates.
            %
            % The criterion for a "potential point to drop" is as follows
            % - the line between the mean of the data and a candidate point in drawn
            % - we check whether this line intersects any line segment connecting two consecutive 
            %   points on the contour
            % - if an intersection is found, the point identified is a "bow-tie" point
           
            
            % Set verbose mode
            Vrb=0;
            
            % Identifier
            Idn=(1:size(X,1))';
            
            % Keep only good points
            X=X(Kep==1,:);
            Idn=Idn(Kep==1);
            nX=size(X,1);
            
            % Find the step vector between points on the approximate contour
            Dlt=nan(nX,2);
            for iX=1:nX
                if iX<nX
                    Dlt(iX,:)=[X(iX+1,1)-X(iX,1) X(iX+1,2)-X(iX,2)];
                else
                    Dlt(nX,:)=[X(1,1)-X(nX,1) X(1,2)-X(nX,2)];
                end
            end %iX
            
            % Define an origin
            Org=mean(X)';
            
            for iX=1:nX
                
                % Candidate point to drop
                C=X(iX,:)'-Org;
                
                if Vrb==1
                    clf; hold on;
                    plot(X(:,1),X(:,2),'ko-');
                    plot(Org(1),Org(2),'b*-');
                    plot(X(iX,1),X(iX,2),'r*');
                end
                
                %Loop over all line segments formed by consecutive points on the
                %contour
                tPrm=nan(nX,2);
                for jX=1:nX
                    % Start point
                    A=X(jX,:)'-Org;
                    if jX<nX
                        % Increment
                        D=Dlt(jX,:)';
                    else
                        % Increment
                        D=Dlt(1,:)';
                    end
                    if jX~=iX-1 && jX~=iX % Cannot match with self
                        tPrm(jX,:)=[-D(1) C(1);-D(2) C(2)]\A; % This solves for the intersection between two vectors
                    end
                end %jX
                % Find any line segments that make the point a valid point to drop
                if sum(tPrm(:,1)>0 & tPrm(:,1)<1 & tPrm(:,2)>0 & tPrm(:,2)<1)>0
                    t=find(tPrm(:,1)>0 & tPrm(:,1)<=1 & tPrm(:,2)>0 & tPrm(:,2)<1);
                    if Vrb==1
                        plot(X(t,1),X(t,2),'g*');
                        fprintf(1,'%g Hit',Idn(iX));
                        fprintf(1,' %g',Idn(t==1));
                        fprintf(1,'\n');
                    end
                    % Drop the offending point in Kep
                    Kep(Idn(iX))=0;
                    % Break the loop because you've found a point to drop
                    break;
                end
                
            end %iX
            
        end %CleanHusebyDropOne
        
        function Yq=exceedanceSmooth(X,Y,Xq)
            % Smooth the exceedance contour
            %% INPUTS
            % X    X locations, possibly with repeated values (ie values are not unique)  
            % Y    Y locations corresponding to Y
            % Xq   Points at which interpolated value is desired
            % Yq:  Interpolated values 
            Yq=nan(size(X,2),1);
            for ti=1:size(X,2)
                [tX,iIndUnq]=unique(X(:,ti),'last');
                tY=Y(iIndUnq);
                t2Y=interp1(tX,tY,Xq);
                Yq(ti)=t2Y;
            end %ti
        end %exceedanceSmooth
        
    end %methods (Static)
    
end %end class Contour


