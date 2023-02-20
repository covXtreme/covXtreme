classdef GaussianExpected<HeffernanTawn
    %Heffernan and Tawn (2004) conditional extreme value model
    %Fitting the data using a Gaussian distribution
    %Expected Hessian 
    %Working with the Variance
            
    properties
        %SigmaTransform - if variance should be a one to one mapping
        %                - if standard deviation take the square root
        SigmaTransform=@(x) sqrt(x); 
        %SigmaTransform - starting solution for the fourth parameter of the
        %Heffernan and Tawn model
        VarStartSol=@(x) var(x,1)
    end
    
    methods        
        %includes those methods that are different for the different
        %inference methods
        %Gradient        
        function [G,H]=Gradient(HT,iP,X,Y,A,p,L)
            %[G,H]=Gradient(iP,X,Y,A,p,L,TwoPrmFlg)
            %compute gradients for the Hefferenan and Tawn model
            %INPUT
            % iP        either 1 for [alp and mu] or 2 for [beta and sigma]
            %-X         n x 1 conditioned value
            %-Y         n x 1 x (nDmn-1) conditioning value
            %-A         n x 1 x (nDmn-1) bin allocation
            %-p         nBin+3 x (nDmn-1) parameter values
            %-L         smoothness parameter (alpha)
            %
            %OUTPUT
            % G       p x (nDmn-1)
            % H       p x p x (nDmn-1)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Title:    Gradient and Expected Fisher information. Vectorized and
            %           suitable for splines.
            % Author:   Thijs Willems (graduate intern TU Delft)
            % Date:     17-10-2016
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin<=5
                L=0; %no penalty
            end
            
            
            [Alp,Bet,M,S]=HT.p2Prm(p);
            V=S.^2;
            
            if HT.nBin==1
                L=0; %switch off penalty
            end
            
            % Get paramters at data
            Alp_Dat = Alp(A(:,1),:);
            Bet_Dat = Bet(A(:,2),:);
            if HT.TwoParamFlag
                M_Dat=0;
                V_Dat=1;
            else
                M_Dat = M(A(:,3),:);
                V_Dat = V(A(:,4),:);  %B*Alp
            end
                        
            Xb=X.^Bet_Dat;
            MeanTerm=Alp_Dat.*X+(Xb).*M_Dat; %mean term of the distribution
                        
            %% First Derivatives
            switch iP
                case 1 %alpha and mu
                    G = NaN(sum(HT.nPrm([1,3])),HT.nDmn-1);
                    %dL/dalp
                    t1=-1./V_Dat.*X.^(1-2.*Bet_Dat).*(Y - MeanTerm);
                    for iDmn=1:HT.nDmn-1
                        %Penalty term
                        if HT.NonStat(1)
                            P=Alp(:,iDmn);
                        else
                            P=0;
                        end                            
                        G(1:HT.nPrm(1),iDmn) = accumarray(A(:,1),t1(:,iDmn),[HT.nPrm(1),1],@sum,0)+L.*P;
                    end
                    
                    %dL/dmu
                    if HT.TwoParamFlag
                        G(HT.nPrm(1)+1,:)=[];
                    else
                        t1=-(1./V_Dat).*X.^(-Bet_Dat).*(Y - MeanTerm);
                        for iDmn=1:HT.nDmn-1
                            %Penalty term
                            if HT.NonStat(3)
                                P=M(:,iDmn);
                            else
                                P=0;
                            end    
                            G(HT.nPrm(1)+1:end,iDmn) = accumarray(A(:,3),t1(:,iDmn),[HT.nPrm(3),1],@sum,0)+L.*P;
                        end                        
                    end
                case 2 %beta and var
                    G = NaN(sum(HT.nPrm([2,4])),HT.nDmn-1);
                    %dL/dbeta
                    t1= -(1./V_Dat).*X.^(-2.*Bet_Dat).*log(X).*(Y - Alp_Dat.*X).*(Y - MeanTerm) + log(X);
                    
                    for iDmn=1:HT.nDmn-1
                        %Penalty term
                        if HT.NonStat(2)
                            P=Bet(:,iDmn);
                        else
                            P=0;
                        end   
                        G(1:HT.nPrm(2),iDmn) = accumarray(A(:,2),t1(:,iDmn),[HT.nPrm(2),1],@sum,0)+L.*P;
                    end
                    
                    %dL/dV
                    if HT.TwoParamFlag
                        G(HT.nPrm(2)+1,:)=[];
                    else
                        t1=(1./(2.*V_Dat)).*(1-(1./V_Dat).*((Y - MeanTerm)./(X.^Bet_Dat)).^2);
                        for iDmn=1:HT.nDmn-1
                            %Penalty term
                            if HT.NonStat(4)
                                P=V(:,iDmn);
                            else
                                P=0;
                            end   
                            G(HT.nPrm(2)+1:end,iDmn) = accumarray(A(:,4),t1(:,iDmn),[HT.nPrm(4),1],@sum,0)+L.*P;
                        end
                    end
            end
                        
            %% Expected Fisher information matrix
            if nargout==2
                H=NaN(size(G,1),size(G,1),size(G,2));
                switch iP
                    case 1 %alpha and mu                      
                        for iDmn=1:HT.nDmn-1 %loop iDmn
                            %Penalty term
                            if HT.NonStat(1)
                                P=1;
                            else
                                P=0;
                            end   
                            t1=1./V_Dat(:,iDmn).*X.^(2-2*Bet_Dat(:,iDmn));
                            E11 = diag(accumarray(A(:,1),t1,[HT.nPrm(1),1],@sum,0)+L.*P);
                            if  HT.TwoParamFlag
                                H(:,:,iDmn)=E11;
                            else
                                t1=1./V_Dat(:,iDmn);
                                %Penalty term
                                if HT.NonStat(3)
                                    P=1;
                                else
                                    P=0;
                                end   
                                E33 =diag(accumarray(A(:,3),t1,[HT.nPrm(3),1],@sum,0)+L.*P);
                                t1=(1./V_Dat(:,iDmn)).*X.^(1-Bet_Dat(:,iDmn));
                                if HT.NonStat(1) &&  HT.NonStat(3) %both non stationary
                                    E13 = diag(accumarray(A(:,1),t1,[HT.nPrm(1),1],@sum,0)); % nBin  x nBin
                                elseif HT.NonStat(1)
                                    E13 = accumarray(A(:,1),t1,[HT.nPrm(1),1],@sum,0); % nBin  x 1
                                elseif HT.NonStat(3)
                                    E13 = accumarray(A(:,3),t1,[HT.nPrm(3),1],@sum,0)'; % 1 x nBin
                                else %both stationary
                                    E13 = sum(t1); % 1 x 1
                                end
                                H(:,:,iDmn)=[E11,E13;E13',E33]; %[nBin +1 x nBin +1]
                            end
                        end %iDmn
                    case 2 %beta and V
                        %d2L/DLbeta2
                        for iDmn=1:HT.nDmn-1 %loop iDmn
                            t1=(log(X).^2).*(2 + M_Dat(:,iDmn).^2./V_Dat(:,iDmn));
                            %Penalty term
                            if HT.NonStat(2)
                                P=1;
                            else
                                P=0;
                            end   
                            E22 = diag((accumarray(A(:,2),t1,[HT.nPrm(2),1],@sum,0))+L.*P);
                            
                            if  HT.TwoParamFlag
                                H(:,:,iDmn)=E22;
                            else
                                %Penalty term
                                if HT.NonStat(4)
                                    P=1;
                                else
                                    P=0;
                                end   
                                %d2L/DLV2
                                t1=1./(2.*V_Dat(:,iDmn).^2);
                                E44 =diag((accumarray(A(:,4),t1,[HT.nPrm(4),1],@sum,0))+L.*P);
                                %d2L/DLbetaV
                                t1=(1./V_Dat(:,iDmn)).*log(X);
                                if HT.NonStat(2) &&  HT.NonStat(4) %both non stationary
                                    E24 = diag(accumarray(A(:,2),t1,[HT.nPrm(2),1],@sum,0)); % nBin  x nBin
                                elseif HT.NonStat(2)
                                    E24 = (accumarray(A(:,2),t1,[HT.nPrm(2),1],@sum,0)); % nBin  x 1
                                elseif HT.NonStat(4)
                                    E24 = (accumarray(A(:,4),t1,[HT.nPrm(4),1],@sum,0))'; % 1 x nBin
                                else %both stationary
                                    E24 = (sum(t1)); % 1 x 1
                                end
                                
                                H(:,:,iDmn)=[E22,E24;E24',E44]; %[2 x 2]
                            end
                        end %iDmn
                end
            end
        end%Gradient        
    end %methods
end %classdef
