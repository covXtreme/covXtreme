classdef DeltaLaplaceObserved<HeffernanTawn
    %Heffernan and Tawn (2004) conditional extreme value model
    %Generalised Gaussian (extra parameter compared to the Normal
    %distribution) - delta=1 => Laplace; delta=2=>Normal;
    %Working with the standard deviation
    %Type with delta=2
    %observed Hessian
    
    properties
        %SigmaTransform - if variance should be a one to one mapping
        %                - if standard deviation take the square root
        SigmaTransform=@(x) x;
        %SigmaTransform - starting solution for the fourth parameter of the
        %Heffernan and Tawn model
        VarStartSol=@(x) std(x,1);
     end%end properties
        
    methods
        
        %includes those methods that are different for the different
        %inference methods
        %Gradient
        %Likelihood
        
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
            
            if HT.nBin==1
                L=0; %switch off penalty
            end
            
            % Get paramters at data
            Alp_Dat = Alp(A(:,1),:);
            Bet_Dat = Bet(A(:,2),:);
            if HT.TwoParamFlag
                M_Dat=0;
                S_Dat=1;
            else
                M_Dat = M(A(:,3),:);
                S_Dat = S(A(:,4),:);  %B*Alp
            end
            Kappa=sqrt(gamma(1./HT.Delta)./gamma(3./HT.Delta));
            Xb=X.^Bet_Dat;  %X^b
            MeanTerm=Alp_Dat.*X+(Xb).*M_Dat; %mean term of the distribution
            StdTerm=Kappa.*(Xb).*S_Dat;%standard deviation of the distribtuion            
            R=Y-MeanTerm;
            
            
            %% First Derivatives
            switch iP
                case 1 %alpha and mu
                    G = NaN(sum(HT.nPrm([1,3])),HT.nDmn-1);
                    %dL/dalp
                    t1=-X.*HT.Delta.*R.*(abs(R).^(HT.Delta-2))./((StdTerm).^HT.Delta);
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
                        t1=-Xb.*HT.Delta.*R.*abs(R).^(HT.Delta-2)./((StdTerm).^HT.Delta);
                        
                        for iDmn=1:HT.nDmn-1
                            %Penalty term
                            if HT.NonStat(2)
                                P=M(:,iDmn);
                            else
                                P=0;
                            end   
                            G(HT.nPrm(1)+1:end,iDmn) = accumarray(A(:,3),t1(:,iDmn),[HT.nPrm(3),1],@sum,0)+L.*P;
                        end
                        
                    end
                case 2 %beta and sigma
                    G = NaN(sum(HT.nPrm([2,4])),HT.nDmn-1);
                    %dL/dbeta
                    t1=log(X)-HT.Delta.*R.*abs(R).^(HT.Delta-2).*M_Dat.*Xb.*log(X)./((StdTerm).^HT.Delta)...
                        -HT.Delta.*log(X).*abs(R).^(HT.Delta)./((StdTerm).^HT.Delta);
                    
                    for iDmn=1:HT.nDmn-1
                        %Penalty term
                        if HT.NonStat(2)
                            P=Bet(:,iDmn);
                        else
                            P=0;
                        end   
                        G(1:HT.nPrm(2),iDmn) = accumarray(A(:,2),t1(:,iDmn),[HT.nPrm(2),1],@sum,0)+L.*P;
                    end
                    
                    %dL/dsigma
                    if HT.TwoParamFlag
                        G(HT.nPrm(2)+1,:)=[];
                    else
                        t1=(1./S_Dat)-(HT.Delta.*abs(R).^(HT.Delta)./(S_Dat.*(StdTerm).^HT.Delta));
                        for iDmn=1:HT.nDmn-1
                            %Penalty term
                            if HT.NonStat(4)
                                P=S(:,iDmn);
                            else
                                P=0;
                            end   
                            G(HT.nPrm(2)+1:end,iDmn) = accumarray(A(:,4),t1(:,iDmn),[HT.nPrm(4),1],@sum,0)+L.*P;
                        end
                    end
            end
            
            tDel=HT.Delta;
            epsilon=1e-1;%1e-1; %deal with delta==1 issues
            %
            if abs(HT.Delta-1)<epsilon
                tDel=1+epsilon;
            end
            
            
            %% Expected Fisher information matrix
            if nargout==2
                H=NaN(size(G,1),size(G,1),size(G,2));
                switch iP
                    case 1 %alpha and mu
                        %t1=1./V_Dat.*X.^(2-2*Bet_Dat);
                        %t1=(tDel-1).*tDel.*X.^2.*abs(R).^(tDel-2)./((StdTerm).^tDel);
                        %observed Hesian
                        for iDmn=1:HT.nDmn-1 %loop iDmn
                            t1=(tDel-1).*tDel.*X.^2.*abs(R(:,iDmn)).^(tDel-2)./((StdTerm(:,iDmn)).^tDel);
                            %Penalty term
                            if HT.NonStat(1)
                                P=1;
                            else
                                P=0;
                            end   
                            E11 = diag(accumarray(A(:,1),t1,[HT.nPrm(1),1],@sum,0)+L.*P);
                            if  HT.TwoParamFlag
                                H(:,:,iDmn)=E11;
                            else
                                %Penalty term
                                if HT.NonStat(3)
                                    P=1;
                                else
                                    P=0;
                                end   
                                t1=(tDel-1).*tDel.*X.^(2*Bet_Dat(:,iDmn)).*abs(R(:,iDmn)).^(tDel-2)./((StdTerm(:,iDmn)).^tDel);
                                E33 =diag(accumarray(A(:,3),t1,[HT.nPrm(3),1],@sum,0)+L.*P);
                                t1=(tDel-1).*tDel.*X.^(Bet_Dat(:,iDmn)+1).*abs(R(:,iDmn)).^(tDel-2)./((StdTerm(:,iDmn)).^tDel);
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
                    case 2 %beta and sigma
                        %d2L/DLbeta2
                        %t1=(log(X).^2).*(2 + M_Dat.^2./V_Dat);
                        Y_m_aX=(Y-Alp_Dat.*X);
                        for iDmn=1:HT.nDmn-1 %loop iDmn
                            t1=((HT.Delta.*(log(X).^2).*Y_m_aX(:,iDmn))./((StdTerm(:,iDmn)).^HT.Delta)).*...
                                (HT.Delta.*(Y_m_aX(:,iDmn))-Xb(:,iDmn).*M_Dat(:,iDmn)).*abs(R(:,iDmn)).^(HT.Delta-2);
                            %Penalty term
                            if HT.NonStat(2)
                                P=1;
                            else
                                P=0;
                            end   
                            E22 = diag(accumarray(A(:,2),t1,[HT.nPrm(2),1],@sum,0)+L.*P);
                            
                            if  HT.TwoParamFlag
                                H(:,:,iDmn)=E22;
                            else
                                
                                %d2L/DLsigma2
                                %Penalty term
                                if HT.NonStat(4)
                                    P=1;
                                else
                                    P=0;
                                end   
                                t1=(HT.Delta+1).*HT.Delta.*abs(R(:,iDmn)).^(HT.Delta)./(((StdTerm(:,iDmn)).^HT.Delta).*(S_Dat(:,iDmn).^2))-1./(S_Dat(:,iDmn).^2);
                                E44 =diag((accumarray(A(:,4),t1,[HT.nPrm(4),1],@sum,0)+L.*P));
                                
                                %d2L/DLbetasigma
                                t1=zeros(size(X));%HT.Delta.^2.*log(X).*Y_m_aX.*R.*abs(R).^(HT.Delta-2)./((StdTerm).^(HT.Delta).*V_Dat);
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
        
    
        
        %
    end%methods
end
