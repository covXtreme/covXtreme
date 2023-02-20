classdef UT_MarginalModel< matlab.unittest.TestCase
    %UT_MarginalModel
    % Unit testing for the marginal model 
    
    properties
        Mrg
        Dat
    end
    
    properties (ClassSetupParameter)
        %%simple case no parameter uncertainty single return period
        Sgm={1,5};
        Xi={-0.1,0};
        Omg={1,10};
        Kpp={0.1,10};
        GmmLct={0,5};
        Tau={0.5,0.9};
        Rat={1,10};
        nBin={10};
        nBoot={5};
    end%ClassSetupParameter
    
    properties(TestParameter)
        X={linspace(0,20,100)'};
        P={linspace(0,1,100)'};
    end %TestParameter
    
    methods (TestClassSetup)
        function ClassSetup(obj,Tau,Omg,Kpp,GmmLct,Xi,Sgm,Rat,nBin,nBoot)
            obj.Mrg=MarginalModel;
            
            if nBin>1 || nBoot>1 %add "bootstrap" uncertainty
                Xi=Xi+zeros(nBoot,1);
                Sgm=exp(log(Sgm)+randn(nBin,nBoot)*0.1);
                Omg=exp(log(Omg)+randn(nBin,nBoot)*0.1);
                Kpp=exp(log(Kpp)+randn(nBin,nBoot)*0.1);
                Rat=exp(log(Rat)+randn(nBin,nBoot)*0.1);
                GmmLct=exp(log(GmmLct)+randn(nBin,1)*0.1);
                Tau=Tau+(rand(nBoot,1)-0.5)*0.1;
                Thr=MarginalModel.gaminv(Tau',Omg,Kpp,GmmLct);
            else
                Thr=MarginalModel.gaminv(Tau',Omg,Kpp,GmmLct);
            end
            
            obj.Mrg.Bn=CovariateBinning;
            obj.Mrg.Bn.nBin=nBin;
            obj.Mrg.nBoot=nBoot;
            obj.Mrg.Shp=Xi;
            obj.Mrg.Scl=Sgm;
            obj.Mrg.Omg=Omg;
            obj.Mrg.Kpp=Kpp;
            obj.Mrg.Thr=Thr;
            obj.Mrg.GmmLct=GmmLct;
            obj.Mrg.NEP=Tau;
            obj.Mrg.Rat=Rat;
        end
    end
    
    methods ( Test )
        function UT_CDF(obj,X)
            Pr=obj.Mrg.CDF(permute(X,[3,2,1]));
            %Check range
            obj.assertGreaterThanOrEqual(Pr,0);
            obj.assertLessThanOrEqual(Pr,1);
            obj.assertSize(Pr,[obj.Mrg.Bn.nBin,obj.Mrg.nBoot,numel(X)]);
            %Check monotonic
            dP=diff(Pr,1,3);
            obj.assertGreaterThanOrEqual(dP,0);
        end
        
        function UT_PDF(obj,X)
            f=obj.Mrg.PDF(permute(X,[3,2,1]));
            %Check positive
            obj.assertGreaterThanOrEqual(f,0);
            obj.assertSize(f,[obj.Mrg.Bn.nBin,obj.Mrg.nBoot,numel(X)]);                        
            
        end
        
        function UT_PIT(obj)
            %check proability integral transform go from orignal to standard margins and back
            
            nRls=1000;
            
            I=randi(obj.Mrg.nBoot,nRls,1);%decide which bootstrap sample to use over all bootstraps HT.nBoot
            A=randi(obj.Mrg.Bn.nBin,nRls,1);  %bin allocation
            
            U= rand(nRls,1); %sample uniform value
            G=obj.Mrg.INV_Standard(U); %transform to standard margins;
            
            U2=obj.Mrg.CDF_Standard(G); %transform back to uniform;
                        
            obj.assertEqual(U,U2,'AbsTol',1e-2);
            
            %transform from uniform to gamgp                       
            Z=obj.Mrg.INV(U,I,A);
            
            %transform from gamgp to uniform                       
            U3=obj.Mrg.CDF(Z,A,I);
            
            obj.assertEqual(U,U3,'AbsTol',1e-2);
            
        end
        
        %PIT
        
        %Fit
        
    end
end

