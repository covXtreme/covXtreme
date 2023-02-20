        function MM=SimulateMarginal(MM,nDat,DrcEdg,parTrue,Thr)
            %Simulate marginal GP data: piecewise sigma, constant xi
            validateattributes(nDat, {'numeric','scalar'},{},'SimulateMarginal','nDat',2);
            validateattributes(DrcEdg, {'numeric'},{'ncols', 1},'SimulateMarginal','DrcEdg',3);
            MM.nBin=length(DrcEdg);
            validateattributes(parTrue, {'numeric'},{'size', [MM.nBin+1,1]},'MarginalModel','parTrue',4);
            validateattributes(Thr, {'numeric'},{'size', [MM.nBin,1]},'MarginalModel','Thr',5);
            MM.nDat=nDat;
            MM.DrcEdg=DrcEdg;
            %simulate
            MM.DrcX=sort(rand(MM.nDat,1)*360);
            MM.XiTrue=parTrue(1,1);
            MM.SigTrue=parTrue(2:(MM.nBin+1),1);
            MM.BinAlc=discretize(MM.DrcX,[0;MM.DrcEdg;360]);
            MM.BinAlc(MM.BinAlc==MM.nBin+1)=1;
            MM.Y=gprnd(MM.XiTrue,MM.SigTrue(MM.BinAlc),Thr(MM.BinAlc)); %constant xi, pieceise sig, thr
            
        end %SimulateMarginal