function [Rsp,Cvr,Asc] = PreProcessData(RspOrg,CvrOrg,AscOrg, VslHd)
   %Function to pre-process case study data: remove "0" responses due to 
   %...heading measurement error giving 90deg
   
   KpInd = ~any(VslHd==90 | AscOrg== 0,2);
   Asc = AscOrg(KpInd,:); %main response
   Rsp = RspOrg(KpInd); %covariate (direction)
   Cvr = CvrOrg(KpInd,:); %associated variable(s)

end
