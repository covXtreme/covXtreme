clc; clear;
%non uniform g

n=10000;

%% Model paramters
a=1;
b=0.2;
m=0;
s=0.3;

xi=[-0.3,-0.1];
sig=[2,2];
alp=[2,2];
bet=[5,5];
v=[0,0];
tau=[0.5,0.5];

thr=MarginalModel.gaminv(tau,alp,bet,v);

HTNEP=0.5;

nRsd=100;
Rsd=sort(trnd(5,nRsd,1));

IS.EndPoint=thr-sig./xi;

%% Monte Carlo direct
MC.U=sort(rand(n,1))*(1-HTNEP)+HTNEP;  %Unif [NEP,1]
MC.G=-log(-log(MC.U)); %gumbel margins
%MC.R=randn(n,1); %generate residuals
MC.R=Rsd(randi(nRsd,n,1)); %generate residuals
MC.G=[MC.G,a*MC.G+MC.G.^b.*(m+s.*MC.R)];
MC.U=[MC.U,exp(-exp(-MC.G(:,2)))];
MC.Z=MarginalModel.gamgpinv(MC.U,xi,sig,thr,alp,bet,v,tau);

%% IS version
UL=MarginalModel.gamgpinv(1-1e-9,xi(1),sig(1),thr(1),alp(1),bet(1),v(1),tau(1));
%original scale
IS.Z=rand(n,1).*(UL(1)-thr(1))+thr(1);
IS.logg_Z=-log(UL(1)-thr(1)).*ones(n,1); %deinsty of sampled points on GamGP scale
IS.logf_ZMrg=MarginalModel.loggamgppdf(IS.Z,xi(1),sig(1),thr(1),alp(1),bet(1),v(1),tau(1));
%uniform
IS.U=MarginalModel.gamgpcdf(IS.Z,xi(1),sig(1),thr(1),alp(1),bet(1),v(1),tau(1));
IS.logg_U=IS.logg_Z-IS.logf_ZMrg;

IS.G=-log(-log(IS.U));


IS.BndRsd=[min(Rsd),max(Rsd)];
%IS.
tmu=a*IS.G(:,1)+m.*IS.G(:,1).^b;
tsig=s.*IS.G(:,1).^b;
IS.R=rand(n,1).*(IS.BndRsd(2)-IS.BndRsd(1))+IS.BndRsd(1);
IS.G(:,2)=IS.R.*tsig+tmu;
IS.logf_gMrg=-(IS.G +exp(-IS.G)); %log Gumbel density
IS.logg_g(:,1)=IS.logg_U+IS.logf_gMrg(:,1);
IS.logg_g(:,2)=-log((IS.BndRsd(2)-IS.BndRsd(1)).*tsig);
%transfrom back
IS.U(:,2)=exp(-exp(-IS.G(:,2)));
IS.logg_U(:,2)=IS.logg_g(:,2)-IS.logf_gMrg(:,2);


IS.Z(:,2)=MarginalModel.gamgpinv(IS.U(:,2),xi(2),sig(2),thr(2),alp(2),bet(2),v(2),tau(2));
IS.logf_ZMrg(:,2)=MarginalModel.loggamgppdf(IS.Z(:,2),xi(2),sig(2),thr(2),alp(2),bet(2),v(2),tau(2));
IS.logg_Z(:,2)=IS.logg_U(:,2)+IS.logf_ZMrg(:,2);
UL(2)=max(IS.Z(:,2));
% gumbel scale density
%IS.logf_g2cndg1=-log(2*pi)/2-log(tsig)-0.5*IS.R.^2; %conditional density Gumbel
IS.logf_g2cndg1=log(ksdensity(Rsd,IS.R)); %conditional density Gumbel
IS.logf_g2cndg1(isnan(IS.logf_g2cndg1))=-Inf;
IS.logf_g=IS.logf_g2cndg1+IS.logf_gMrg(:,1); %joint density Gumbel
%Uniform scale density
IS.logJ_Stn_u=-sum(IS.logf_gMrg,2); %Gumbel --> Uniform transform Jacobian
IS.logf_u=IS.logJ_Stn_u+IS.logf_g; %joint density uniform
IS.logf_u(isnan(IS.logf_u))=-Inf;
% Gam GP scale denisty
IS.logJ_U_GamGP=sum(IS.logf_ZMrg,2); %Uniform --> GamGP transform Jacobian
IS.logf_Z=IS.logJ_U_GamGP + IS.logf_u;

%% Compute regression lines on original scale

nG=40;
Edg=linspace(thr(1),UL(1),nG+1);
Edg2=linspace(v(2),UL(2),nG+1);
EdgStn=linspace(0,10,nG+1);
Md=(Edg(2:end)+Edg(1:end-1))/2;
Md2=(Edg2(2:end)+Edg2(1:end-1))/2;
MdStn=(EdgStn(2:end)+EdgStn(1:end-1))/2;
A=discretize(IS.Z(:,1),Edg);
AStn=discretize(IS.G(:,1),EdgStn);
[Q,Q2,Q3,QStn]=deal(NaN(nG,3));

[tY2,srtI]=sort(IS.Z(:,2));
Dst=(Md-IS.Z(srtI,1))./(Edg(end)-Edg(1));
nu=0.01;
tW=exp(IS.logf_Z(srtI)-Dst.^2./(nu.^2));    
tP3=cumsum(tW,1)./sum(tW,1);

for iA=1:nG

    %Standard Margins CAse
    
    tI=AStn==iA;
    tY=IS.G(tI,2);
    tW=exp(IS.logf_g2cndg1(tI)-IS.logg_g(tI));
    
    [tY,srtI]=sort(tY);
    tW=tW(srtI);
    tP=cumsum(tW)./sum(tW);
    
    if any(tP)
        [sPStn,JStn]=unique(tP);
        if numel(sPStn)>1
            
            QStn(iA,:)=interp1(sPStn,tY(JStn),[0.025,0.5,0.975]);
        end
    end
    
    %% Other cases
    tI=A==iA;
    tY=IS.Z(tI,2);
    tW=exp(IS.logf_Z(tI));    
    tW2=exp(IS.logf_g2cndg1(tI)-IS.logg_g(tI));    
    
    [tY,srtI]=sort(tY);
    tW=tW(srtI);
    tW2=tW2(srtI);
            
    tP=cumsum(tW)./sum(tW);
    tP2=cumsum(tW2)./sum(tW2);
    
    if any(tP)
        [sP,J]=unique(tP);
        [sP2,J2]=unique(tP2);
        [sP3,J3]=unique(tP3(:,iA));
      %  [sPStn,JStn]=unique(tPStn(:,iA));
        if numel(sP)>1
            Q(iA,:)=interp1(sP,tY(J),[0.025,0.5,0.975]);
            Q2(iA,:)=interp1(sP2,tY(J2),[0.025,0.5,0.975]);
            Q3(iA,:)=interp1(sP3,tY2(J3),[0.025,0.5,0.975]);
           % QStn(iA,:)=interp1(sPStn,tYStn(JStn),[0.025,0.5,0.975]);
        end
    end
end

%%
xlmTGmb=-log(-log([0.5,1-1e-7]));
RngPrd=linspace(xlmTGmb(1),xlmTGmb(2),100);
%Figure
figure(1);
clf;
mu=a*RngPrd+m.*RngPrd.^b;
sig=s.*RngPrd.^b;
axG(1)=subplot(4,2,1);
scatter(IS.G(:,1),IS.G(:,2),10,exp(IS.logf_g),'filled');
hold on
plot(RngPrd,mu,'g-','linewidth',2);
plot(RngPrd,mu+1.96*sig,'g--','linewidth',2);
plot(RngPrd,mu-1.96*sig,'g--','linewidth',2);
colorbar
title('IS Density Gumbel Margins')
axG(2)=subplot(4,2,2);
plot(MC.G(:,1),MC.G(:,2),'k.')
hold on

plot(RngPrd,mu,'g-','linewidth',2);
plot(RngPrd,mu+1.96*sig,'g--','linewidth',2);
plot(RngPrd,mu-1.96*sig,'g--','linewidth',2);
hold on
plot(MdStn,QStn(:,2),'r-','linewidth',1.5);
plot(MdStn,QStn(:,[1,3]),'r--','linewidth',1);

title('MC simulation Gumbel Margins')
linkaxes(axG,'xy');
xlim(xlmTGmb)

subplot(4,2,3)
scatter(IS.U(:,1),IS.U(:,2),10,exp(IS.logf_u),'filled');
xlim([0,1]);
ylim([0,1]);
caxis([0,2])
colorbar
title('IS: Density  Uniform Margins')
subplot(4,2,4)
plot(MC.U(:,1),MC.U(:,2),'k.')
xlim([0,1]);
ylim([0,1]);
title('MC Simulation  Uniform Margins')
colorbar
axZ(1)=subplot(4,2,5);
scatter(IS.Z(:,1),IS.Z(:,2),10,exp(IS.logf_Z),'filled');
title('IS: Density GamGP Margins')
colorbar
caxis(quantile(exp(IS.logf_Z),[0,0.99]))
axZ(2)=subplot(4,2,6);
plot(MC.Z(:,1),MC.Z(:,2),'k.');
hold on
plot(Md,Q(:,2),'r-','linewidth',1.5);
plot(Md,Q(:,[1,3]),'r--','linewidth',1);
plot(Md,Q2(:,2),'g-','linewidth',1.5);
plot(Md,Q2(:,[1,3]),'g--','linewidth',1);
plot(Md,Q3(:,2),'b-','linewidth',1.5);
plot(Md,Q3(:,[1,3]),'b--','linewidth',1);
title('MC Simulation GamGP Margins')
linkaxes(axZ,'xy');
%colormap(jet)
subplot(4,2,7);
scatter(IS.Z(:,1),IS.Z(:,2),10,exp(IS.logf_g),'filled');
title('Uncorrected Density GamGP Margins')
colorbar

[Xi,Yi]=meshgrid(Md,Md2);
xi=[Xi(:),Yi(:)];

PtVlm=1e-2.*(UL-[thr(1),v(2)]);
%BW=0.2;
f = ksdensity(IS.Z,xi,  'Weights' ,exp(IS.logf_Z), 'Bandwidth',PtVlm);
subplot(4,2,8)

Ind=[8,15,25];
Lck=[Md(Ind)',Q(Ind,2)];
plot(MC.Z(:,1),MC.Z(:,2),'k.');
[~,LckInd]=min(sum((permute(Lck,[3,2,1])-xi).^2,2));
hold on
contour(unique(xi(:,1)),unique(xi(:,2)),reshape(f,[nG,nG]),f(LckInd),'linecolor','r','linewidth',1.5);

title('Density Contour')

%%

figure(2);
clf;
scatter(IS.G(:,1),IS.G(:,2),40,exp(IS.logf_g2cndg1),'filled');
hold on

plot(RngPrd,mu,'g-','linewidth',2);
plot(RngPrd,mu+1.96*sig,'g--','linewidth',2);
plot(RngPrd,mu-1.96*sig,'g--','linewidth',2);
%xlim(xlmTGmb)
colormap('jet')
%ylim(xlmTGmb)
hold on
plot(MdStn,QStn(:,2),'m-','linewidth',2);
plot(MdStn,QStn(:,[1,3]),'m--','linewidth',2);

%%
figure(3);
clf;
PtVlm=1e-2.*(UL-[thr(1),v(2)]);
%BW=0.2;
f = ksdensity(IS.Z,xi,  'Weights' ,exp(IS.logf_Z), 'Bandwidth',PtVlm);
%f = ksdensity(IS.Z,xi,  'Weights' ,exp(IS.logf_Z), 'Bandwidth',[0.1,0.5]);
Ind=[8,15,25];
Lck=[Md(Ind)',Q(Ind,2)];
plot(MC.Z(:,1),MC.Z(:,2),'k.');
[~,LckInd]=min(sum((permute(Lck,[3,2,1])-xi).^2,2));
hold on
contour(unique(xi(:,1)),unique(xi(:,2)),reshape(f,[nG,nG]),f(LckInd),'linecolor','r','linewidth',1.5);
hold on
plot(Lck(:,1),Lck(:,2),'g.','markersize',20)
title('Density Contour')