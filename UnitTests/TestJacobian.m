clc; clear;

n=3000;

xi=-0.1;
sig=3;
alp=1;
bet=1;

v=0;
tau=0.7;
thr=gaminv(tau,alp,bet);

MrgCs=2; %case 1 Gmb, case 2 Lap

U=sort(rand(n,1));
switch MrgCs
    case 1
        G=-log(-log(U));
    case 2
        G = sign(0.5-U).*log(1-2*abs(U-0.5));         
end

P=MarginalModel.gamgpinv(U,xi,sig,thr,alp,bet,v,tau);


switch MrgCs
    case 1 
        f_Stn=exp(-(G +exp(-G))); %Gumbel density
        J_Stn_u=(-1./(U.*log(U))); %Gumbel --> Uniform transform Jacobian
    case 2
        f_Stn =  0.5.*exp(-abs(G)); %Laplace density
        J_Stn_u=sign(0.5-U)./min(U,(1-U)); %Laplace --> Uniform transform Jacobian
end

f_u=abs(J_Stn_u).*f_Stn;

V=(U-tau)./(1-tau);
J_u_gp=sig.*((U-1)./(tau-1)).^(-xi)./(1-U); %Uniform --> GP transform Jacobian




f_gp=abs(1./J_u_gp).*f_u;
f_gp2=MarginalModel.gamgppdf(P,xi,sig,thr,alp,bet,v,tau);



clf;
subplot(3,1,1)
plot(G,f_Stn,'k-','linewidth',2);
hold on
histogram(G,'normalization','pdf');

hold on

subplot(3,1,2)
plot(U,f_u,'k-','linewidth',2)
hold on
plot(U,1,'g-','linewidth',2);

histogram(U,'normalization','pdf');

subplot(3,1,3);
hold on
plot(P,f_gp,'k-','linewidth',2)
hold on
plot(P,f_gp2,'g--','linewidth',2);

histogram(P,100,'normalization','pdf');
