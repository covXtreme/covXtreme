clc; clear;

%% simulate data
n=100000;
x=sort(rand(n,1));
y=sin(x*2)*5+randn(n,1);

%% Model Form
q=0.9; %quantile level 
X=[ones(n,1),x,x.^2];  %%design matrix

%% Fit QR
betahat=ipqr(q, y, X);

%% Compute Prediction
yhat=X*betahat;

%% Plot
clf;
plot(x,y,'k.')
hold on
plot(x,yhat,'r-','linewidth',2);

