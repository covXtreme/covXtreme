clc; clear;

n=3000;

rho=0.4;
mu=[0,0];
Sgm=[1,rho;rho,1];

Z=mvnrnd(mu,Sgm,n);

f=mvnpdf(Z,mu,Sgm);

U=normcdf(Z);

f_u=1/sqrt(det(Sgm))*exp(-0.5.*diag(Z*(inv(Sgm)-eye(2))*Z'));

f_u2=abs(1./prod(normpdf(Z),2)).*f;


clf;
subplot(2,2,1)
scatter(Z(:,1),Z(:,2),30,f,'filled');
title('Original Margins')

subplot(2,2,2)
scatter(U(:,1),U(:,2),30,log10(f_u),'filled');
title('Uniform Margins')
colorbar

subplot(2,2,3)
scatter(U(:,1),U(:,2),30,log10(f_u2),'filled');
title('Uniform Margins')
colorbar

subplot(2,2,4)
scatter(U(:,1),U(:,2),30,f_u./f_u2,'filled');
title('Uniform Margins')
colorbar