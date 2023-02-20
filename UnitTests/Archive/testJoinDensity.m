clc; clear; 
n=10000;

S1=[1,0.8;0.8,1];
M1=[1,1];
Z1=mvnrnd(M1,S1,n);

S2=[2,0.95;0.95,0.5];
M2=[0,0];
Z2=mvnrnd(M2,S2,n);

Z=[Z1;Z2];


LockX=4;
nX=30;
X=linspace(0,LockX,nX);
Y=linspace(0,6,nX);

[XG,YG]=meshgrid(X,Y);


Dns=NaN(nX,nX,2);
Dns(:,:,1)=reshape(mvnpdf([XG(:),YG(:)],M1,S1),[nX,nX]);
Dns(:,:,2)=reshape(mvnpdf([XG(:),YG(:)],M2,S2),[nX,nX]);
Dns(:,:,3)=sum(Dns,3); %omni

[FMx,I]=max(Dns(:,end,:));
LockY=Y(squeeze(I));
FMx=squeeze(FMx);


%%
%subplot(1,2,1)
clf;
plot(Z(:,1),Z(:,2),'k.')
axis equal
hold on
for i=1:2
contour(XG,YG,Dns(:,:,i),FMx(i)*[1,1],'linewidth',2,'color','g')
end
contour(XG,YG,Dns(:,:,end),FMx(end)*[1,1],'linewidth',2,'color','r','linestyle','-')
xlim([min(X),max(X)])
ylim([min(Y),max(Y)])
% 
% subplot(2,2,2)
% imagesc(X,Y,Dns(:,:,1))
% axis xy
% axis equal
% 
% subplot(2,2,4)
% imagesc(X,Y,Dns(:,:,2))
% axis xy
% axis equal


