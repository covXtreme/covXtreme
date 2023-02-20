%tFil='NORA10_6176N_1434W.nc';
%tFil='NORA10_6176N_1594W.nc';
%tFil='NORA10_6176N_1708W.nc'; 
%tFil='NORA10_6600N_1079W.nc';
%tFil='NORA10_5500N_2069W.nc';
%tFil='NORA10_7000N_0154W.nc';

%%
tFil='NORA10_5277N_0346E.nc'; %clear change in slope of Hs-T2, not much swell!
VrbNms={'Hs';'T2';'PeakWaveDirection'}; 

%% 
ncinfo(tFil).Variables.Name

%%
[X,Tim,LngLtt]=pGetNetCdfDat(tFil,VrbNms);
save Dat X Tim LngLtt;

figure(1); clf; hold on;
load coastlines;
worldmap([45 75],[-30 15])
plotm(coastlat, coastlon);
plotm(LngLtt(2),LngLtt(1),'r*')

figure(2); clf;
subplot(2,3,1); plot(Tim,X(:,1),'ko'); datetick;
subplot(2,3,4); plot(Tim,X(:,2),'ko'); datetick;
subplot(2,3,2); plot(X(:,3),X(:,1),'ko'); 
subplot(2,3,5); plot(X(:,3),X(:,2),'ko'); 
subplot(1,3,3); plot(X(:,1),X(:,2),'ko');

figure(3); clf;
pPltDrcSct(X(:,3),X(:,1:2));

%% Peak pick



%% PPC
