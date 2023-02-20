figure(3);
clf;
subplot(2,3,1)
scatter(wave.hs,wave.tp_sea,20,vessel.maxHeave,'filled')
title('Max Heave')
subplot(2,3,2)
scatter(wave.hs,wave.tp_sea,20,vessel.maxPitch,'filled')
title('Max Pitch')
subplot(2,3,3)
scatter(wave.hs,wave.tp_sea,20,vessel.maxRoll,'filled')
title('Max Roll')
subplot(2,3,4)
scatter(wave.hs,wave.tp_sea,20,vessel.maxOffset,'filled')
title('Max OffSet')
subplot(2,3,5)
scatter(wave.hs,wave.tp_sea,20,jacket.bsr,'filled')
title('Base Shear')
subplot(2,3,6)
scatter(wave.hs,wave.tp_sea,20,jacket.otm,'filled')
title('Overturning Moment')

%%
Q=0.999;
clf;
subplot(2,3,1)
I=vessel.maxHeave>quantile(vessel.maxHeave,Q);
hold on
plot(wave.hs,wave.tp_sea,'.','color',[1,1,1]*0.8)
plot(wave.hs(I),wave.tp_sea(I),'k.')
title('Max Heave')
ylim([0,20]); xlim([0,15])
subplot(2,3,2)
hold on
I=vessel.maxPitch>quantile(vessel.maxPitch,Q);
plot(wave.hs,wave.tp_sea,'.','color',[1,1,1]*0.8)
plot(wave.hs(I),wave.tp_sea(I),'k.')
title('Max Pitch')
ylim([0,20]); xlim([0,15])
subplot(2,3,3)
I=vessel.maxRoll>quantile(vessel.maxRoll,Q);
hold on
plot(wave.hs,wave.tp_sea,'.','color',[1,1,1]*0.8)
plot(wave.hs(I),wave.tp_sea(I),'k.')
title('Max Roll')
ylim([0,20]); xlim([0,15])
subplot(2,3,4)
hold on
I=vessel.maxOffset>quantile(vessel.maxOffset,Q);
plot(wave.hs,wave.tp_sea,'.','color',[1,1,1]*0.8)
plot(wave.hs(I),wave.tp_sea(I),'k.')
title('Max OffSet')
ylim([0,20]); xlim([0,15])
subplot(2,3,5)
hold on
I=jacket.bsr>quantile(jacket.bsr,Q);
plot(wave.hs,wave.tp_sea,'.','color',[1,1,1]*0.8)
plot(wave.hs(I),wave.tp_sea(I),'k.')
title('Base Shear')
ylim([0,20]); xlim([0,15])
subplot(2,3,6)
hold on
I=jacket.otm>quantile(jacket.otm,Q);
plot(wave.hs,wave.tp_sea,'.','color',[1,1,1]*0.8)
plot(wave.hs(I),wave.tp_sea(I),'k.')
title('Overturning Moment')
ylim([0,20]); xlim([0,15])

