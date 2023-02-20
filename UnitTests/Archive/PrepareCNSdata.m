%% CNS Data Preparation script

clear variables
close all

disp('Select an option:')
disp('1. Create a mat file for CNS for metocean and vessel response.')
iOption = input('Selection: ');


%% Create a mat file for CNS for metocean and vessel response
if iOption==1
    
    %Load file
    tso = tsoget('mtbfile','F:\10 - Ad-Hoc and Project Work\Central North Sea\Curlew\2016-06 Curlew GER - development of SI\05. Analysis\tsoStormCMB_Job28.dlp');
    
    %Extract parameters of interest
    [t, wave.hs]     = tsogetpar(tso,'Hs','-inclDate');
    wave.dm          = tsogetpar(tso,'MeanWaveDirection');
    wind.spd         = tsogetpar(tso,'WindSpeed');
    wind.dirn        = tsogetpar(tso,'WindDirection');
    curr.spd         = tsogetpar(tso,'CurrentSpeed');
    curr.dirn        = tsogetpar(tso,'CurrentDirection');
    wave.hs_sea      = tsogetpar(tso,'Hs_Sea');
    wave.tp_sea      = tsogetpar(tso,'Tp_Sea');
    wave.tz_sea      = tsogetpar(tso,'Tz_Sea');
    wave.dm_sea      = tsogetpar(tso,'Dm_Sea');
    wave.hs_swell    = tsogetpar(tso,'Hs_Swell1');
    wave.tp_swell    = tsogetpar(tso,'Tp_Swell1');
    wave.tz_swell    = tsogetpar(tso,'Tz_Swell1');
    wave.dm_swell    = tsogetpar(tso,'Dm_Swell1');
    vessel.maxOffset   = tsogetpar(tso,'offset_omni');
    vessel.heading     = tsogetpar(tso,'vesselHeading');
    vessel.maxHeave    = tsogetpar(tso,'heave');
    vessel.maxPitch    = tsogetpar(tso,'pitch');
    vessel.maxRoll     = tsogetpar(tso,'roll');
    
    %Calculate BS and OTM for a certain water depth. Look at Stokes Vth in
    %100m of water
    BSA = [57.16152	1.81967	8.63239	1.20383	14.71889	0.01269];
    OMA = [3729.51500	99.54618	987.45150	114.33710	1637.85400	0.59597];

    %Estimate certain parameters for the loading equation just to get some
    %reasonable values
    a = 1.71*wave.hs/2;  % wave amplitude
    phi = 0.96;     % spreading estimate
    u = curr.spd;         % depth-average current speed
    cos_th = cosd(wave.dm-curr.dirn);  % cos of angle between currents and waves
    T = wave.tz_sea;    % mean wave period estimate
    
    %Calculate each load term for BS and OTM
    BS1 = BSA(1)*u.^2;
    BS2 = BSA(2)*u.*a.*T.*phi.*cos_th;
    BS3 = BSA(3)*u.*phi.*a.^2.*cos_th./T;
    BS4 = BSA(4)*(phi*a).^2;
    BS5 = BSA(5)*phi^2*a.^3./T.^2;
    BS6 = BSA(6)*(phi*a.*T).^2;
    jacket.bsr = BS1+BS2+BS3+BS4+BS5+BS6;
    
    OM1 = OMA(1)*u.^2;
    OM2 = OMA(2)*u.*a.*T.*phi.*cos_th;
    OM3 = OMA(3)*u.*phi.*a.^2.*cos_th./T;
    OM4 = OMA(4)*(phi*a).^2;
    OM5 = OMA(5)*phi^2*a.^3./T.^2;
    OM6 = OMA(6)*(phi*a.*T).^2;
    jacket.otm = OM1+OM2+OM3+OM4+OM5+OM6;
    
    %Only take data where currents exist
    igood = find(isfinite(curr.spd));
    t                =      t(igood);
    wave.hs          =	    wave.hs(igood);          
    wave.dm          =	    wave.dm(igood);          
    wind.spd         =	    wind.spd(igood);         
    wind.dirn        =	    wind.dirn(igood);        
    curr.spd         =	    curr.spd(igood);         
    curr.dirn        =	    curr.dirn(igood);        
    wave.hs_sea      =	    wave.hs_sea(igood);      
    wave.tp_sea      =	    wave.tp_sea(igood);      
    wave.tz_sea      =	    wave.tz_sea(igood);      
    wave.dm_sea      =	    wave.dm_sea(igood);      
    wave.hs_swell    =	    wave.hs_swell(igood);    
    wave.tp_swell    =	    wave.tp_swell(igood);    
    wave.tz_swell    =	    wave.tz_swell(igood);    
    wave.dm_swell    =      wave.dm_swell(igood);    
    vessel.maxOffset =	    vessel.maxOffset(igood);   
    vessel.heading   =	    vessel.heading(igood);     
    vessel.maxHeave  =	    vessel.maxHeave(igood);    
    vessel.maxPitch  =	    vessel.maxPitch(igood);    
    vessel.maxRoll   =	    vessel.maxRoll(igood);     
    jacket.bsr       =      jacket.bsr(igood);
    jacket.otm       =      jacket.otm(igood);
    
    
    %save useful info
    save CNS_mo_response t wave wind curr vessel jacket
    
    
end
