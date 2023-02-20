%Combine data from individual locations into one structure for analysis
clear all;
DatFld='R:\Public Folder\Extremes\Development\Students\Shiraz\Data';
Nms={'Anga', 'Brent', 'Clipper', 'Sean'};
nLct=numel(Nms);
Drc = cell(nLct,1);

S=cell(nLct,1);
for iL=1:nLct;
    tStr=sprintf('S{%g}=load(fullfile(DatFld,''%s'',''D_%s.mat''));',iL,Nms{iL},Nms{iL});
    eval(tStr);
end;

clear D;

Inc=90; %sector size
DrcBnd=[(0:Inc:360-Inc)' (0:Inc:360-Inc)'+Inc];
ns=size(DrcBnd,1)+1; %number of sectors
DatAll = cell(nLct,1);
for iL=1:4  %loop over locations ANGA, BRENT, CLIP, SEAN   
    Hs=cell(ns,1);
    SrgMxm=cell(ns,1);
    SrgRng=cell(ns,1);
    SrgNgtMnm=cell(ns,1);
    SrgMdn=cell(ns,1);
    t= ones(S{iL}.D.Mxm.n,1);
    t=find(t==1); %locate members of subset
    
   
    ymax = NaN(length(t),1);
    for i=1:length(t)
        ymax(i) = max(S{iL}.D.Mxm.SrgTrj{t(i)},[],'omitnan');
    end
    IMax = ~isnan(ymax);

    ymin = NaN(length(t),1);
    for i=1:length(t)
        ymin(i) = - min(S{iL}.D.Mxm.SrgTrj{t(i)},[],'omitnan');
    end
    IMin = ~isnan(ymin);

    
    yrange = NaN(length(t),1);
    for i=1:length(t)
        yrange(i) = range(S{iL}.D.Mxm.SrgTrj{t(i)},'omitnan');
    end
    IRng = ~isnan(yrange);
    
    ymedian = NaN(length(t),1);
    for i=1:length(t)
        ymedian(i) = median(S{iL}.D.Mxm.SrgTrj{t(i)},'omitnan');
        if isnan(ymedian(i))
            warning('.')
            sprintf('Nan in median calc for Lct %d, value %d',iL,i)
            %Note: Emma added this check - shiraz was carrying NaNs around
            %which SrgTrj data for given iL and t(i) only included NaNs
        end
    end
    IMdn = ~isnan(ymedian);
    
    tHs = S{iL}.D.Mxm.Wav;
    IHs = ~isnan(tHs);
    Yrs = range(S{iL}.D.Data.Tim)/365.25;
    tDrc=S{iL}.D.Mxm.WavDrc;
    IDrc = ~isnan(tDrc);
    
    
       
    DatAll{iL}.Hs=tHs(IDrc);  %don't take Nans out yet, will be diff for each Analysis type (surge etc)
    DatAll{iL}.SrgMxm=ymax(IDrc);
    DatAll{iL}.SrgNgtMnm=ymin(IDrc);
    DatAll{iL}.SrgRng=yrange(IDrc);
    DatAll{iL}.SrgMdn=ymedian(IDrc);
    DatAll{iL}.Nam=Nms{iL};
    DatAll{iL}.Yrs=Yrs;
    DatAll{iL}.WavDrc = tDrc(IDrc);
    


end;
if 0
    %% Check sector
    tI = DatAll{1}.WavDrc>=0 & DatAll{1}.WavDrc<=90;
    tI=find(tI);
    figure(1); clf;
    plot(DatAll{1}.Hs(tI),DatAll{1}.SrgMdn(tI),'.')
    grid on
    xlim([2,9])
    ylim([-0.6,0.6])
    
    
    %% -----------
    figure(1); clf;
    plot(D{1}.Hs{5},D{1}.SrgMdn{5},'k.')
    grid on
    xlim([2,9])
    ylim([-0.6,0.6])
    
    tIn=find(D{1}.Hs{3}==2.7)
    Hs27=[D{1}.Hs{3}(tIn),D{1}.SrgMdn{3}(tIn)];
end


%% Put into PPC format

% Dat struct with fields Y:[nObs x 2], X:[nOBs x 1], Lbl = {'.' '.'}
AnlNms={'NgtMnm' 'Mdn' 'Mxm' 'Rng'}; 
nAnl = numel(AnlNms); %4 analyses: min, max, median, rng

Fld = 'R:\Emma.Ross\Extremes\CE_Surge\Analysis';
for iL = 1:nLct  %Loop over locations
    FldLct= fullfile(Fld,Nms{iL});
    if ~exist(FldLct,'dir')
        mkdir(FldLct)
    end
    for iA = 1:nAnl  %loop over analysis types / surge summary
        FldLctAnl=fullfile(FldLct,AnlNms{iA});
        if ~exist(FldLctAnl,'dir')
             mkdir(FldLctAnl)
        end
        Dat.Lbl = {'Hs' ['Srg',AnlNms{iA}]};
        tStr=sprintf('DatAll{iL}.Srg%s',AnlNms{iA});
        tY2=eval(tStr);
        
        IY2=~isnan(tY2);
        IHs=~isnan(DatAll{iL}.Hs);
        IDrc=~isnan(DatAll{iL}.WavDrc);
        IDat=logical(IHs.*IDrc.*IY2);
        
        Dat.Y= [DatAll{iL}.Hs(IDat), tY2(IDat)];  %some nans in Y
        Dat.X = DatAll{iL}.WavDrc(IDat);
        Dat.Yrs = DatAll{iL}.Yrs;
        
        save(fullfile(FldLctAnl,'Dat'),'Dat')

    end
end



%% Plot to compare against Shiraz's data 
load('R:/Emma.Ross/Extremes/CE_Surge/Analysis/Brent/Mxm/Dat.mat');
tI = Dat.X>=0 & Dat.X<=90;
tI=find(tI);
figure(3); clf;
plot(Dat.Y(tI,1),Dat.Y(tI,2),'.')

