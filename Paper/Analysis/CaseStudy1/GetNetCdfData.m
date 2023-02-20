function pGetNetCdfDat(tFil,VrbNms);

if nargin==0
    tFil='NORA10_6176N_1434W.nc';
    VrbNms={'Hs';'Tp'};
end;

nVrb=size(VrbNms,1);

%% Always get the time stamp, called dateser
tX=double(ncread(tFil,'dateser'));
nX=size(tX,1);

X=nan(nX,nVrb+1);
X(:,1)=tX; %populate with time stamp
for iV=1:nVrb;
    
    fprintf(1,'File %s, variable %s\n',tFil,VrbNms{iV});
    
    %% Check input
    tStr={ncinfo(tFil).Variables.Name}';
    tIdn1=find(strcmp(tStr,VrbNms{iV})==1);
    if isempty(tIdn1)==1;
        fprintf(1,'Error: variable not found. Terminating\n');
        return;
    end;
    tStr={ncinfo(tFil).Variables(tIdn1).Attributes.Name}';
    tIdn2=find(strcmp(tStr,'scale')==1);
    if isempty(tIdn2)==1;
        fprintf(1,'Error: scale factor not found. Terminating\n');
        return;
    end;
    tScl=ncinfo(tFil).Variables(tIdn1).Attributes(tIdn2).Value;
    X(:,iV+1)=double(ncread(tFil,VrbNms{iV}))*tScl;
    
end;

return;