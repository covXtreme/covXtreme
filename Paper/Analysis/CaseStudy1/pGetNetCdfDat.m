function [X,Tim,LngLtt]=pGetNetCdfDat(tFil,VrbNms);
% function [X,Tim,LngLtt]=pGetNetCdfDat(tFil,VrbNms);
%
% Get hindcast data (X) from a netcdf file
% Always get the dateser (Tim), longitude and latitude (LngLtt)

if nargin==0
    tFil='NORA10_6176N_1434W.nc';
    VrbNms={'Hs';'Tp'};
end;

nVrb=size(VrbNms,1);

%% Always get the time stamp, called dateser
tStr={ncinfo(tFil).Variables.Name}';
tIdnTim=find(strcmp(tStr,'dateser')==1);
if isempty(tIdnTim)==1;
        fprintf(1,'Error: dateser not found. Terminating\n');
else;
    Tim=double(ncread(tFil,'dateser'));
    nX=size(Tim,1);
    X=nan(nX,nVrb);
end;

%% Always get the longitude and latitude, called dateser
tStr={ncinfo(tFil).Variables.Name}';
tIdnTim=find(strcmp(tStr,'longitude')==1);
if isempty(tIdnTim)==1;
        fprintf(1,'Error: longitude not found. Terminating\n');
else;
    LngLtt(1,1)=double(ncread(tFil,'longitude'));
end;
tIdnTim=find(strcmp(tStr,'latitude')==1);
if isempty(tIdnTim)==1;
        fprintf(1,'Error: latitude not found. Terminating\n');
else;
    LngLtt(2,1)=double(ncread(tFil,'latitude'));
end;

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
    X(:,iV)=double(ncread(tFil,VrbNms{iV}))*tScl;
    
end;

return;