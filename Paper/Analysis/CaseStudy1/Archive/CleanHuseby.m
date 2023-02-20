function XNew=CleanHuseby(X);
% function X=CleanHuseby(X);
%
% Clean up Huseby contours around "cusps"
%
% Algorithm works by searching for and eliminating points that cannot be
% accessed from the contour median without crossing a line segment (between
% adjacent points on the contour)
%
% Philip Jonathan 20221103

% Verbose mode
Vrb=0;

% Test mode
if nargin==0;
    Om=(0:0.1:2*pi)';
    nOm=size(Om,1);
    X=[cos(Om)+0.1*randn(nOm,1) sin(Om)+0.1*randn(nOm,1)];
end;

% Drop missing values
X=X(sum(~isnan(X),2)==2,:);
nX=size(X,1);

% Kep is an indicator to good points (to KEeP)
Kep=ones(nX,1);
nKep=sum(Kep);
for iX=1:nX;
    
    % Drop a single point
    Kep=CleanHusebyDropOne(X,Kep);
    if sum(Kep)==nKep-1;
        nKep=nKep-1;
    else;
        % Stop if no points for dropping
        break;
    end;
    
end;

% Return the cleaned contour, padded with NaN if necessary
XNew=nan(nX,2);
XNew(1:sum(Kep==1),:)=X(Kep==1,:);

return;

function Kep=CleanHusebyDropOne(X,Kep);

% Set verbose mode
Vrb=0; 

% Identifier
Idn=(1:size(X,1))';

% Keep only good points
X=X(Kep==1,:);
Idn=Idn(Kep==1);
nX=size(X,1);

% Find the step vector between points on the approximate contour
Dlt=nan(nX,2);
for iX=1:nX;
    if iX<nX;
        Dlt(iX,:)=[X(iX+1,1)-X(iX,1) X(iX+1,2)-X(iX,2)];
    else;
        Dlt(nX,:)=[X(1,1)-X(nX,1) X(1,2)-X(nX,2)];
    end;
end;

% Define an origin
Org=median(X)';

for iX=1:nX;
     
    % Candidate point to drop
    C=X(iX,:)'-Org;
    
    if Vrb==1;
        clf; hold on;
        plot(X(:,1),X(:,2),'ko-');
        plot(Org(1),Org(2),'b*-');
        plot(X(iX,1),X(iX,2),'r*');
    end;
    
    %Loop over all line segments formed by consecutive points on the
    %contour
    tPrm=nan(nX,2);
    for jX=1:nX;
        % Start point
        A=X(jX,:)'-Org;
        if jX<nX;
            % Increment
            D=Dlt(jX,:)';
        else;
            % Increment
            D=Dlt(1,:)';
        end;
        if jX~=iX-1 && jX~=iX; % Cannot match with self
            % This solves for the intersection between two vectors
            tPrm(jX,:)=[-D(1) C(1);-D(2) C(2)]\A;
        end;
    end;
    % Find any line segments that make the point a valid point to drop
    if sum(tPrm(:,1)>0 & tPrm(:,1)<1 & tPrm(:,2)>0 & tPrm(:,2)<1)>0;
        t=find(tPrm(:,1)>0 & tPrm(:,1)<=1 & tPrm(:,2)>0 & tPrm(:,2)<1);
        if Vrb==1;
            plot(X(t,1),X(t,2),'g*');
            fprintf(1,'%g Hit',Idn(iX));
            fprintf(1,' %g',Idn(t==1));
            fprintf(1,'\n');
        end;
        % Drop the offending point in Kep
        Kep(Idn(iX))=0;
        % Break the loop because you've found a point to drop
        break;
    end;
    
end;

return;