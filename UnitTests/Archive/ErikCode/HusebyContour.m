function O=HusebyContour(x, varargin)
%O=HusebyContour(x, varargin) Drawing Huseby contours from a bivariate sample
%
% INPUT
% x: should be on the format of a matrix with two columns, one for each variable
% (optional)
% pf: default  [0.01,0.001,0.0001]
% n: default 60  %number of points to compute contour at?
% Standardize: default false
% scatter: default false-
% discarded: default 0 conatins the number of discarded data-points if importance sampling has been applied
% main:  plot title default ''Environmental contours'
% plotOn:  flag to make plot
%
%OUTPUT
%
%O.C_theta = C_theta   %??
%O.xy = xy   %xy points of contour
%O.chull = chull   %convex hull
%O.type='Huseby';  %contour type

%% Parse inputs
% TODO add validate arguement1
p = inputParser;
addOptional(p,'pf',[0.01,0.001,0.0001]);
addOptional(p,'n',60);  %number of angles
addOptional(p,'Standardize',false);
addOptional(p,'scatter',false);
addOptional(p,'discarded',0)
addOptional(p,'main','Environmental contours')
addOptional(p,'parname',{'Hs','Tz'})
addOptional(p,'rotate',false)
addOptional(p,'plotOn',true)

parse(p,varargin{:});
IN=p.Results;
  
%% use bivariate standardization if Standardize = TRUE
if IN.Standardize      
    y = x;  % Storing the the data on the original scale
    mu = mean(x, 2);  %mean
    sigma = std(x, 2);   %standard deviation
    rho=corr(x); %correlation

    x(:,1) = (y(:,1)-mu(1))./(sigma(1));
    x(:,2) = ((y(:,2)-mu(2))./sigma(2) - rho*x(:,1))./sqrt(1-rho.^2);
end


k = ceil((length(x(:,1))+IN.discarded)*(1-IN.pf)-IN.discarded); %number of points above probiblity level?
angles = linspace(0, 2*pi, IN.n); %set of angles

%% Calculate the function C_theta
C_theta=NaN(IN.n, numel(k)); %preallocate
for i=1:length(angles)
    X = x(:,1)*cos(angles(i)) + x(:,2)*sin(angles(i));
    X=sort(X);
    C_theta(i,:) = X(k);
end

%% Calculate the intersection points
xy=NaN(IN.n,2,length(IN.pf));
for i =1:IN.n-1
    for j=1:numel(k)
        xy(i, 1, j) = (sin(angles(i+1))*C_theta(i, j)-sin(angles(i))*C_theta(i+1, j))./(sin(angles(i+1))*cos(angles(i))-sin(angles(i))*cos(angles(i+1)));
        xy(i, 2, j) = (-cos(angles(i+1))*C_theta(i, j)+cos(angles(i))*C_theta(i+1, j))./(sin(angles(i+1))*cos(angles(i))-sin(angles(i))*cos(angles(i+1)));
    end
end
xy(end,:,:) = xy(1,:,:); %final point is start of contour again

%% If standardization has been applied, transfer the intersection points back to the original space
if IN.Standardize                
    x = y;
    xy_stand = xy;                  %% Storing the normalized intersection points
    xy(:,1,:) = sigma(1).*xy_stand(:,1,:)+mu(1);
    xy(:,2,:) = sigma(2).*(rho.*xy_stand(:,1,:) + xy_stand(:,2,:).*sqrt(1-rho^2)) + mu(2);
end

if IN.rotate
    xy = xy(:,[2,1],:)  ;
    IN.parname=IN.parname([2,1]);
end %rotate if

%% Calculate the convex hull of all intersection points
chull = cell(length(IN.pf),1);
for i=1:length(IN.pf)
     chull{i} = convhull(xy(:,1,i),xy(:,2,i));
 end

%% Make plot


C=lines(length(k));
if IN.plotOn
    clf; hold on;
    if IN.scatter
        plot(x(:,1),x(:,2), 'k.','handlevisibility','off')            
    end        
    for i=1:length(k)
        plot(xy(:,1,i),xy(:,2,i),'.', 'color',C(i,:),'linewidth',2,'markersize',10,'handlevisibility','off')
        hold on
        plot(xy(chull{i},1,i),xy(chull{i},2,i),'-', 'color',C(i,:),'linewidth',2)
    end
    title(IN.main)
    xlabel(IN.parname(1));
    ylabel(IN.parname(2));
    axis tight
    legend(num2str(IN.pf(:)),'location','northwest')
end

%% Output
O.C_theta = C_theta;
O.xy = xy;
O.chull = chull;
O.type='Huseby';

