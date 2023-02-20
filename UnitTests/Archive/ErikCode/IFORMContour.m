function O=IFORMContour(varargin)
% Traditional I-FORM contours.
%%%%% NB! Assuming the conditional model of 3-p Weibull and log-normal. With other models, the transformations become different.
%%%%% To change axis (putting Tz on the x-axis), set rotate to T
%
% INPUT
% (optional)
% pf: default  (0.01,0.001,0.0001)
% n: default 60  %number of points to compute contour at?
% discarded: default 0 conatins the number of discarded data-points if importance sampling has been applied
% theta:  default [1.471,2.776,0.8888] weibull parameter values shape, scale, thres
% a:  default [0.1,1.489,0.1901] lognormal mean function parameters values a(1) + a(2).*X ^ a(3);
% b: default [0.04,0.1748,-.2243] lognormal std function parameters values b(1) + b(2).*exp(b(3) X);
% main:  plot title default 'IFORM contours'
% parname: parameter names for plot
% rotate:  flag to make rotate x y axis to put Hs on Y axis
% plotOn:  flag to make plot
%
%OUTPUT
%
%O.xy = xy   %xy points of contour
%O.type='Huseby';  %contour type

%% Parse inputs
% TODO add validate arguement1
p = inputParser;
addOptional(p,'pf',[0.01,0.001,0.0001]); %probability levels
addOptional(p,'n',60);  %number of angles
addOptional(p,'theta',[1.471,2.776,0.8888]) %weibull shap, scale, thres
addOptional(p,'a',[0.1,1.489,0.1901]) %conditioal lognormal mean function a(1) + a(2).*X ^ a(3);
addOptional(p,'b',[0.04,0.1748,-.2243])  %conditioal lognormal std function b(1) + b(2).*exp(b(3) X);
addOptional(p,'main','IFORM contours') %plot title
addOptional(p,'parname',{'Hs','Tz'})
addOptional(p,'rotate',false)
addOptional(p,'plotOn',true)

parse(p,varargin{:});
IN=p.Results;

nPrb=numel(IN.pf); %numer of probability levels

angles = linspace(0, 2*pi, IN.n)';
Prj=[cos(angles), sin(angles)];
%log normal conditional function
cond_function_mu = @(x)  IN.a(1) + IN.a(2)*x.^IN.a(3);
cond_function_sigma = @(x) IN.b(1) + IN.b(2)*exp(IN.b(3)*x);


u =NaN(IN.n,2,nPrb);
xy=NaN(IN.n,2,nPrb);

for i=1:nPrb
    u(:,:,i) = norminv(1-IN.pf(i)).*Prj;
end
U=normcdf(u);

xy(:,1,:) = wblinv(U(:,1,:),IN.theta(2),IN.theta(1))+IN.theta(3);
xy(:,2,:) = logninv(U(:,2,:), cond_function_mu(xy(:,1,:)), cond_function_sigma(xy(:,1,:)));
if IN.rotate
    xy = xy(:,[2,1],:)  ;
    IN.parname=IN.parname([2,1]);
end %rotate if
if IN.plotOn
    C=lines(nPrb);
    
    for i=1:nPrb
        plot(xy(:,1,i),xy(:,2,i),'color',C(i,:), 'linewidth',2)
        hold on
    end
    
    xlabel(IN.parname(1))
    ylabel(IN.parname(2));
    
    legend(num2str(IN.pf(:)),'location','northwest')
end %plot if
O.xy=xy;
O.type='IFORM';

end %IFORMContour
