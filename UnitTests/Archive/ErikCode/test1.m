% Dear all, 
% 
% As promised, please find attached my R-scripts for calculating the contours. Arne’s code is presumably much more professional, and this has been more for “personal use”. There are probably much more elegant ways of implementing this. But this one is in R and seems to work. I should also note that these contours are not calculated in the exact same way as Arne’s implementation. I have implemented Method 1, whereas Arne has used Method 2 in his code (with references to our Structural Safety paper). However, differences should be small with large enough sample size. 
% 
% The Huseby-contours (in lack of a better name when I made this) takes as input a bivariate sample. In addition the exceedance probabilities should be specified (pf), number of tangent lines to calculate (n), and whether or not to standardize the variables (Standardize). If importance sampling has been used to generate the input sample, the number of discarded samples should be specified (discarded), and you can decide whether to plot the contours (plot) and whether to include the scatter (scatter). You may also change the title of the plot (main). 
% 
% The output is a list that includes the estimated C-function (C_theta), the intersection points that makes up the contours (xy) and the smallest convex hull that contains all the contour-points (chull – used to get smooth contours). In addition the type of contour is kept. The resulting contours that could be plotted or used in subsequent analyses are either “xy” or “chull”. IF sufficient number of samples have been used compared to the exceedance probability, these should be quite similar. 
% 
% I also include my script for the IFORM contours. This particular script assumes a 3-p Weibull and a conditional log-normal, but of course the IFORM contours can be used for any other distribution in principle. The input here are the exceedance probabilities and the distribution parameters as well as the number of points where the contours are evaluated (uniformly distributed between 0 and 360). The output is on the same format as the Huseby-contours, except that there is no C-function and no convex hull. 
% 
% Comments and questions are welcome if you want to try out the functions. 
% 
% Erik 


%%%%% Drawing Huseby contours from a bivariate sample
%%%%% X should be on the format of a matrix with two columns, one for each variable
%%%%% NB: should include colnames for the variables
%%%%% Discarded conatins the number of discarded data-points if importance sampling has been applied

clear; clc; 

nRls=100000;

if 1
    Mu=[0,0];
    Rho=0.5;
    Sigma=[1,Rho;Rho,1];
    Z=exp(-mvnrnd(Mu, Sigma,nRls));   
else
    
    Alp=0.4;
    Z=rndsymlgs(Alp,2,nRls)';
    Z=log(Z); %frechet to gumbel margins
end
figure(1);
clf;
O=HusebyContour(Z,'scatter',true,'n',60);


% % savePics('HusebyContour', 8, 8*2/3,'jpg');    % Height in inches)
% % 
% figure(2);
% clf;
% O2=IFORMContour('rotate',true);
% savePics('IFORMExample', 8, 8*2/3,'jpg');    % Height in inches)
% 
% % 

