##### Drawing Huseby contours from a bivariate sample
##### X should be on the format ov a matrix with two columns, one for each variable
##### NB: should include colnames for the variables
##### Discarded conatins the number of discarded data-points if importance sampling has been applied

HusebyContour = function(x, pf = c(0.01, 0.001, 0.0001), n=60, Standardize=F, scatter=F, discarded=0, main="Environmental contours") {
  if(identical(colnames(x),NULL)) {colnames(x) = c("v1", "v2")}
  mu = apply(x, 2, mean); sigma = apply(x, 2, sd); rho=cor(x)[1,2]
  if(Standardize) {      ### use bivariate standardization if Standardize = TRUE
    y = x                ### Storing the the data on the original scale
    mu = apply(x, 2, mean); sigma = apply(x, 2, sd); rho=cor(x)[1,2]
    x[,1] = (y[,1]-mu[1])/(sigma[1])
    x[,2] = ((y[,2]-mu[2])/sigma[2] - rho*x[,1])/sqrt(1-rho^2)
  }
  k = ceiling((length(x[,1])+discarded)*(1-pf)-discarded)   
  delta=2*pi/n
  angles = seq(0, 2*pi, by=delta)
  C_theta=array(NA, dim=c(length(angles), length(k)), dimnames=list(seq(0, 360, by=360/n), pf))

  for(i in 1:length(angles)){                   ### Calculate the function C_theta
    X = x[,1]*cos(angles[i]) + x[,2]*sin(angles[i])
    C_theta[i,] = sort(X)[k];  
    }

  xy=replicate(length(pf), matrix(NA, ncol=2, nrow=length(angles), dimnames=list(seq(0, 360, by=360/n), colnames(x))))
  dimnames(xy)[[3]] = paste("Pf = ", pf)
  
  for(i in 1:(length(angles)-1)) {    ##### Calculate the intersection points
    for(j in 1:length(k)) {
      xy[i, 1, j] = (sin(angles[i+1])*C_theta[i, j]-sin(angles[i])*C_theta[i+1, j])/(sin(angles[i+1])*cos(angles[i])-sin(angles[i])*cos(angles[i+1]))
      xy[i, 2, j] = (-cos(angles[i+1])*C_theta[i, j]+cos(angles[i])*C_theta[i+1, j])/(sin(angles[i+1])*cos(angles[i])-sin(angles[i])*cos(angles[i+1]))
    }
  }
  xy[length(angles), , ] = xy[1,,]

  if(Standardize) {               ### If standardization has been applied, transfer the intersection points back to the original space
    x = y
    xy_stand = xy                  ### Storing the normalized intersection points
    xy[,1,] = (sigma[1]*xy_stand[,1,])+mu[1]
    xy[,2,] = sigma[2]*(rho*xy_stand[,1,] + xy_stand[,2,]*sqrt(1-rho^2)) + mu[2]
  }
  if(scatter) {plot(x, col="grey", main = main)} else 
    plot(NULL, xlim=range(x[,1]), ylim=range(x[,2]), xlab=colnames(x)[1], ylab=colnames(x)[2], main=main)
  for(i in 1:length(k)) {
    lines(xy[,,i], col=i+1, lwd=2)
  }
  legend("topleft", dimnames(xy)[[3]], col=seq(2, i+1) , text.col=seq(2, i+1), bty="n")
list(C_theta = C_theta, xy = xy, type="Huseby")
}

### To change the axes and rotate the contours, change Tzaxis from 1 to 2

# Plotting a Contour object that has already been calculated
# May be either type "IFORM" or "Huseby"
PlotContour = function(cont, add=F, main="Environmental contours", col.start=2, levels=NULL, lty=1, Tzaxis=1) {       
  points=cont$xy
  if(Tzaxis==2) points = points[,c(2,1),]
   if(is.null(levels)) {					###If levels are not specified - plot all levels
	levels = seq(1: length(points[1,1,]))     
  }
  col = seq(col.start, col.start+length(levels)-1)		### The colours to use for plotting
  if(!add) {
    plot(NULL, xlim=range(points[,1,]), ylim=range(points[,2,]), xlab=colnames(points)[1], ylab=colnames(points)[2], main=main)
    legend("topleft", dimnames(points)[[3]][levels], col=col, text.col=col, bty="n")
  }
  j=1
  for(i in levels) {
    lines(points[,,i], col=col[j], lwd=2, lty=lty)
    j=j+1
  }
	legend("bottomright", cont$type, bty="n", lty=lty, lwd=2)
}

### Function to compare two sets of contours - typically estiamted by IFORM and Huseby methods, respectively. 
### NB: Should have same exceedance probabilities. Otherwise, comparison is not valid and an error message is returned
### Tzaxis can change the exis of each of the sets of contours
### Levels refers to which Pf-levels to plot, not hte actual level
### ExtUpperRange refers to the factor to extend the upper limits of the x- and y-axes (default to no extension)

CompareContours = function(cont1, cont2, levels=NULL, Tzaxis=c(1,1), col.start=2, main="Environmental contours", ExtUpperRange=1) {
 ### Check if sets of contours are compareable
  points1=cont1$xy; points2=cont2$xy
  if(!all(dimnames(points1)[[3]]==dimnames(points2)[[3]])) warning("Exceedance probabilities does not match")
  if(Tzaxis[1]==2) points1 = points1[,c(2,1),]
  if(Tzaxis[2]==2) points2 = points2[,c(2,1),]
   if(is.null(levels)) {					###If levels are not specified - plot all levels
	levels = seq(1: length(points1[1,1,]))     
  }
  col = seq(col.start, col.start+length(levels)-1)		### The colours to use for plotting
    plot(NULL, xlim=c(range(points1[,1,], points2[,1,])[1], range(points1[,1,], points2[,1,])[2]*ExtUpperRange), ylim=c(range(points1[,2,], points2[,2,])[1], range(points1[,2,], points2[,2,])[2]*ExtUpperRange), xlab=colnames(points1)[1], ylab=colnames(points2)[2], main=main)
    legend("topleft", dimnames(points1)[[3]][levels], col=col, text.col=col, bty="n")
  j=1
  for(i in levels) {
    lines(points1[,,i], col=col[j], lwd=2)
    lines(points2[,,i], col=col[j], lwd=2, lty=2)
    j=j+1
  }
	legend("bottomright", c(cont1$type, cont2$type), bty="n", lty=c(1, 2), lwd=2)
}


########## sampling from a conditional model, as specified in the paper ##########################################
########## Assuming a marginal 3-parameter Weibull and a conditional log-normal with a, b parameters #############

Simulate_ConditionalModel = function(theta = list(shape=1.471, scale=2.776, thres=0.8888), a = c(0.1000, 1.489, 0.1901), b =c(0.0400, 0.1748, -0.2243) , N = 100000, parname=c("Hs", "Tz")) {
  require(FAdist)
  cond_function_mu = function(x, a1, a2, a3)  a1 + a2*x^a3
  cond_function_sigma = function(x, b1, b2, b3) b1 + b2*exp(b3*x)
  SimData = matrix(ncol=4, nrow=N)     ### Data with four columns: Hs, conditional function for mu, conditional function for sigma and Tz
  colnames(SimData)=c(parname[1], "mu(h)", "sigma(h)", parname[2])
  SimData[,1] = rweibull3(N, shape=theta[["shape"]], scale=theta[["scale"]], thres=theta[["thres"]])
  for(i in 1:N) {
    SimData[i,2] = cond_function_mu(SimData[i,1], a[1], a[2], a[3])
    SimData[i,3] = cond_function_sigma(SimData[i,1], b[1], b[2], b[3])
    if(SimData[i,3]>0)          ### Checking if the conditional standard deviation is positive and non-zero. Otherwise, a small positive sigma is assigned
      SimData[i,4] = rlnorm(1, meanlog=SimData[i,2], sdlog=SimData[i,3])
    else
      SimData[i,4] = rlnorm(1, meanlog=SimData[i,2], sdlog=0.0001)
    }

  SimData[, c(1, 4)]  
}

### Importance sampling in order to improve the smoothness of the contours

ImportanceSampling_ConditionalModel = function(theta = list(shape=1.471, scale=2.776, thres=0.8888), a = c(0.1000, 1.489, 0.1901), b =c(0.0400, 0.1748, -0.2243) , N = 100000, Pd = 0.99, 
							parname = c("Hs", "Tz")) {
  require(FAdist)
  sample=matrix(ncol=2, nrow=N); colnames(sample) = parname
  k=0         ### Keeping track of how many datapoints that have been discarded
  n=0         ### Keeping track of how many datapoints have been kept
  while(n <= N) {
    h = rweibull3(1, shape=theta[["shape"]], scale=theta[["scale"]], thres=theta[["thres"]])
    cond_mu = a[1] + a[2]*h^a[3]
    cond_sigma = b[1] + b[2]*exp(b[3]*h)
    if(cond_sigma<=0) cond_sigma=0.0001          ### Checking if the conditional standard deviation is positive and non-zero. Otherwise, a small positive sigma is assigned
    t = rlnorm(1, meanlog=cond_mu, sdlog=cond_sigma)
    #### Rosenblatt transform into standard normal space (h, t) -> (u1, u2)
    u1 = qnorm(pweibull3(h, shape=theta[["shape"]], scale=theta[["scale"]], thres=theta[["thres"]]))
    u2 = qnorm(plnorm(t, meanlog=cond_mu, sdlog=cond_sigma))
    if((u1^2 + u2^2) > qchisq(Pd, 2)){        ### Checking if the sample should be kept according to specified discard probability
      sample[n, ] = c(h, t)
      n=n+1
    } else
      k=k+1      
  }
  o = list(sample=sample, discarded=k)
  o
}

##### Traditional I-FORM contours. 
##### NB! Assuming the conditional model of 3-p Weibull and log-normal. With other models, the transformations become different. 
##### To change axis (putting Tz on the x-axis), set rotate to T

IFormContour_ConditionalModel = function(Pf = c(0.01, 0.001, 0.0001), n=60, main="IFORM contours", parname = c("Hs", "Tz"),
			theta = list(shape=1.471, scale=2.776, thres=0.8888), a = c(0.1000, 1.489, 0.1901), b =c(0.0400, 0.1748, -0.2243),
			rotate = F){
	require(FAdist)
	angles = seq(0, 2*pi, by=2*pi/n)
 	cond_function_mu = function(x)  a[1] + a[2]*x^a[3]
  	cond_function_sigma = function(x) b[1] + b[2]*exp(b[3]*x)
	u =replicate(length(Pf), matrix(NA, ncol=2, nrow=length(angles), dimnames=list(seq(0, 360, by=360/n), c("u1", "u2"))))
	xy = replicate(length(Pf), matrix(NA, ncol=2, nrow=length(angles), dimnames=list(seq(0, 360, by=360/n), parname)))
	dimnames(xy)[[2]] = parname
	dimnames(xy)[[3]] = paste("Pf = ", Pf)
	for(i in 1:length(Pf)){
		u[,,i] = qnorm(1-Pf[i])*cbind(cos(angles), sin(angles))
		}	
	xy[,1,] = qweibull3(pnorm(u[,1,]), shape=theta[["shape"]], scale=theta[["scale"]], thres=theta[["thres"]])
	xy[,2,] = qlnorm(pnorm(u[,2,]), meanlog=cond_function_mu(xy[,1,]), sdlog=cond_function_sigma(xy[,1,]))
 	if(rotate) xy = xy[,c(2,1),]
	plot(NULL, xlim=range(xy[,1,]), ylim=range(xy[,2,]), xlab=dimnames(xy)[[2]][1], ylab=dimnames(xy)[[2]][2], main=main)
  	for(i in 1:length(Pf)) {
    		lines(xy[,,i], col=i+1, lwd=2)
  	}
  	legend("topleft", dimnames(xy)[[3]], col=seq(2, i+1) , text.col=seq(2, i+1), bty="n")
list(xy=xy, type="IFORM")
}

#### Rosenblatt transform of points in X-space to points in U-space
#### NB! Assuming a conditional model of 3-par Weibull and conditional log-normal. Other models gives different transformations
#### Transforms point by point and assumes x is on the form c(x1, x2). 
#### NB! x1 should be Hs and x2 should be Tp
RosenblattTransform_CondMod = function(x, theta = list(shape=1.471, scale=2.776, thres=0.8888), a = c(0.1000, 1.489, 0.1901), b =c(0.0400, 0.1748, -0.2243)){
 	cond_function_mu = function(x)  a[1] + a[2]*x^a[3]
  	cond_function_sigma = function(x) b[1] + b[2]*exp(b[3]*x)
	x1 = x[1]; x2 = x[2]
	u1 = qnorm(pweibull3(x1, shape=theta[["shape"]], scale=theta[["scale"]], thres=theta[["thres"]]))
	u2 = qnorm(plnorm(x2, meanlog=cond_function_mu(x1), sdlog=cond_function_sigma(x1)))
	return = c(u1, u2)
}


