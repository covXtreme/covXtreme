##### Drawing Huseby contours from a bivariate sample
##### X should be on the format of a matrix with two columns, one for each variable
##### NB: should include colnames for the variables
##### Discarded conatins the number of discarded data-points if importance sampling has been applied

HusebyContour = function(x, pf = c(0.01, 0.001, 0.0001), n=60, Standardize=F, scatter=F, discarded=0, main="Environmental contours", plot=T) {
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
  dimnames(xy)[[3]] = paste("Pf = ", format(pf, scientific = TRUE, digits = 3))
  
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
  #### Calculate the convex hull of all intersection points
  chull = list()
  for(i in 1:length(pf)){
	chull[[i]] = xy[c(chull(xy[,,i]), chull(xy[,,i])[1]),,i]
  }
  if(plot){
	if(scatter) {plot(x, col="grey", main = main)} else 
	    plot(NULL, xlim=range(x[,1]), ylim=range(x[,2]), xlab=colnames(x)[1], ylab=colnames(x)[2], main=main)
  	for(i in 1:length(k)) {
    		lines(xy[,,i], col=i+1, lwd=2)
  	}
  	legend("topleft", dimnames(xy)[[3]], col=seq(2, i+1) , text.col=seq(2, i+1), bty="n")
  }
list(C_theta = C_theta, xy = xy, chull = chull, type="Huseby")
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
	dimnames(xy)[[3]] = paste("Pf = ", format(Pf, scientific = TRUE, digits = 3))
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




