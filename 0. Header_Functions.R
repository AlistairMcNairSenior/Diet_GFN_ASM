
# Function to return a palette of hot cold colors that I think look good

cold.hot.colors<-function(n.col){
	
	pallette<-rainbow(n.col, start = 3/6, end = 4/6, s = seq(0.005, 1, (0.995 / n.col))[1:n.col])
	cold.cols<-pallette[length(pallette):1] 

	pallette<-rainbow(n.col, start = 0, end = 1/6, s = seq(1, 0.005, -(0.995 / n.col))[1:n.col])
	warm.cols<-pallette[length(pallette):1] 

	total.cols<-c(cold.cols, warm.cols)
	pallette<-total.cols[seq(1, length(total.cols), 2)]
	
	return(pallette)

}

# Function to return a paleltte with color similar to the classic colors from the fields package

fields.colors<-function (n, alpha = 1) 
{
    if ((n <- as.integer(n[1L])) > 0) {
        j <- n%/%3
        k <- n%/%3
        i <- n - j - k
        c(if (i > 0) hsv(h = seq.int(from = 43/60, to = 30/60, 
            length.out = i), alpha = alpha), if (j > 0) hsv(h = seq.int(from = 12/60, 
            to = 9/60, length.out = j), alpha = alpha), if (k > 
            0) hsv(h = seq.int(from = 5/60, to = 0, length.out = k), 
            alpha = alpha, 
            v = 1))
    }
    else character()
}



# A function created to find the outer perimeter over which the surface should be fitted

findConvex<-function(x,y,rgnames,res=101){
	hull<-cbind(x,y)[chull(cbind(x,y)),]
	px<-pretty(x)
	py<-pretty(y)
	x.new<-seq(min(px),max(px),len=res)
	y.new<-seq(min(py),max(py),len=res)
	ingrid<-as.data.frame(expand.grid(x.new,y.new))                                                              
	Fgrid<-ingrid
	Fgrid[(point.in.polygon(ingrid[,1], ingrid[,2], hull[,1],hull[,2])==0),]<-NA
	names(Fgrid)<-rgnames
	return(Fgrid)
}

# Function to find isocaloric vectors across a specfied polygon of intakes in nutrient space

findIso<-function(intake.x, intake.y, iso.x, iso.y){
	
	intake.x<-summary.stats[,surfaces.list[[i]][1]]
	intake.y<-summary.stats[,surfaces.list[[i]][2]]
	
	point.in.polygon(point.x = iso.vec[,1], point.y = iso.vec[,2], pol.x = hull[,1], pol.y = hull[,2])
}


# Funtion to return the brass coefficients for a specific diet reltaive to the standard
# The regression coeffs are for y = a + b*x, where x is the logit survival on the standard schedule, and y gives the corresponding survival on the specific diet

brass_coefficients<-function(ages_diet, ages_standard, return_standard = F, return_fitted = F){
	
	# Needs the arm package
	require(arm)
	require(survival)
	require(mgcv)
	
	# Generate the standard survival schedule
	obj<-survfit(Surv(ages_standard) ~ 1)
	
	# logit transform survival
	ls<-logit(obj$surv)
	x<-obj$time
	
	# remove the maximal values which generate inf
	ls<-ls[-length(ls)]
	x<-x[-length(x)]
	
	# get a smooth curve for the schedule
	model<-gam(ls ~ s(x))
	
	# now predict lsx for each timepoint over the whole age range
	x_new<-seq(min(x), max(ages_standard, na.rm=T), 1)
	lsx<-predict(model, newdata=data.frame(x=x_new))
	
	# Create the standard
	standard.schedule<-data.frame(age = x_new, logit_lsx = lsx, lsx = invlogit(lsx))
	
	# Now get the patchy schedule for the specific diet
	obj<-survfit(Surv(ages_diet) ~ 1)
	
	# Again logit transform etc
	ls<-logit(obj$surv)
	x<-obj$time
	ls<-ls[-length(ls)]
	x<-x[-length(x)]
	
	# match times to standard
	pred<-standard.schedule$logit_lsx[match(x, standard.schedule$age)]
	
	# fit the model to get the coefficients	
	model<-lm(ls ~ pred)
	
	# Return the coefficients and SEs
	out<-summary(model)$coef[,c(1:2)]
	row.names(out)<-c("alpha", "beta")
	
	if(return_standard == T){
		out<-standard.schedule
	}
	
	if(return_fitted == T){
		predicted<-invlogit(model$coef[1] + model$coef[2] * pred)
		out<-cbind(invlogit(ls), predicted)
		out<-as.data.frame(out)
		names(out)<-c("Observed", "Predicted_Model")
	}
	
	
	return(out)

}


#---------------------------------
#---- special functions
#---------------------------------
#--- function for making subdiagonals (from Bill Venables)
  subdiag <- function (v, k) {
    n <- length(v) + abs(k)
    x <- matrix(0, n, n)
    if (k == 0)
        diag(x) <- v
    else if (k < 0)
      { ## sub-diagonal
        j <- 1:(n+k)
        i <- (1 - k):n
        x[cbind(i, j)] <- v
} 
x } 

# Function to return sensitivity of index to mortality - lifted from van Raalte and Caswell 2013 Demograph 50(5):1615-1640 Supp matts

return_sensitivity<-function(q, x, type){
	
	require(Matrix)
	require(signal)
	
	#------------------------------------------
	# Preliminaries
	#------------------------------------------
	# defining survival probabilities
	  p <- 1-q
	# number of transient states + 1, (in construction of U first row is zero)
	  s <- length(p)+1
	# number of transient states only
	  s2 <- length(p)
	# other things we will need later
	  I <- diag(rep(1,s))  # identity matrix
	  I <- as(I,"sparseMatrix")
	  e <- rep(1,s) # vector of ones for summations
	  e1 <- c(1,rep(0,s-1))
	  age <- 0:s2
	# C matrix is for calculating cumulative sums
	  C <- Matrix(0,nrow=s,ncol=s)
	    for (i in 1:s){
	      C[,i] <- c(rep(0,i),rep(1,s-i))
	      }
	  C <- as(C,"sparseMatrix")
	#--------------------------------------------------
	# Markov chain formulation of longevity
	#--------------------------------------------------
	# U matrix describes transient states in the Markov chain
	  U <- subdiag(p,-1)
	  U <- as(U,"sparseMatrix")
	# N matrix, where (i,j) entry is the mean time spent in each age class i,
	# conditional upon starting in age class j
	  N <- solve(I-U)
	  N <- as(N,"sparseMatrix")
	# M matrix has probability of death at each age on diagonal
	  M <- diag(c(q,1))
	  M <- as(M,"sparseMatrix")
	# The distribution of ages at death
	  B <- M %*% N  # the complete distribution of ages at death
	  f <- B %*% e1 # the distribution of ages at death from birth (1st age class)
	  B <- as(B,"sparseMatrix")
	# survivorship (alternatively ell <- e-C%*%f )
	  ell <- N %*% e1
	  ell <- as(ell,"sparseMatrix")
	# remaining life expectancy at each age
	  mean_eta <- colSums(N)-0.5
	# life expectancy at birth (or first age class)
	  eta <- mean_eta[1]
	# NB: in Markov chain formulations, the life expectancy at birth is always
	# 0.5 years higher than that found by conventional life table methods,
	# which is why we subtract 0.5 years
	#------------------------------------
	# Indices of lifespan variation
	#------------------------------------
	# variance in lifespan
	  V <- colSums(N) %*% (2*N-I) - mean_eta*mean_eta
	# standard deviation in lifespan
	  S <- sqrt(V)	
	
	# derivatives of U with respect to mortality change
	   dvecU_dmu <- Matrix(0,nrow=s*s,ncol=s-1)
	   r <- seq(2,s*s,s+1) # rows that will contain the elements once stacked
	     for (i in 1:(s-1)){
	      dvecU_dmu[r[i],i] <- -p[i]
	      }
	   dvecU_dmu <- as(dvecU_dmu,"sparseMatrix")
	# derivatives of M with respect to mortality change
	   dvecM_dmu <- Matrix(0,nrow=s*s,ncol=s-1)
	   r2 <- seq(1,s*s,s+1) # rows that will contain the elements once stacked
	    for (i in 1:(s-1)){
	     dvecM_dmu[r2[i],i] <- p[i]
	     }
	   dvecM_dmu <- as(dvecM_dmu,"sparseMatrix")
	# derivative of f with respect to mortality change
	   df_dmu <- t(kronecker(N[,1],I)) %*% dvecM_dmu +
	     t(kronecker(N[,1],t(B))) %*% dvecU_dmu
	# sensitivity of expected longevity with respect to mortality change
	   deta_dmu <- t(kronecker(N[,1],colSums(N))) %*% dvecU_dmu
	# sensitivity of variance in longevity with respect to mortality change
	   dV_dmu <- (2*kronecker(t(N) %*% t(N),t(mean_eta))  +
	     2*kronecker(t(N),mean_eta %*% N) -
	     (I + 2*(diag(mean_eta))) %*% kronecker(t(N),t(mean_eta))) %*% dvecU_dmu
	# sensitivity of standard deviation with respect to mortality change
	   dS_dmu <- 0.5*diag(as.vector(1/S)) %*% dV_dmu	
	
	out<-list()
	  
	  if(type=="S"){
	  	out[[1]]<-S[1]
	  	out[[2]]<-dS_dmu[1,]
	  }
	  	  
	  if(type == "e0"){
	  	out[[1]]<-eta
	  	out[[2]]<-as.numeric(deta_dmu[1,])
	  }
	  
	  names(out)<-c("index", "sensitvity")
	  
	  return(out)    
    
}    

# Function to return and index calculated from life table parameters - it is a subset of the above code 

return_index<-function(q, x, type){
	
	#------------------------------------------
	# Preliminaries
	#------------------------------------------
	# defining survival probabilities
	  p <- 1-q
	# number of transient states + 1, (in construction of U first row is zero)
	  s <- length(p)+1
	# number of transient states only
	  s2 <- length(p)
	# other things we will need later
	  I <- diag(rep(1,s))  # identity matrix
	  I <- as(I,"sparseMatrix")
	  e <- rep(1,s) # vector of ones for summations
	  e1 <- c(1,rep(0,s-1))
	  age <- 0:s2
	# C matrix is for calculating cumulative sums
	  C <- Matrix(0,nrow=s,ncol=s)
	    for (i in 1:s){
	      C[,i] <- c(rep(0,i),rep(1,s-i))
	      }
	  C <- as(C,"sparseMatrix")
	#--------------------------------------------------
	# Markov chain formulation of longevity
	#--------------------------------------------------
	# U matrix describes transient states in the Markov chain
	  U <- subdiag(p,-1)
	  U <- as(U,"sparseMatrix")
	# N matrix, where (i,j) entry is the mean time spent in each age class i,
	# conditional upon starting in age class j
	  N <- solve(I-U)
	  N <- as(N,"sparseMatrix")
	# M matrix has probability of death at each age on diagonal
	  M <- diag(c(q,1))
	  M <- as(M,"sparseMatrix")
	# The distribution of ages at death
	  B <- M %*% N  # the complete distribution of ages at death
	  f <- B %*% e1 # the distribution of ages at death from birth (1st age class)
	  B <- as(B,"sparseMatrix")
	# survivorship (alternatively ell <- e-C%*%f )
	  ell <- N %*% e1
	  ell <- as(ell,"sparseMatrix")
	# remaining life expectancy at each age
	  mean_eta <- colSums(N)-0.5
	# life expectancy at birth (or first age class)
	  eta <- mean_eta[1]
	# NB: in Markov chain formulations, the life expectancy at birth is always
	# 0.5 years higher than that found by conventional life table methods,
	# which is why we subtract 0.5 years
	#------------------------------------
	# Indices of lifespan variation
	#------------------------------------
	# variance in lifespan
	  V <- colSums(N) %*% (2*N-I) - mean_eta*mean_eta
	# standard deviation in lifespan
	  S <- sqrt(V)
	  
	  if(type=="S"){
	  	out<-S[1]
	  }
	  
	  if(type == "e0"){
	  	out<-eta
	  }
	  
	  return(out)
	  
}

# Function to copnvert brass coefficients life expectancy or SD in age at death based on a standard pattern

convert_brass<-function(alpha, beta, lsx, type = "e0", age1 = 0, return_lt = F){
	
	require(arm)
	
	if(is.na(alpha) == T | is.na(beta) == T){
		out<-NA
	}else{
			
		# Calculate the life table based on Brass Coeffs
		lx<-c(1, invlogit(alpha + beta * logit(lsx)))
		px<-c(lx[-1], 0) / lx
		x<-seq(age1, age1+length(lx)-1, 1)
		qx<-(1 - px)
		lt<-data.frame(age = x, lx, qx, px)
		rownames(lt)<-seq(1, length(lx), 1)
		
		# From supp mats van Raalte and Caswell 2013 - convert q values to other various stats
		
		q<-lt$qx
		x<-lt$age+0.5
	  	x<-c(x, max(x)+1)
		out<-return_index(q = q, x = x, type = type)
			
		if(return_lt == T){
			out<-lt
		}	
	}
		
	return(out)
}

