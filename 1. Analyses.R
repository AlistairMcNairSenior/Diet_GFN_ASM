
# Code written by AM Senior @ the Univeristy of Sydney in 2018/2019
# Performs the analyses in Senior et al. 

# Clear the environment
rm(list=ls())

# Load the packages
library(plyr)
library(fields)
library(mgcv)
library(sp)
library(flexsurv)
library(arm)
library(pracma)
library(MALDIquant)
library(lattice)
library(gridExtra)
library(signal)
library(foreach)
library(doParallel)

# Slow things down to let packages load (helps when running in terminal)
Sys.sleep(5)

# Set the wd
setwd("/Users/alistairsenior/Dropbox (Sydney Uni)/Work-Home Share/DECRA/Aim 3/Analyses Main Text")

# Functions written for these analyses
source("0. Header_Functions.R")

# Figure 1 - conceptual figure for effects of alpha and beta on survival in population relative to stnadard 

# Open the file for plotting
pdf("Fig1.pdf", height=10, width=10)

# Set up the plotting region
par(mfrow=c(2,2), mar=c(5,6,5,2))

# Logit survival on the standard - assume simplistic relationship between 6 (close to 100%) and -6 (close to 0%)
lsx<-seq(-6, 6, 1)

# A series of alpha and beta values to plot
a<-c(0, 1, 0, -1)
b<-c(1, 1, 1.5, 0.5)

# Labels for panels
labels<-c("A", "B", "C", "D")

# open loop to plot each series
for(i in 1:4){
	
	# Logit survival in target population
	lx<-a[i] + b[i] * lsx
	
	# Plot the data
	plot(lsx, lx, pch=16, xlab = expression(logit(italic(l[x]^s))), ylab=expression(logit(italic(l[x]))), ylim=c(-6, 6), xlim=c(-6, 6), cex.lab=1.25)
	
	# Add some reference lines
	abline(v=0, lty=2)
	abline(h=0, lty=2)
	abline(a=0, b=1, lty=2)
	
	# Add the values to the panel
	text(-5, 5, substitute(italic(alpha) == a, list(a = a[i])), cex=1.5)
	text(-5, 4, substitute(italic(beta) == b, list(b = b[i])), cex=1.5)
	mtext(labels[i], at=-6.5, line=2, cex=1.5)
	
	# Add some text to axis on panel A
	if(i == 1){
		mtext("survival at
young ages", side=1, at = 4, line=3.5)
		mtext("survival at
old ages", side=1, at = -4, line=3.5)
		mtext("survival at
young ages", side=2, at = 4, line=3)
		mtext("survival at
old ages", side=2, at = -4, line=3)
	}
}

# Close the pdf
dev.off()

# Reading in data and some processing

# Mouse self selected diets is poorly defined for fat Sorenson et al. give an estimate for P:C based on a diet around 7% fat. Solon-biet et al. using the half-life of an assymptotic curve give an estimate of the effect of diet macronutrient composition on feeding. Taking the average of these two estimates we get
p_s<-0.2241
c_s<-0.4706
f_s<-0.3054

# Read in the raw data
raw.data<-read.csv("Raw_Mouse.csv")

# Drop data from the poor diets, which were terminated in the first few weeks
drop.diets<-c("SF09-046", "SF09-049", "SF09-050", "SF09-058", "SF09-059", "STANDARD")
raw.data<-raw.data[which(is.na(match(raw.data$diet, drop.diets)) == T),]

## Get the paramters, alpha and beta from each diet by fitting the brass_coefficients function to each diet, using the whole dataset as a standard
summary.stats<-ddply(raw.data, .(diet), summarise, 
brass.alpha = brass_coefficients(ages_diet = ageatdeath_w, ages_standard = raw.data$ageatdeath_w)[1,1],
brass.beta = brass_coefficients(ages_diet = ageatdeath_w, ages_standard = raw.data$ageatdeath_w)[2,1],
SE.brass.alpha = brass_coefficients(ages_diet = ageatdeath_w, ages_standard = raw.data$ageatdeath_w)[1,2],
SE.brass.beta = brass_coefficients(ages_diet = ageatdeath_w, ages_standard = raw.data$ageatdeath_w)[2,2],
P_kJ.g = P_kJ.g[1], 
C_kJ.g = C_kJ.g[1],
F_kJ.g = F_kJ.g[1])

# Transform SEs to sampling variances
summary.stats$V.brass.alpha<-summary.stats$SE.brass.alpha^2
summary.stats$V.brass.beta<-summary.stats$SE.brass.beta^2

# Write the coefficients to a table
write.table(summary.stats, file="Table.S1.csv", sep=",", row.names=F, col.names=names(summary.stats))

# Figure 2 - visualise the standard life table and check how the model lifetables fit the data

# Open the file for plotting
pdf("Fig2.pdf", height=5, width=10)

# Set up the plotting region
par(mfrow=c(1,2), mar=c(5,5,3,2))

# survival curve for whole dataset
curve<-survfit(Surv(raw.data$ageatdeath_w) ~ 1)
plot(curve, xlab="Age (weeks)", lwd=1.25)
mtext(expression(italic(l[x])), side=2, line=2, cex=2)

# Get the logit survival back for the whole population and plot
standard<-brass_coefficients(ages_diet = raw.data$ageatdeath_w, ages_standard = raw.data$ageatdeath_w, return_standard = T)
lines(standard$age, standard$lsx, col=2, lwd=1.5, lty=2)
mtext("A", at=-1, line=1, cex=2)

# The model on the dashed line (same as is fitted internally in the brass_coefficients function)
outcome<-logit(curve$surv)[-length(curve$surv)]
predictor<-curve$time[-length(curve$surv)]
GAM<-gam(outcome ~ s(predictor))
summary(GAM)

# Fit the brass model to each diet individually and return the observed vs predicted
diets<-unique(raw.data$diet)
diet<-raw.data[which(raw.data$diet == diets[1]),]
results<-brass_coefficients(ages_diet = diet$ageatdeath_w, ages_standard=raw.data$ageatdeath_w, return_fitted=T)
results$diet<-diets[1]
for(i in 2:length(diets)){
	diet<-raw.data[which(raw.data$diet == diets[i]),]
	r_i<-brass_coefficients(ages_diet = diet$ageatdeath_w, ages_standard=raw.data$ageatdeath_w, return_fitted=T)
	r_i$diet<-diets[i]
	results<-rbind(results, r_i)
}

# Plot the predicted vs observed
plot(results$Observed, results$Predicted_Model, xlab=expression(Observed ~ italic(l[x])), ylab=expression(Predicted ~ italic(l[x])), cex.lab=1.2, pch=16, cex=0.5, ylim=c(0, 1), xlim=c(0, 1))
mtext("B", at=-0.01, line=1, cex=2)
abline(a=0, b=1)

# Add on the r squared
r2<-round(cor(results$Observed, results$Predicted_Model)^2, 2)
text(0.1, 1, substitute(italic(r^2) == r2, list(r2 = r2)))

# Close the plotting file
dev.off()

# Figure 3 - GAM surfaces for alpha and beta

# File for plotting
pdf("Fig3.pdf", height=13, width=15)

# File to write results from each GAM too
res.file<-"GAM_LS.csv"

# Set up the plotting region
par(mar=c(7,5,11,1), mfrow=c(2, 3))	

# Specify the model formula to use (passed to GAM function)
model.form<-as.formula("this.outcome ~ s(P_kJ.g, C_kJ.g, F_kJ.g, k = 10)")

# Set the resolution of the surface (note that a higher value produces a nicer figure, but is larger in size and is slower to run)
surface.resolution<-501

# colour pallette
pall<-fields.colors(30)

# Color to show line for self-selected diet
line.col<-"black"

# How many values to round surface
round.surf<-3

# This specifies the color scheme for surface - it is actually a function that returns a function
rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")

# How many different colours should we use on the plot
no.cols<-256

# Get the colors to use from the pallette specified above
map<-rgb.palette(no.cols)

# How many levels should there be on the surface
nlev<-4

# The outcomes we are interested in modelling
outcomes<-c("brass.alpha", "brass.beta")

# surfaces in order of presentation for x, y, and median slices
surfaces.list<-list()
surfaces.list[[1]]<-c("P_kJ.g", "C_kJ.g", "F_kJ.g")
surfaces.list[[2]]<-c("P_kJ.g", "F_kJ.g", "C_kJ.g")
surfaces.list[[3]]<-c("C_kJ.g", "F_kJ.g", "P_kJ.g")

# Ratios of y to x on each panel
self.selected.diets<-c(c_s/p_s, f_s/p_s, f_s/c_s)

# Labels for each panel
labels<-list()
labels[[1]]<-c("Protein (kJ/g)", "Carbohydrate (kJ/g)", "Fat (kJ/g)")
labels[[2]]<-c("Protein (kJ/g)", "Fat (kJ/g)", "Carbohydrate (kJ/g)")
labels[[3]]<-c("Carbohydrate (kJ/g)", "Fat (kJ/g)", "Protein (kJ/g)")

titles<-list()
titles[[1]]<-expression(italic(alpha))
titles[[2]]<-expression(italic(beta))

# Panel labels
panels<-list()
panels[[1]]<-c("A", "", "")
panels[[2]]<-c("B", "", "")

# Lists to hold the models and fitted values on surfaces - they are useful for later
models.list<-list()
plots.list<-list()

# Do the surface for each item in the outcomes
for(k in 1:length(outcomes)){
	
	# Write some soace to an empty file
	if(k == 1){write.table("", file=res.file, sep=",", row.names=F, col.names=F)}
	
	# Find the right outcome
	summary.stats$this.outcome<-summary.stats[,outcomes[k]]
	
	# Get the weights based on the sampling variance for each estimate
	weights.i<-paste("V.", outcomes[k], sep="")		
	W<-1/summary.stats[,weights.i]
	
	# Fit the model
	model<-gam(model.form, data=summary.stats, weights=W)
	print(summary(model))
	
	# Write the results to the output
	write.table("", file=res.file, sep=",", row.names=F, col.names=F, append=T)
	write.table(outcomes[k], file=res.file, sep=",", row.names=F, col.names=F, append=T)
	res.k<-round(summary(model)$s.table, 4)
	res.k<-cbind(row.names(res.k), res.k)
	colnames(res.k)[1]<-"Coef."
	write.table(res.k, file=res.file, sep=",", row.names=F, col.names=colnames(res.k), append=T)
	
	# Save the model and also a prediction for each of the specific diets
	models.list[[k]]<-model
	
	# List to save the three surfaces
	sub.surf.list<-list()
	
	# Generate the three panels for the slices
	for(i in 1:3){
	
		# Predicted values of intake for surface
		Pred.Values<-findConvex(summary.stats[,surfaces.list[[i]][1]], summary.stats[,surfaces.list[[i]][2]], c(surfaces.list[[i]][1], surfaces.list[[i]][2]), res=surface.resolution)
		Pred.Values$med<-median(summary.stats[,surfaces.list[[i]][3]])
		names(Pred.Values)[3]<-surfaces.list[[i]][3]
		
		# Get the prediction from the GAM 
		out<-predict(model, newdata = Pred.Values)
		
		# Round to nicer values
		surf<-matrix(out, nrow=sqrt(dim(Pred.Values)[1]))
		surf<-round(surf, round.surf)
		
		# Save the surfaces for the kth outcome
		sub.surf.list[[i]]<-surf
		
		# Remove surf to be safe
		rm(surf)
		
	}	
	
	# We need to know the minimum and maximum predicted values to make the scale of the plots sensible
	mn<-min(unlist(sub.surf.list), na.rm=TRUE)
	mx<-max(unlist(sub.surf.list), na.rm=TRUE)	
	
	# Now make the surfaces
	for(i in 1:3){
		
		# surface i
		surf<-sub.surf.list[[i]]
		locs<-(range(surf, na.rm=TRUE) - mn) / (mx-mn) * no.cols
			
		# Pretty comes up with nice values of x and y over which to fit
		px<-pretty(summary.stats[,surfaces.list[[i]][1]])
		py<-pretty(summary.stats[,surfaces.list[[i]][2]])
		
		# Uses px and py to generate the x and y axes
		x.new<-seq(min(px), max(px), len=surface.resolution)
		y.new<-seq(min(py), max(py), len=surface.resolution)
				
		# Actually plots the surface using all of the above info above
		image(x.new, y.new, surf, col=map[locs[1]:locs[2]], xlab="", ylab="", axes=FALSE, main="")
		# Adds some axes
		axis(1, cex.axis = 1.5)
		axis(2, cex.axis = 1.5)
		# Adds a contour over the top (can add a title using main)
		contour(x.new, y.new, surf, add=TRUE, levels=pretty(range(mn,mx), nlev), labcex=1)

		# Add a title to the middle panel
		if(i == 2){mtext(titles[[k]], line=4, cex=5, font=2)}
		
		# Add the labels	
		mtext(labels[[i]][1], side=1, cex=1.25, line=3)
		mtext(labels[[i]][2], side=2, cex=1.25, line=3)
		mtext(paste(median(summary.stats[,surfaces.list[[i]][3]]), labels[[i]][3], sep=" "), side=1, cex=1, line=5)
		mtext(panels[[k]][i], side=2, cex=3, at=16, las=2)
		
		# Remove surf to be safe
		rm(surf)
	}
	
	# Save the three surfaces
	plots.list[[k]]<-sub.surf.list
	
	# Remove a bunch of paramters so they do not affect the analysis of the next outcome
	rm(model)
	rm(sub.surf.list)
	rm(Pred.Values)
}

# Remove the "this.outcome" column
summary.stats<-summary.stats[,-which(names(summary.stats) == "this.outcome")]

# Close the plotting file
dev.off()


# Figure 4 - converting the predicted brass coefficients in figure 3 in to lifetables, and estimates of e3 (note for analytical purposes this is actually e0) and S

# Note this first step i slow as it requires a calculation for each data point on the surface (high resolution surfaces take longer to convert)

# # We need to make sure we have the standard for estimation 
lsx<-brass_coefficients(ages_diet = raw.data$ageatdeath_w, ages_standard = raw.data$ageatdeath_w, return_standard = T)$lsx

# Open a loop to run 3 times - once for each slice/surface - note I am using paralellisation here on a multicore machine
registerDoParallel(3)
conversion<-foreach(i = 1:3) %dopar% {
	
	# Grab an already formed plot to use a template and make e0 and SD matrices in to which calculated values will be added
	e0<-plots.list[[1]][[i]]
	SD<-e0
	
	# We only need to consider those values on the surface that are not NAs (make things a bit quicker)
	calculate<-which(is.na(plots.list[[1]][[i]]) == F)
	
	# Grab the alphas and betas on the ith slice
	alphas<-plots.list[[1]][[i]][calculate]
	betas<-plots.list[[2]][[i]][calculate]
	
	# Loop to run through each value
	for(j in 1:length(calculate)){
		
		# Convert the jth alpha and beta value to e0 and standard deviation, and add in to the corresponding location on the surface
		e0[calculate[j]]<-convert_brass(alpha = alphas[j], beta = betas[j], lsx = lsx, age1 = 3, type = "e0")
		SD[calculate[j]]<-convert_brass(alpha = alphas[j], beta = betas[j], lsx = lsx, age1 = 3, type = "S")

	}
	
	# Save e0 and SD
	output<-list(e0, SD)
	return(output)
	
}

# The life expecs are the first object for each slice, and the SDs the second
e0s<-list(conversion[[1]][[1]], conversion[[2]][[1]], conversion[[3]][[1]])
SDs<-list(conversion[[1]][[2]], conversion[[2]][[2]], conversion[[3]][[2]])
rm(conversion)

# Save the surfaces - they take so long to generate you dont want to loose them - I recommend commenting out this section after a successful run 
save(e0s, file="Life_expectancy.Rdata")
save(SDs, file="SD_Life_expectancy.Rdata")

# Open the file for plotting
pdf("Fig4.pdf", height=13, width=15)

# Set the plotting region
par(mar=c(7,5,11,1), mfrow=c(2, 3))	

# Reload the results, in case you are running this section independently
load("Life_expectancy.Rdata")
load("SD_Life_expectancy.Rdata")

# Specify the isocaloric vector and resolution to test across
kJ.iso<-12
vector.res<-100
iso.vec<-cbind(seq(0, kJ.iso, length = vector.res), seq(kJ.iso, 0, length = vector.res))

# The outcomes we will plot on surfaces
outcomes<-list()
outcomes[[1]]<-e0s
outcomes[[2]]<-SDs

# Some nice titles to add
titles<-list()
titles[[1]]<-expression(italic(e)[3])
titles[[2]]<-expression(italic(s))

# Open the loop to makes the plots for each outcome
for(k in 1:length(outcomes)){
	
	# We need to know the minimum and maximum predicted values to make the scale of the plots sensible
	rng<-range(unlist(outcomes[[k]]), na.rm = TRUE)
	mn<-rng[1]
	mx<-rng[2]
	
	# List to hold the results across the isocaloric vector and the lifetables across the vectors
	iso.vec.list<-list()
	
	# Generate the three panels for the slices
	for(i in 1:3){
		
		# Get the ith surface for the kth outcome and round nicely 			
		out<-outcomes[[k]][[i]]
		surf<-matrix(out, nrow = surface.resolution)
		surf<-round(surf, round.surf)
		
		# Scale the colors
		locs<-(range(out, na.rm=T) - mn) / (mx-mn) * no.cols
		
		# Pretty comes up with nice values of x and y over which to fit
		px<-pretty(summary.stats[,surfaces.list[[i]][1]])
		py<-pretty(summary.stats[,surfaces.list[[i]][2]])
		
		# Uses px and py to generate the x and y axes
		x.new<-seq(min(px), max(px), len=surface.resolution)
		y.new<-seq(min(py), max(py), len=surface.resolution)
				
		# Actually plots the surface using all of the above info above
		image(x.new, y.new, surf, col=map[locs[1]:locs[2]], xlab="", ylab="", axes=FALSE, main="")
		# Adds some axes
		axis(1, cex.axis = 1.5)
		axis(2, cex.axis = 1.5)
		# Adds a contour over the top (can add a title using main)
		contour(x.new, y.new, surf, add=TRUE, nlevels=nlev, labcex=1)

		# Add a title to the middle panel
		if(i == 2){mtext(titles[[k]], line=4, cex=5, font=2)}
		
		# Add the labels	
		mtext(labels[[i]][1], side=1, cex=1.25, line=3)
		mtext(labels[[i]][2], side=2, cex=1.25, line=3)
		mtext(paste(median(summary.stats[,surfaces.list[[i]][3]]), labels[[i]][3], sep=" "), side=1, cex=1, line=5)
		mtext(panels[[k]][i], side=2, cex=3, at=16, las=2)
		abline(a = 0, b = self.selected.diets[i], lwd=3, col=line.col)
		
		# find the isocaloric vector across the polygon on the surface - we dont want ratios that fall outside the area of the visualised surface
		hull<-cbind(summary.stats[,surfaces.list[[i]][1]], summary.stats[,surfaces.list[[i]][2]])[chull(x = summary.stats[,surfaces.list[[i]][1]], y = summary.stats[,surfaces.list[[i]][2]]),]	
		search.space<-point.in.polygon(point.x = iso.vec[,1], point.y = iso.vec[,2], pol.x = hull[,1], pol.y = hull[,2])
		
		# Add in to the list as a dataframe
		iso.vec.list[[i]]<-as.data.frame(iso.vec[which(search.space == 1),])
		iso.vec.list[[i]]$nut3<-median(summary.stats[,surfaces.list[[i]][3]])
		names(iso.vec.list[[i]])<-surfaces.list[[i]]
		
		# Add the points across the surface
		points(iso.vec.list[[i]][,1], iso.vec.list[[i]][,2], pch=16, col="purple", cex=0.9)
		
	}
}

# Close the plotting file
dev.off()


# Figures 5 through 7 - LTRE analysis. NOTE this step can be a bit slow as calculating the sensitivities is a bit slow

# We need to make sure we have the standard for estimation 
lsx<-brass_coefficients(ages_diet = raw.data$ageatdeath_w, ages_standard = raw.data$ageatdeath_w, return_standard = T)$lsx

# Labels to go on the axes and the self selected ratios to plot
xlabs<-c("Protein / Carbohydrate", "Protein / Fat", "Carbohydrate / Fat")
self_select<-c(p_s / c_s, p_s / f_s, c_s / f_s)

# Some details about the size of plotting features
line.width=2.5
labels.size<-1.5

# Coloration and breaks for the legends on the level plots
breaks<-c(-1,-0.1,-0.01,-0.001,-0.0001,0.0001,0.001,0.01,0.1,1)
my.col<-cold.hot.colors(length(breaks)-1)
color.labels<-as.character(round(breaks, 4))
color.labels[c(3, 5, 6, 8)]<-""		


# Open the loop to do analyses for each vector
for(i in 1:length(iso.vec.list)){
			
	# Get predictions for alpha and beta across the iso vec
	brass.alpha<-predict(models.list[[1]], newdata = iso.vec.list[[i]])
	brass.beta<-predict(models.list[[2]], newdata = iso.vec.list[[i]])
	
	# List to hold the lifetables across the vector
	lifetables.list<-list()
	
	# Calculate the lifetable at each point along the isocaloric vector
	for(j in 1:dim(iso.vec.list[[i]])[1]){
		lifetables.list[[j]]<-convert_brass(alpha = brass.alpha[j], beta = brass.beta[j], lsx = lsx, type = "e0", age1 = 3, return_lt = T)
	}			

	# Get x ages at death in each category
	x<-lifetables.list[[1]]$age + 0.5
	x<-c(x, max(x)+0.5)
	
	# Convert the q values in liftables list in to a matrix
	qs<-sapply(lifetables.list, "[[", 3)
	rownames(qs)<-(x - 0.5)[-length(x)]
	colnames(qs)<-round(iso.vec.list[[i]][,1] / iso.vec.list[[i]][,2], 3)
	
	# Get the sensitivity of e0 and SD to qx (note the first object is the predicted e0/s at each point)
	e0_sens<-apply(qs, 2, return_sensitivity, x=x, type="e0")
	sd_sens<-apply(qs, 2, return_sensitivity, x=x, type="S")

	# List to hold the plotting objects
	plots.list<-list()
	
	# Now work out the change in ASM over the vector
	q<-qs#[-1,]
	
	# The ages and also the ratios we have here
	ages<-as.numeric(rownames(q))
	ratios<-as.numeric(colnames(q))	
	
	# First work out the change in each e3 and s over the vector
	change<-sapply(e0_sens, "[[", 1)
	delta_e0<-change - change[1]
	change<-sapply(sd_sens, "[[", 1)
	delta_s<-change - change[1]
	change.data<-as.data.frame(cbind(as.numeric(names(change)), delta_e0, delta_s))
		
	# Now make the plot of change over the vector
	plots.list[[1]]<-print(xyplot(delta_e0 ~ V1, data=change.data, 
	xlim=range(change.data$V1), ylim=c(-25, 25), type="l", xlab="", ylab = list(label=expression(Delta), cex=labels.size), lwd=line.width, col="black",
	scales=list(x=list(at=NULL)),
	key = list(space = "top", lines = list(col=c("black", "red"), lwd=line.width), cex=1.25, text = list(c(expression(italic(e[3])), expression(italic(s))))),
	par.settings = list(layout.heights=list(xlab.key.padding = -3.1)),
	panel=function(...){
		panel.xyplot(...)
		panel.abline(h=0, lty=2, lwd=line.width)
		panel.abline(v=self_select[i], col="grey", lwd=line.width)
		panel.lines(change.data$V1, change.data$delta_s, col="red", lwd=line.width)	
	}))
	
	# Note these labels are useful for all plots
	at.age<-seq(min(ages), max(ages), 10)
	at.ratio<-round(seq(min(ratios), max(ratios), length = 5), 1)
	
	# Data for plotting change in ASM
	plot.delta.q<-t(q - q[, 1])	
	
	# Transform in to dataframe for levelplot
	plot.data<-data.frame(ratios, plot.delta.q)
	names(plot.data)[1]<-c("ratios")
	plot.data<-reshape(plot.data, idvar = "ratios", direction="long", varying = list(2:dim(plot.data)[2]), v.names="z")
	plot.data$age<-plot.data$time + min(ages) - 1
	plot.data$type<-"delta q[x]"
	
	# Make the level plot for delta q			
	plots.list[[2]]<-print(levelplot(z ~ ratios * age | type, data = plot.data, at=breaks, col.regions=my.col, 
	ylab=list(label = "Age", cex=labels.size), xlab=list(label=xlabs[i], cex=labels.size), 
	colorkey=FALSE,
	scales=list(x=list(at=at.ratio, labels = at.ratio), y=list(at=at.age, labels=at.age)),
	strip=strip.custom(factor.levels = c(expression(Delta~italic(q[x]))), par.strip.text=list(cex = 1.5)),
	par.settings = list(layout.heights=list(strip=1.5, xlab.key.padding=3))))	
			
	# Get the differential for mortality over the vector	
	dq_dv<-t(apply(q, 1, gradient))
				
	# Do the LTRE for ith vectors life expectancy
	sens<-sapply(e0_sens, "[[", 2)
	#sens<-sens[-1,]
	C<-(sens * dq_dv)
	int<-apply(C, 1, cumsum)
	o<-t(int) * (1 / e0_sens[[1]]$index)
	o<-o[-dim(o)[1],]
	plot.o<-t(o)
	
	# Transform in to dataframe for levelplot
	plot.data<-data.frame(ratios, plot.o)
	names(plot.data)[1]<-c("ratios")
	plot.data<-reshape(plot.data, idvar = "ratios", direction="long", varying = list(2:dim(plot.data)[2]), v.names="z")
	plot.data$age<-plot.data$time + min(ages) - 1
	plot.data$type<-"e[3]"
	
	# Do the LTRE for ith vectors SD in age at death 
	sens<-sapply(sd_sens, "[[", 2)
	#sens<-sens[-1,]	
	C<-(sens * dq_dv)
	int<-apply(C, 1, cumsum)
	o<-t(int) * (1 / sd_sens[[1]]$index)
	o<-o[-dim(o)[1],]
	plot.o<-t(o)
	
	# Transform in to dataframe for levelplot
	plot.data2<-data.frame(ratios, plot.o)
	names(plot.data2)[1]<-c("ratios")
	plot.data2<-reshape(plot.data2, idvar = "ratios", direction="long", varying = list(2:dim(plot.data2)[2]), v.names="z")
	plot.data2$age<-plot.data2$time + min(ages) - 1
	plot.data2$type<-"s"	
	plot.data<-rbind(plot.data, plot.data2)
	
	# Make the plot	
	plots.list[[3]]<-print(levelplot(z ~ ratios * age | type, data = plot.data, at=breaks, col.regions=my.col, layout=c(2,1),
	ylab=list(label = "", cex=labels.size), xlab=list(label=xlabs[i], cex=labels.size),
	colorkey=list(col=my.col, at=seq(-1, 10, length=length(breaks)), space = "top", labels=color.labels),
	scales=list(x=list(at=at.ratio, labels = at.ratio), y=list(at=at.age, labels=at.age)),
	strip=strip.custom(factor.levels = c(expression(contributions~to~italic(Delta~e)[3]), expression(contributions~to~italic(Delta~s))), par.strip.text=list(cex = 1.2)),
	par.settings = list(layout.heights=list(strip=1.5, xlab.key.padding=3, top.padding = 0))))
	
	# A few bits for plotting the key
	x<-c(0.1, 0.3)
	y1<-c(0.75, 0.75)
	y2<-c(0.5, 0.5)
	y3<-c(0.25, 0.25)
	text.size<-1.5
	
	# Add the key
	plots.list[[4]]<-print(xyplot(y1 ~ x, type="l", col="black", lwd=line.width, ylim=c(0, 1), xlim=c(0, 1), ylab="", xlab="", bty="n",
	scales=list(x=list(at=NULL), y=list(at=NULL)),
	par.settings = list(layout.heights=list(top.padding = 0.5, xlab.key.padding = -5), axis.line = list(col=0)),
	panel=function(...){
		panel.xyplot(...)
		panel.lines(x, y3, col="grey", lwd=line.width)
		panel.lines(x, y2, col="red", lwd=line.width)
		panel.text(x[2]+0.05, y1[1], expression(italic(e[3])), cex=text.size, pos=4)
		panel.text(x[2]+0.05, y2[1], expression(italic(s)), cex=text.size, pos=4)
		panel.text(x[2]+0.05, y3[1], "self selected ratio", cex=text.size, pos=4)	
	}))
	
	# Open the plotting file	
	pdf(paste("Fig", i+4, ".pdf", sep=""), height=10, width=10)
	
	# Lay it all out nicely
	lattice.options(layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=-4)))
	
	hlay<-rbind(c(1,4),
				c(1,3),
				c(2,3),
				c(2,3),
				c(2,3),
			    c(2,3),
			    c(2,3),
			    c(2,3))
	
	# Print it out using grid arrange 
	grid.arrange(plots.list[[1]], plots.list[[2]], plots.list[[3]], plots.list[[4]], layout_matrix=hlay)
	
	# Close the plotting file
	dev.off()

}

