


rm(list=ls())

library(lmodel2)
library(PerformanceAnalytics)
library(nlme)
library(factoextra)
library(missMDA)
library(RColorBrewer)
library(scico)
library(plyr)
library(tagcloud)
library(smatr)
library(ape)
library(phytools)
library(pracma)
library(caper)
library(ggtree, quietly = T)
library(ggnewscale)
library(ggtreeExtra, quietly = T)
library(lme4)
library(lmerTest)



local_path = dirname(rstudioapi::getSourceEditorContext()$path)



############################# [1] formatting #############################

flight_dataset <- read.csv(paste0(local_path,'/hoverflies_species_averaged_flight_dataset.csv'))

#### build a species averaged morphological dataset matching the flight data to have the full dataset
morphol <- read.csv(paste0(local_path,'/hoverflies_weight_and_wing_data_ind_level.csv'), sep=',')
morphol <- morphol[which(morphol$filmed=='yes'),] ## select only the filmed specimen at this point
## compute mean parameter value per species
mean_morphol = as.data.frame(matrix(NA, length(levels(factor(morphol$genus_species))), length(colnames(morphol)[7:17])+1))
colnames(mean_morphol) = c('genus_species', colnames(morphol)[7:17])
mean_morphol$genus_species = as.factor(levels(factor(morphol$genus_species)))
for (i in 1:(length(colnames(mean_morphol))-1) ){
  mean_morphol[,(1+i)] = tapply(morphol[,which(colnames(morphol)==colnames(mean_morphol)[1+i])], morphol$genus_species, mean)
}
## make the row order to match between the flight and morpho dataset:
### re-order the rows in the flight data so as to match the morphol data
flight_dataset_ordered = flight_dataset
rownames(flight_dataset_ordered) = as.character(mean_morphol$genus_species)
for (i in 1:nrow(flight_dataset_ordered)){flight_dataset_ordered[i,] = flight_dataset[which(flight_dataset$genus_species==mean_morphol$genus_species[i]),]}

full_dataset = cbind(flight_dataset_ordered,mean_morphol[2:ncol(mean_morphol)])
rm(flight_dataset_ordered, flight_dataset, mean_morphol, morphol)

# sapply(full_dataset, class) # useful to check the class of all column  



## set character variables into factors
full_dataset$genus_species = as.factor(full_dataset$genus_species)

## add a variable angular_velo * S2
full_dataset$angvelo_time_S2 = (full_dataset$abs_wing_speed)^2*full_dataset$secondMoment

## locate numerical column in "full_dataset"
full_dataset_num_cols <- unlist(lapply(full_dataset, is.numeric)) # identify numeric columns

## compute wing loading
full_dataset$wing_loading <- full_dataset$fresh_weight / (full_dataset$wing_area_cm2*2)

## compute advance ratio (J = V / R*ang_velo, Ellington 1984 kinematics) with V the forward speed and R the wingspan
full_dataset$advance_ratio <- full_dataset$mean_body_speed / (full_dataset$length_cm * full_dataset$peak_ang_velo_fs)




##### *** double check S2 and S2* calculation

## the dimensionlss second-moment-of-area, S2* is define as
## S2* = S2 / ( S * b^2 )
## but in Ty Hedrick MATLAB code, it's computed as
## S2* = sqrt(S2) / (S * b^2)
## which is actually the r2*

## we compute our-self the S2*
full_dataset$nd_2ndMoment_custom = full_dataset$secondMoment / (full_dataset$wing_area_cm2 * (full_dataset$length_cm)^2)
# plot(full_dataset$nd_2ndMoment_custom , full_dataset$non_dimensional_2ndMoment , pch=16, cex = 2)


## --> WE SHALL USE ( WingImageProcessor S2* values )^2 IN OUR ANALYSIS !!
##     because the actual definition for S2* (S2* = S2 / ( S * b^2 )) correspond to (Ty Hedrick value)^2
full_dataset$non_dimensional_2ndMoment = (full_dataset$non_dimensional_2ndMoment)^2

### save "full_dataset"
# write.csv(full_dataset, paste0(local_path, '/full_dataset_species_level.csv'), row.names = F)





############ ___ set color code ############ 

##### color code matching weight
weight_color <- scico(25, palette = 'bamako')[c(1,3,6,8,10,12,13,15,17,19)] #choose gradient among: scico_palette_show()
weight_color = weight_color[length(weight_color):1] # reverse the gradient
# plot(rep(1,length(weight_color)), pch=15, cex=15,col=weight_color, bty='n', axes=F)
# set a gradient
variable_col_gradient = full_dataset$fresh_weight
gradient_now <- smoothPalette(variable_col_gradient, weight_color)
t <- colorRampPalette(c(weight_color[1], weight_color[length(weight_color)]))
seq_now = seq(min(variable_col_gradient), max(variable_col_gradient))

##### colors for species
levels(full_dataset$genus_species)

colors_species = c('#D1BA15', # Episyrphus_viridaureus  ## darker (greener) color are heavier species here
                   '#2B3A21', # Eristalis_tenax
                   '#78C155', # Eupeodes_nielseni
                   '#A63D40', # Melanostoma_mellinum
                   '#7BCEEF', # Platycheirus_clypeatus
                   '#674ACE', # Sphaerophoria_scripta
                   '#287C34', # Syrphus_nigrilinearus
                   '#F07E15') # Tropidia_scita

colo_species <- colors_species[match(full_dataset$genus_species, levels(full_dataset$genus_species))]

# plot legend only
plot(0, col='white',axes=F, xlab='', ylab='');legend('center', pch=16, pt.cex=2.2, c(levels(full_dataset$genus_species)), col=colors_species, bty='n', cex=0.9)




### plot variation in body mass among filmed species
body_mass_data <- as.data.frame(matrix(NA, nrow(full_dataset), 4))
colnames(body_mass_data) = c('species', 'weight', 'weight_SE','weight_SD')
body_mass_data$species = full_dataset$genus_species
body_mass_data$weight = full_dataset$fresh_weight
se <- function(x) sd(x) / sqrt(length(x))
### get standard errors:
weight_data <- read.csv(paste0(local_path,'/hoverflies_weight_and_wing_data_ind_level.csv'), h=T) # species order is different in those 
weight_data <- weight_data[which(weight_data$origin=='field'),]
weight_data$genus_species = as.factor(weight_data$genus_species)
for (i in 1:length(levels(full_dataset$genus_species))){
  body_mass_data$weight_SE[i] = se(weight_data$fresh_weight[which(weight_data$genus_species == body_mass_data$species[i])])
  body_mass_data$weight_SD[i] = sd(weight_data$fresh_weight[which(weight_data$genus_species == body_mass_data$species[i])])
}
## order following phylogeny:
body_mass_data$species <- factor(body_mass_data$species, levels = c('Tropidia_scita','Eristalis_tenax','Melanostoma_mellinum',
                                                     'Platycheirus_clypeatus','Episyrphus_viridaureus','Sphaerophoria_scripta',
                                                     'Syrphus_nigrilinearus','Eupeodes_nielseni'))
colors_species_tip_order = c('#F07E15', # Tropidia_scita
                   '#2B3A21', # Eristalis_tenax
                   '#A63D40', # Melanostoma_mellinum
                   '#7BCEEF', # Platycheirus_clypeatus
                   '#D1BA15', # Episyrphus_viridaureus
                   '#674ACE', # Sphaerophoria_scripta
                   '#287C34', # Syrphus_nigrilinearus
                   '#78C155') # Eupeodes_nielseni
ggplot(body_mass_data, aes(x=species, y=weight, fill=species)) + scale_fill_manual(values=colors_species_tip_order) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=weight-weight_SE, ymax=weight+weight_SE), width=.2) +
  theme_classic() + theme(legend.position = "none")

## define function to plot detailed log scale while keeping real value:
log10Tck <- function(side, type){
  lim <- switch(side, 
                x = par('usr')[1:2],
                y = par('usr')[3:4],
                stop("side argument must be 'x' or 'y'"))
  at <- floor(lim[1]) : ceil(lim[2])
  return(switch(type, 
                minor = outer(1:9, 10^(min(at):max(at))),
                major = 10^at,
                stop("type argument must be 'major' or 'minor'")
  ))
} 


############ ___ load phylogenies ############ 

### load the hoverfly phylogeny (from Wong et al. 2023. The phylogeny and evolutionary ecology of hoverflies)
### the tree was pruned and formatted in the code "format_hoverfly_phylogeny.R"
### tree matching the large morphological dataset (8 filmed species + 20 Leiden museum species:
load(paste0(local_path, '/phylogeny from Daniel Wong et al. 2023/hoverfly_tree_(morphology_dataset).RData'))
plot(phy, cex=0.6, edge.width = 0.5)
### tree matching the small flight dataset:
load(paste0(local_path, '/phylogeny from Daniel Wong et al. 2023/hoverfly_tree_(flight_dataset).RData'))
# plot(phy_filmed_only, cex=0.6, edge.width = 0.5)



############ [2] body mass vs wb kinematics [species lvl] ############

## variation in body vs kinematic can only be analyze at the species level because
## we don't know which weight data match which kinematic data/

# chart.Correlation(cbind(full_dataset$fresh_weight,full_dataset[,c(2:18)]), method='pearson', histogram=T, pch=16)

# plot(full_dataset$abs_wing_speed ~ full_dataset$fresh_weight, pch=16, bty='n', col=colo_species, cex=3, xlab='weight [mg]', ylab='wing speed [rad/sec]')
# plot(full_dataset$WB_freq ~ full_dataset$fresh_weight, pch=16, bty='n', col=colo_species, cex=3, xlab='weight [mg]', ylab='wing speed [rad/sec]')

plot(full_dataset$angvelo_time_S2 ~ full_dataset$fresh_weight, pch=16, bty='n', col=colo_species, cex=3, xlab='weight [mg]', ylab='S2 [cm^4] * (wing speed)^2 [rad/sec]')
plot(log10(full_dataset$angvelo_time_S2) ~ log10(full_dataset$fresh_weight), pch=16, bty='n', col=colo_species, cex=3, xlab='log10(weight)', ylab='log10(S2 [cm^4] * (wing speed)^2 [rad/sec])')
abline(lm(log10(full_dataset$angvelo_time_S2) ~ log10(full_dataset$fresh_weight)), col='black', lwd=2)

layout(matrix(1:2,ncol=2), width = c(2,1), height = c(1,1))
plot(full_dataset$angvelo_time_S2 ~ full_dataset$fresh_weight, pch=16, bty='n', col=gradient_now, cex=3, xlab='weight [mg]', ylab='S2 [cm^4] * (wing speed)^2 [rad/sec]')
legend_image <- as.raster(matrix(t(length(seq_now)), ncol=1))
plot(c(0,2), c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'weight [mg]', cex.main = 1)
text(x=1.4, y = seq(0,1,l=5), labels = seq(5, 100, l=5))
rasterImage(legend_image[length(seq_now):1], 0.3, 0, 0.7, 1) # inverse color direction here with "[length(seq_now):1]"
layout(matrix(1,1))

plot(log10(full_dataset$stroke_amplitude) ~ log10(full_dataset$fresh_weight), pch=16, bty='n', col=colo_species, cex=3, xlab='log10(weight)', ylab='log10(stroke amplitude)')
plot(log10(full_dataset$rotation_amplitude) ~ log10(full_dataset$fresh_weight), pch=16, bty='n', col=colo_species, cex=3, xlab='log10(weight)', ylab='log10(rotation amplitude)')
plot(log10(full_dataset$deviation_amplitude) ~ log10(full_dataset$fresh_weight), pch=16, bty='n', col=colo_species, cex=3, xlab='log10(weight)', ylab='log10(deviation amplitude)')

plot(log10(full_dataset$peak_ang_accel) ~ log10(full_dataset$fresh_weight), pch=16, bty='n', col=colo_species, cex=3, xlab='log10(weight)', ylab='log10(peak ang. acceleration)')






############ ___ WB freq vs mass ############

WB_freq_data <- as.data.frame(matrix(NA, nrow(full_dataset), 5))
colnames(WB_freq_data) = c('species', 'WB_freq', 'WB_freq_SE', 'weight', 'weight_SE')
WB_freq_data$species = full_dataset$genus_species
WB_freq_data$WB_freq = full_dataset$WB_freq
WB_freq_data$weight = full_dataset$fresh_weight
se <- function(x) sd(x) / sqrt(length(x))

### get standard errors for wb frequency:
wb_parameters <- read.csv(paste0(local_path,'/wb_parameters.csv'), h=T)
wb_parameters$wing_speed_ftimesA = wb_parameters$WB_freq * wb_parameters$stroke_amplitude
weight_data <- read.csv(paste0(local_path,'/hoverflies_weight_and_wing_data_ind_level.csv'), h=T) # species order is different in those 
weight_data <- weight_data[which(weight_data$filmed=='yes'),] ## select only the filmed specimen at this point
weight_data$genus_species = as.factor(weight_data$genus_species)
for (i in 1:length(levels(full_dataset$genus_species))){
  WB_freq_data$WB_freq_SE[i] = se(wb_parameters$WB_freq[which(wb_parameters$genus_species == WB_freq_data$species[i])])
  WB_freq_data$weight_SE[i] = se(weight_data$fresh_weight[which(weight_data$genus_species == WB_freq_data$species[i])])
}

### make graph with both vertical and horizontal error bar
# https://ggplot2.tidyverse.org/reference/geom_errorbarh.html

ggplot(WB_freq_data, aes(x = weight, y = WB_freq, label = species)) + #geom_text(hjust=0, vjust=0) + 
  geom_point(size=7, color=colors_species) + labs(y = 'wingbeat frequency [Hz]', x='body mass [mg]') +
  geom_errorbar(aes(ymax = WB_freq + WB_freq_SE, ymin = WB_freq - WB_freq_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 5) + ylim(0,250) +
  theme_classic() + theme(legend.position = "none")
lmodel2(log10(WB_freq_data$WB_freq) ~ log10(WB_freq_data$weight))


### re-order the rows in 'WB_freq_data' to match the order of the tips in the phylogeny 'phy_filmed_only'
WB_freq_data <- WB_freq_data[match(phy_filmed_only$tip.label, WB_freq_data$species), ]
### re-order the colors 
colors_species_2 = c('#78C155', # Eupeodes_nielseni  ## darker (greener) color are heavier species here
                   '#287C34', # Syrphus_nigrilinearus
                   '#674ACE', # Sphaerophoria_scripta
                   '#D1BA15', # Episyrphus_viridaureus
                   '#7BCEEF', # Platycheirus_clypeatus
                   '#A63D40', # Melanostoma_mellinum
                   '#2B3A21', # Eristalis_tenax
                   '#F07E15') # Tropidia_scita


## PGLS using nlme on LOG TRANSFORMED DATA
## (following https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(WB_freq_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log10(WB_freq_data$WB_freq); names(y_now) = phy_filmed_only$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)

## Obtain the adjusted R-squared (not given through gls() function so we compute it ourselves) 
residuals <- residuals(pglsModel) # extract residuals
RSS <- sum(residuals^2) # compute RSS (Residual Sum of Squares)
## Compute the Total Sum of Square (TSS)
fitted_values <- pglsModel$fitted # extract fitted values
response_values <- fitted_values + residuals # Reconstruct the response variable
TSS <- sum((response_values - mean(response_values))^2) # Compute Total Sum of Squares (TSS)
# Get number of observations (n) and predictors (k)
n <- length(residuals) # Get number of observations (n)
k <- length(coef(pglsModel)) - 1  # Get number of predictors (k)
## Compute Adjusted R-squared
R2_adj <- 1 - ((RSS / (n - k - 1)) / (TSS / (n - 1))) ; R2_adj ## sounds uncorrect
summary(lm(log10(WB_freq_data$WB_freq) ~ log10(WB_freq_data$weight)))

## Obtain IC 95% of the pgls regression slope
coef_est <- coef(summary(pglsModel))["x_now", "Value"] # extract coefficient estimates
se_est <- coef(summary(pglsModel))["x_now", "Std.Error"]  # extract standard errors
df <- summary(pglsModel)$dims$N - summary(pglsModel)$dims$p  # get residual degrees of freedom (df) (nb of observation - nb of estimated parameters)
## Compute t-value for 95% CI (alpha = 0.05, two-tailed)
t_value <- qt(0.975, df)  # 97.5th percentile for two-tailed test
## Compute confidence interval
CI_lower <- coef_est - t_value * se_est
CI_upper <- coef_est + t_value * se_est
## Display the Confidence Interval of the slope
c(CI_lower, CI_upper)
## OR SIMPLY USE conint()
confint(pglsModel)



#### plot in log scale
#### including the scaling hypotheses and (non-significant) fitted lines onto the graph

iso_intercept <- mean(log10(WB_freq_data$WB_freq)) - 0 * mean(log10(WB_freq_data$weight))
intercept_W <- mean(log10(WB_freq_data$WB_freq)) - (-1/6) * mean(log10(WB_freq_data$weight))

ggplot(WB_freq_data, aes(x = log10(weight), y = log10(WB_freq), label = species)) + 
  geom_point(size=7, color=colors_species_2) + 
  labs(y = 'log10(wingbeat frequency)', x='log10(body mass)') +
  geom_errorbar(aes(ymax = log10(WB_freq + WB_freq_SE), ymin = log10(WB_freq - WB_freq_SE)), width = 0.025) +
  geom_errorbarh(aes(xmax = log10(weight + weight_SE), xmin = log10(weight - weight_SE)), height = 0.005) + 
  theme_classic() + 
  theme(legend.position = "none") +
  ## isometry:
  geom_abline(intercept = iso_intercept, slope = (0), color = '#ba4c36', linetype = "dashed", lwd=1.5) + # isometry
  geom_abline(intercept = 2.3722573, slope = -0.0848472, color = 'black', linetype = "solid", lwd=1.5) + # pgls
  geom_abline(intercept = intercept_W, slope = (-1/6), color = '#A0A0A0', linetype = "dotted", lwd=1.5) # required for weight support if everything else scales isometrically:


# Plot in log-log scale while showing real values using Base R
plot(WB_freq_data$weight, WB_freq_data$WB_freq, 
     log='xy', pch=16, cex=3, bty='n', col=colors_species_2, 
     ylab='wingbeat frequency [Hz]', xlab='body mass [mg]', axes=F, ylim=c(80, 350))
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

## recompute the PGLS with the log 10 values:
x_log <- log10(WB_freq_data$weight)
y_log <- log10(WB_freq_data$WB_freq)
pglsModel <- gls(y_log ~ x_log, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)
abline(2.3722573, -0.0848472 , col='black', lwd=2) # PGLS
abline(2.3722573,(0), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(2.3722573, (-1/6), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically

## add error bars on the log scale plot showing real values
# Convert standard errors to log space (calculate upper and lower limits)
log_weight_SE_upper <- log10(WB_freq_data$weight + WB_freq_data$weight_SE)
log_weight_SE_lower <- log10(WB_freq_data$weight - WB_freq_data$weight_SE)
log_WB_freq_SE_upper <- log10(WB_freq_data$WB_freq + WB_freq_data$WB_freq_SE)
log_WB_freq_SE_lower <- log10(WB_freq_data$WB_freq - WB_freq_data$WB_freq_SE)
# Make sure lower limits are valid (avoid log of zero or negative values)
log_weight_SE_lower[WB_freq_data$weight - WB_freq_data$weight_SE <= 0] <- NA
log_WB_freq_SE_lower[WB_freq_data$WB_freq - WB_freq_data$WB_freq_SE <= 0] <- NA
# Add horizontal error bars (Body Mass SE)
arrows(x0 = 10^log_weight_SE_lower, x1 = 10^log_weight_SE_upper,
       y0 = WB_freq_data$WB_freq, y1 = WB_freq_data$WB_freq,
       angle=90, code=3, length=0.05, col="black")

# Add vertical error bars (Wingbeat Frequency SE)
arrows(x0 = WB_freq_data$weight, x1 = WB_freq_data$weight,
       y0 = 10^log_WB_freq_SE_lower, y1 = 10^log_WB_freq_SE_upper,
       angle=90, code=3, length=0.05, col="black")





##### Plot in linear scale (power function)

## instead of linear trends, we should be using power-law relationships for both the PGLS and the weight-support trends.
## In these cases, the relationship is not linear but follows a specific scaling rule :
## f = k*m^(−1/6)
## equivalent to f = k / m^(1/6)


#### find the k value
# Transform the data
x_now <- log(WB_freq_data$weight)
y_now <- log(WB_freq_data$WB_freq)
# Fit pgls model to the transformed data
model <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
# Extract the intercept (log(k)) and slope
intercept <- coef(model)[1]
slope_estimated <- coef(model)[2]
# Calculate k by taking the exponential of the intercept
k_estimated <- exp(intercept)

# Define the scaling function for weight-support (k * m^(-1/6)) (needed in ggplot below)
weight_support_function <- function(m) {
  k <- k_estimated  # Adjust this value to your specific constant 'k'
  return(k / m^(1/6))
}

#### Define the scaling function for PGLS (k * m^(slope)) with the estimated slope
pgls_function <- function(m) {
  k <- k_estimated # Adjust this value to your specific constant 'k'
  slope <- slope_estimated  # Adjust the slope based on your PGLS model
  return(k * m^slope)
}

iso_intercept <- mean(WB_freq_data$WB_freq) - 0 * mean(WB_freq_data$weight)

ggplot(WB_freq_data, aes(x = weight, y = WB_freq, label = species)) + 
  geom_point(size=7, color=colors_species_2) + 
  labs(y = 'wingbeat frequency [Hz]', x = 'body mass [mg]') +
  geom_errorbar(aes(ymax = WB_freq + WB_freq_SE, ymin = WB_freq - WB_freq_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 5) + 
  ylim(0, 250) +
  theme_classic() + 
  theme(legend.position = "none") +
  # Plotting isometry line
  geom_abline(intercept = iso_intercept, slope = (0), color = '#ba4c36', linetype = "dashed", lwd=1.5) + # isometry
  # Plotting the weight-support scaling as a power-law curve
  geom_function(fun = weight_support_function, color = '#A0A0A0', linetype = "dotted", lwd = 1.5) + 
  # Plotting the PGLS scaling as a power-law curve
  geom_function(fun = pgls_function, color = 'black', linetype = "solid", lwd = 1.5)





####### Phylogenetic Reduced Major Axis (phyloRMA):
x_now = log(WB_freq_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log(WB_freq_data$WB_freq); names(y_now) = phy_filmed_only$tip.label
rma_result <- phyl.RMA(x_now, y_now, phy_filmed_only, fixed=T) ; rma_result

### test is it deviate from isometry based on function-embedded test
phyl.RMA(x_now, y_now, phy_filmed_only, fixed=T, h0=0.1)







############ ___ wing speed vs mass ############

wing_speed_data <- as.data.frame(matrix(NA, nrow(full_dataset), 5))
colnames(wing_speed_data) = c('species', 'wing_speed', 'wing_speed_SE', 'weight', 'weight_SE')
wing_speed_data$species = full_dataset$genus_species
wing_speed_data$wing_speed = full_dataset$abs_wing_speed
wing_speed_data$weight = full_dataset$fresh_weight
for (i in 1:length(levels(full_dataset$genus_species))){
  wing_speed_data$wing_speed_SE[i] = se(wb_parameters$abs_wing_speed[which(wb_parameters$genus_species == wing_speed_data$species[i])])
  wing_speed_data$weight_SE[i] = se(weight_data$fresh_weight[which(weight_data$genus_species == wing_speed_data$species[i])])
}

ggplot(wing_speed_data, aes(x = weight, y = wing_speed, label = species)) + #geom_text(hjust=0, vjust=0) + 
  geom_point(size=7, color=colors_species) + labs(y = 'angular speed [rad/sec]', x='body mass [mg]') +
  geom_errorbar(aes(ymax = wing_speed + wing_speed_SE, ymin = wing_speed - wing_speed_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 1000) + ylim(0,80000) +
  theme_classic() + theme(legend.position = "none")
lmodel2(log10(wing_speed_data$wing_speed) ~ log10(wing_speed_data$weight))

## plot for Table summing up scaling expectation:
plot(log10(wing_speed_data$wing_speed) ~ log10(wing_speed_data$weight), pch=16, col='grey80', cex=1.5, bty='L', xlab='log(body mass)',ylab='log(wing angular speed)', asp=F)
abline(4.781954, -0.06868760  , col='black', lwd=4.5)  # observed
iso_intercept <- mean(log10(wing_speed_data$wing_speed)) - 0 * mean(log10(wing_speed_data$weight))
abline(a = iso_intercept, b = 0, col = '#ba4c36', lwd=4.5, lty='dashed') # isometric
intercept_W <- mean(log10(wing_speed_data$wing_speed)) - (-1/6) * mean(log10(wing_speed_data$weight))
abline(a = intercept_W, b = -1/6, col = '#0e9594', lwd=4.5, lty='solid') # weight support





############ ___ wing speed (f x A) vs mass ############


wing_speed_data <- as.data.frame(matrix(NA, nrow(full_dataset), 5))
colnames(wing_speed_data) = c('species', 'wing_speed_fA', 'wing_speed_SE', 'weight', 'weight_SE')
wing_speed_data$species = full_dataset$genus_species
wing_speed_data$wing_speed_fA = full_dataset$WB_freq * full_dataset$stroke_amplitude
wing_speed_data$weight = full_dataset$fresh_weight
for (i in 1:length(levels(full_dataset$genus_species))){
  wing_speed_data$wing_speed_SE[i] = se(wb_parameters$wing_speed_ftimesA[which(wb_parameters$genus_species == wing_speed_data$species[i])])
  wing_speed_data$weight_SE[i] = se(weight_data$fresh_weight[which(weight_data$genus_species == wing_speed_data$species[i])])
}

ggplot(wing_speed_data, aes(x = weight, y = wing_speed_fA, label = species)) + #geom_text(hjust=0, vjust=0) + 
  geom_point(size=7, color=colors_species) + labs(y = 'angular speed (f x A) [deg/sec]', x='body mass [mg]') +
  geom_errorbar(aes(ymax = wing_speed_fA + wing_speed_SE, ymin = wing_speed_fA - wing_speed_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 1000) + ylim(0,27000) +
  theme_classic() + theme(legend.position = "none")
lmodel2(log10(wing_speed_data$wing_speed_fA) ~ log10(wing_speed_data$weight))



### re-order the rows in 'wing_speed_data' to match the order of the tips in the phylogeny 'phy_filmed_only'
wing_speed_data <- wing_speed_data[match(phy_filmed_only$tip.label, wing_speed_data$species), ]


## PGLS using nlme on LOG TRANSFORMED DATA
## (following https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(wing_speed_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log10(wing_speed_data$wing_speed_fA); names(y_now) = phy_filmed_only$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)

## Obtain the adjusted R-squared (not given through gls() function so we compute it ourselves)
residuals <- residuals(pglsModel) # extract residuals
RSS <- sum(residuals^2) # compute RSS (Residual Sum of Squares)
## Compute the Total Sum of Square (TSS)
fitted_values <- pglsModel$fitted # extract fitted values
response_values <- fitted_values + residuals # Reconstruct the response variable
TSS <- sum((response_values - mean(response_values))^2) # Compute Total Sum of Squares (TSS)
# Get number of observations (n) and predictors (k)
n <- length(residuals) # Get number of observations (n)
k <- length(coef(pglsModel)) - 1  # Get number of predictors (k)
## Compute Adjusted R-squared
R2_adj <- 1 - ((RSS / (n - k - 1)) / (TSS / (n - 1))) ; R2_adj
summary(lm(log(wing_speed_data$wing_speed_fA) ~ log(WB_freq_data$weight)))

## Obtain IC 95% of the pgls regression slope 
coef_est <- coef(summary(pglsModel))["x_now", "Value"] # extract coefficient estimates
se_est <- coef(summary(pglsModel))["x_now", "Std.Error"]  # extract standard errors
df <- summary(pglsModel)$dims$N - summary(pglsModel)$dims$p  # get residual degrees of freedom (df) (nb of observation - nb of estimated parameters)
## Compute t-value for 95% CI (alpha = 0.05, two-tailed)
t_value <- qt(0.975, df)  # 97.5th percentile for two-tailed test
## Compute confidence interval
CI_lower <- coef_est - t_value * se_est
CI_upper <- coef_est + t_value * se_est
## Display the Confidence Interval of the slope
c(CI_lower, CI_upper)
## OR SIMPLY USE conint()
confint(pglsModel)



#### plot in log scale
#### including the scaling hypotheses and (non-significant) fitted lines onto the graph

iso_intercept <- mean(log10(wing_speed_data$wing_speed_fA)) - 0 * mean(log10(wing_speed_data$weight))
intercept_W <- mean(log10(wing_speed_data$wing_speed_fA)) - (-1/6) * mean(log10(wing_speed_data$weight))

ggplot(wing_speed_data, aes(x = log10(weight), y = log10(wing_speed_fA), label = species)) + 
  geom_point(size=7, color=colors_species_2) + 
  labs(y = 'log10(angular speed (f x A) )', x='log10(body mass)') +
  geom_errorbar(aes(ymax = log10(wing_speed_fA + wing_speed_SE), ymin = log10(wing_speed_fA - wing_speed_SE)), width = 0.03) +
  geom_errorbarh(aes(xmax = log10(weight + weight_SE), xmin = log10(weight - weight_SE)), height = 0.009) + 
  theme_classic() + 
  geom_abline(intercept = iso_intercept, slope = (0), color = '#ba4c36', linetype = "dashed", lwd=1.5) + # isometry
  geom_abline(intercept = 4.413195, slope = -0.101888, color = 'black', linetype = "solid", lwd=1.5) + # pgls
  geom_abline(intercept = intercept_W, slope = (-1/6), color = '#A0A0A0', linetype = "dotted", lwd=1.5)# required for weight support if everything else scales isometrically
 


# Plot in log-log scale while showing real values using Base R
plot(wing_speed_data$weight, wing_speed_data$wing_speed_fA, 
     log='xy', pch=16, cex=3, bty='n', col=colors_species_2, 
     ylab='angular speed [deg/sec]', xlab='body mass [mg]', axes=F, ylim=c(0.5e+04, 4e+04))
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

## recompute the PGLS with the log 10 values:
x_log <- log10(wing_speed_data$weight)
y_log <- log10(wing_speed_data$wing_speed_fA)
pglsModel <- gls(y_log ~ x_log, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)
abline(4.413195, -0.101888 , col='black', lwd=2) # PGLS
abline(4.413195,(0), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(4.413195, (-1/6), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically

## add error bars on the log scale plot showing real values
# Convert standard errors to log space (calculate upper and lower limits)
log_weight_SE_upper <- log10(wing_speed_data$weight + wing_speed_data$weight_SE)
log_weight_SE_lower <- log10(wing_speed_data$weight - wing_speed_data$weight_SE)
log_wing_speed_SE_upper <- log10(wing_speed_data$wing_speed_fA + wing_speed_data$wing_speed_SE)
log_wing_speed_SE_lower <- log10(wing_speed_data$wing_speed_fA - wing_speed_data$wing_speed_SE)
# Make sure lower limits are valid (avoid log of zero or negative values)
log_weight_SE_lower[wing_speed_data$weight - wing_speed_data$weight_SE <= 0] <- NA
log_wing_speed_SE_lower[wing_speed_data$wing_speed_fA - wing_speed_data$wing_speed_SE <= 0] <- NA
# Add horizontal error bars (Body Mass SE)
arrows(x0 = 10^log_weight_SE_lower, x1 = 10^log_weight_SE_upper,
       y0 = wing_speed_data$wing_speed_fA, y1 = wing_speed_data$wing_speed_fA,
       angle=90, code=3, length=0.05, col="black")

# Add vertical error bars
arrows(x0 = wing_speed_data$weight, x1 = wing_speed_data$weight,
       y0 = 10^log_wing_speed_SE_lower, y1 = 10^log_wing_speed_SE_upper,
       angle=90, code=3, length=0.05, col="black")



##### Plot in linear scale (power function)

## instead of linear trends, we should be using power-law relationships for both the PGLS and the weight-support trends.
## In these cases, the relationship is not linear but follows a specific scaling rule :
## f = k*m^(−1/6)
## equivalent to f = k / m^(1/6)


#### find the k value
# Transform the data
x_now <- log(wing_speed_data$weight)
y_now <- log(wing_speed_data$wing_speed_fA)
# Fit pgls model to the transformed data
model <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
# Extract the intercept (log(k)) and slope
intercept <- coef(model)[1]
slope_estimated <- coef(model)[2]
# Calculate k by taking the exponential of the intercept
k_estimated <- exp(intercept)

# Define the scaling function for weight-support (k * m^(-1/6)) (needed in ggplot below)
weight_support_function <- function(m) {
  k <- k_estimated  # Adjust this value to your specific constant 'k'
  return(k / m^(1/6))
}

#### Define the scaling function for PGLS (k * m^(slope)) with the estimated slope
pgls_function <- function(m) {
  k <- k_estimated # Adjust this value to your specific constant 'k'
  slope <- slope_estimated  # Adjust the slope based on your PGLS model
  return(k * m^slope)
}

iso_intercept <- mean(wing_speed_data$wing_speed_fA) - 0 * mean(wing_speed_data$weight)

ggplot(wing_speed_data, aes(x = weight, y = wing_speed_fA, label = species)) + 
  geom_point(size=7, color=colors_species_2) + 
  labs(y = 'angular speed (f x A) [deg/sec]', x = 'body mass [mg]') +
  geom_errorbar(aes(ymax = wing_speed_fA + wing_speed_SE, ymin = wing_speed_fA - wing_speed_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 700) + 
  ylim(0,27000) +
  theme_classic() + 
  theme(legend.position = "none") +
  # Plotting isometry line
  geom_abline(intercept = iso_intercept, slope = (0), color = '#ba4c36', linetype = "dashed", lwd=1.5) + # isometry
  # Plotting the weight-support scaling as a power-law curve
  geom_function(fun = weight_support_function, color = '#A0A0A0', linetype = "dotted", lwd = 1.5) + 
  # Plotting the PGLS scaling as a power-law curve
  geom_function(fun = pgls_function, color = 'black', linetype = "solid", lwd = 1.5)




####### Phylogenetic Reduced Major Axis (phyloRMA):
x_now = log(wing_speed_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log(wing_speed_data$wing_speed_fA); names(y_now) = phy_filmed_only$tip.label
rma_result <- phyl.RMA(x_now, y_now, phy_filmed_only, fixed=T) ; rma_result

### test is it deviate from isometry based on function-embedded test
phyl.RMA(x_now, y_now, phy_filmed_only, fixed=T, h0=0.1)









############ ___ stroke amplitude vs mass ############

stroke_amp_data <- as.data.frame(matrix(NA, nrow(full_dataset), 5))
colnames(stroke_amp_data) = c('species', 'stroke_amp', 'stroke_amp_SE', 'weight', 'weight_SE')
stroke_amp_data$species = full_dataset$genus_species
stroke_amp_data$stroke_amp = full_dataset$stroke_amplitude
stroke_amp_data$weight = full_dataset$fresh_weight
for (i in 1:length(levels(full_dataset$genus_species))){
  stroke_amp_data$stroke_amp_SE[i] = se(wb_parameters$stroke_amplitude[which(wb_parameters$genus_species == stroke_amp_data$species[i])])
  stroke_amp_data$weight_SE[i] = se(weight_data$fresh_weight[which(weight_data$genus_species == stroke_amp_data$species[i])])
}

ggplot(stroke_amp_data, aes(x = weight, y = stroke_amp, label = species)) + #geom_text(hjust=0, vjust=0) + 
  geom_point(size=7, color=colors_species) + labs(y = 'stroke amplitude [deg]', x='body mass [mg]') +
  geom_errorbar(aes(ymax = stroke_amp + stroke_amp_SE, ymin = stroke_amp - stroke_amp_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 1) + ylim(0,150) +
  theme_classic() + theme(legend.position = "none")
lmodel2(log10(stroke_amp_data$stroke_amp) ~ log10(stroke_amp_data$weight))


### re-order the rows in 'stroke_amp_data' to match the order of the tips in the phylogeny 'phy_filmed_only'
stroke_amp_data <- stroke_amp_data[match(phy_filmed_only$tip.label, stroke_amp_data$species), ]

## PGLS using nlme on LOG TRANSFORMED DATA
## (following https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(stroke_amp_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log10(stroke_amp_data$stroke_amp); names(y_now) = phy_filmed_only$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)

## Obtain the adjusted R-squared (not given through gls() function so we compute it ourselves)  
residuals <- residuals(pglsModel) # extract residuals
RSS <- sum(residuals^2) # compute RSS (Residual Sum of Squares)
## Compute the Total Sum of Square (TSS)
fitted_values <- pglsModel$fitted # extract fitted values
response_values <- fitted_values + residuals # Reconstruct the response variable
TSS <- sum((response_values - mean(response_values))^2) # Compute Total Sum of Squares (TSS)
# Get number of observations (n) and predictors (k)
n <- length(residuals) # Get number of observations (n)
k <- length(coef(pglsModel)) - 1  # Get number of predictors (k)
## Compute Adjusted R-squared
R2_adj <- 1 - ((RSS / (n - k - 1)) / (TSS / (n - 1))) ; R2_adj
summary(lm(log(stroke_amp_data$stroke_amp) ~ log(stroke_amp_data$weight)))

## Obtain IC 95% of the pgls regression slope 
coef_est <- coef(summary(pglsModel))["x_now", "Value"] # extract coefficient estimates
se_est <- coef(summary(pglsModel))["x_now", "Std.Error"]  # extract standard errors
df <- summary(pglsModel)$dims$N - summary(pglsModel)$dims$p  # get residual degrees of freedom (df) (nb of observation - nb of estimated parameters)
## Compute t-value for 95% CI (alpha = 0.05, two-tailed)
t_value <- qt(0.975, df)  # 97.5th percentile for two-tailed test
## Compute confidence interval
CI_lower <- coef_est - t_value * se_est
CI_upper <- coef_est + t_value * se_est
## Display the Confidence Interval of the slope
c(CI_lower, CI_upper)
## OR SIMPLY USE conint()
confint(pglsModel)


#### plot in log scale
#### including the scaling hypotheses and (non-significant) fitted lines onto the graph

iso_intercept <- mean(log10(stroke_amp_data$stroke_amp)) - 0 * mean(log10(stroke_amp_data$weight))
intercept_W <- mean(log10(stroke_amp_data$stroke_amp)) - (-1/6) * mean(log10(stroke_amp_data$weight))

ggplot(stroke_amp_data, aes(x = log10(weight), y = log10(stroke_amp), label = species)) + 
  geom_point(size=7, color=colors_species_2) + 
  labs(y = 'log10(stroke amplitude)', x='log10(body mass)') +
  geom_errorbar(aes(ymax = log10(stroke_amp + stroke_amp_SE), ymin = log10(stroke_amp - stroke_amp_SE)), width = 0.025) +
  geom_errorbarh(aes(xmax = log10(weight + weight_SE), xmin = log10(weight - weight_SE)), height = 0.005) + 
  theme_classic() + 
  theme(legend.position = "none") +
  ## isometry:
  geom_abline(intercept = iso_intercept, slope = (0), color = '#ba4c36', linetype = "dashed", lwd=1.5) + # isometry
  geom_abline(intercept = 2.0409379, slope = -0.0170404, color = 'black', linetype = "solid", lwd=1.5) + # pgls
  geom_abline(intercept = intercept_W, slope = (-1/6), color = '#A0A0A0', linetype = "dotted", lwd=1.5) # required for weight support if everything else scales isometrically




# Plot in log-log scale while showing real values using Base R
plot(stroke_amp_data$weight, stroke_amp_data$stroke_amp, 
     log='xy', pch=16, cex=3, bty='n', col=colors_species_2, 
     ylab='stroke amplitude [deg]', xlab='body mass [mg]', axes=F, ylim=c(60, 150))
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

## recompute the PGLS with the log 10 values:
x_log <- log10(stroke_amp_data$weight)
y_log <- log10(stroke_amp_data$stroke_amp)
pglsModel <- gls(y_log ~ x_log, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)
abline(2.0409379, -0.0170404 , col='black', lwd=2) # PGLS
abline(2.0409379,(0), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(2.0409379, (-1/6), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically

## add error bars on the log scale plot showing real values
# Convert standard errors to log space (calculate upper and lower limits)
log_weight_SE_upper <- log10(stroke_amp_data$weight + stroke_amp_data$weight_SE)
log_weight_SE_lower <- log10(stroke_amp_data$weight - stroke_amp_data$weight_SE)
log_stroke_amp_SE_upper <- log10(stroke_amp_data$stroke_amp + stroke_amp_data$stroke_amp_SE)
log_stroke_amp_SE_lower <- log10(stroke_amp_data$stroke_amp - stroke_amp_data$stroke_amp_SE)
# Make sure lower limits are valid (avoid log of zero or negative values)
log_weight_SE_lower[stroke_amp_data$weight - stroke_amp_data$weight_SE <= 0] <- NA
log_stroke_amp_SE_lower[stroke_amp_data$stroke_amp - stroke_amp_data$stroke_amp_SE <= 0] <- NA
# Add horizontal error bars (Body Mass SE)
arrows(x0 = 10^log_weight_SE_lower, x1 = 10^log_weight_SE_upper,
       y0 = stroke_amp_data$stroke_amp, y1 = stroke_amp_data$stroke_amp,
       angle=90, code=3, length=0.05, col="black")

# Add vertical error bars
arrows(x0 = stroke_amp_data$weight, x1 = stroke_amp_data$weight,
       y0 = 10^log_stroke_amp_SE_lower, y1 = 10^log_stroke_amp_SE_upper,
       angle=90, code=3, length=0.05, col="black")







##### Plot in linear scale (power function)

## instead of linear trends, we should be using power-law relationships for both the PGLS and the weight-support trends.
## In these cases, the relationship is not linear but follows a specific scaling rule :
## f = k*m^(−1/6)
## equivalent to f = k / m^(1/6)


#### find the k value
# Transform the data
x_now <- log(stroke_amp_data$weight)
y_now <- log(stroke_amp_data$stroke_amp)
# Fit pgls model to the transformed data
model <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
# Extract the intercept (log(k)) and slope
intercept <- coef(model)[1]
slope_estimated <- coef(model)[2]
# Calculate k by taking the exponential of the intercept
k_estimated <- exp(intercept)

# Define the scaling function for weight-support (k * m^(-1/6)) (needed in ggplot below)
weight_support_function <- function(m) {
  k <- k_estimated  # Adjust this value to your specific constant 'k'
  return(k / m^(1/6))
}

#### Define the scaling function for PGLS (k * m^(slope)) with the estimated slope
pgls_function <- function(m) {
  k <- k_estimated # Adjust this value to your specific constant 'k'
  slope <- slope_estimated  # Adjust the slope based on your PGLS model
  return(k * m^slope)
}

iso_intercept <- mean(stroke_amp_data$stroke_amp) - 0 * mean(stroke_amp_data$weight)

ggplot(stroke_amp_data, aes(x = weight, y = stroke_amp, label = species)) + 
  geom_point(size=7, color=colors_species_2) + 
  labs(y = 'stroke amplitude [deg]', x = 'body mass [mg]') +
  geom_errorbar(aes(ymax = stroke_amp + stroke_amp_SE, ymin = stroke_amp - stroke_amp_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 1) + 
  ylim(0,150) +
  theme_classic() + 
  theme(legend.position = "none") +
  # Plotting isometry line
  geom_abline(intercept = iso_intercept, slope = (0), color = '#ba4c36', linetype = "dashed", lwd=1.5) + # isometry
  # Plotting the weight-support scaling as a power-law curve
  geom_function(fun = weight_support_function, color = '#A0A0A0', linetype = "dotted", lwd = 1.5) + 
  # Plotting the PGLS scaling as a power-law curve
  geom_function(fun = pgls_function, color = 'black', linetype = "solid", lwd = 1.5)





####### Phylogenetic Reduced Major Axis (phyloRMA):
x_now = log(stroke_amp_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log(stroke_amp_data$stroke_amp); names(y_now) = phy_filmed_only$tip.label
rma_result <- phyl.RMA(x_now, y_now, phy_filmed_only, fixed=T) ; rma_result

### test is it deviate from isometry based on function-embedded test
phyl.RMA(x_now, y_now, phy_filmed_only, fixed=T, h0=0.1)







############ ___ angle-of-attack vs mass ############

AoA_data <- as.data.frame(matrix(NA, nrow(full_dataset), 5))
colnames(AoA_data) = c('species', 'AoA', 'AoA_SE', 'weight', 'weight_SE')
AoA_data$species = full_dataset$genus_species
AoA_data$AoA = full_dataset$AoA_peak_ang_velo_bs
AoA_data$weight = full_dataset$fresh_weight
for (i in 1:length(levels(full_dataset$genus_species))){
  AoA_data$AoA_SE[i] = se(wb_parameters$AoA_peak_ang_velo_bs[which(wb_parameters$genus_species == AoA_data$species[i])])
  AoA_data$weight_SE[i] = se(weight_data$fresh_weight[which(weight_data$genus_species == AoA_data$species[i])])
}

ggplot(AoA_data, aes(x = weight, y = AoA, label = species)) + #geom_text(hjust=0, vjust=0) + 
  geom_point(size=7, color=colors_species) + labs(y = 'angle-of-attack [deg]', x='body mass [mg]') +
  geom_errorbar(aes(ymax = AoA + AoA_SE, ymin = AoA - AoA_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 1) + ylim(0,100) +
  theme_classic() + theme(legend.position = "none")
lmodel2(log10(AoA_data$AoA) ~ log10(AoA_data$weight))


### re-order the rows in 'AoA_data' to match the order of the tips in the phylogeny 'phy_filmed_only'
AoA_data <- AoA_data[match(phy_filmed_only$tip.label, AoA_data$species), ]

## PGLS using nlme on LOG TRANSFORMED DATA
## (following https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(AoA_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log10(AoA_data$AoA); names(y_now) = phy_filmed_only$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)

## Obtain the adjusted R-squared (not given through gls() function so we compute it ourselves)  
residuals <- residuals(pglsModel) # extract residuals
RSS <- sum(residuals^2) # compute RSS (Residual Sum of Squares)
## Compute the Total Sum of Square (TSS)
fitted_values <- pglsModel$fitted # extract fitted values
response_values <- fitted_values + residuals # Reconstruct the response variable
TSS <- sum((response_values - mean(response_values))^2) # Compute Total Sum of Squares (TSS)
# Get number of observations (n) and predictors (k)
n <- length(residuals) # Get number of observations (n)
k <- length(coef(pglsModel)) - 1  # Get number of predictors (k)
## Compute Adjusted R-squared
R2_adj <- 1 - ((RSS / (n - k - 1)) / (TSS / (n - 1))) ; R2_adj
summary(lm(log(AoA_data$AoA) ~ log(AoA_data$weight)))

## Obtain IC 95% of the pgls regression slope 
coef_est <- coef(summary(pglsModel))["x_now", "Value"] # extract coefficient estimates
se_est <- coef(summary(pglsModel))["x_now", "Std.Error"]  # extract standard errors
df <- summary(pglsModel)$dims$N - summary(pglsModel)$dims$p  # get residual degrees of freedom (df) (nb of observation - nb of estimated parameters)
## Compute t-value for 95% CI (alpha = 0.05, two-tailed)
t_value <- qt(0.975, df)  # 97.5th percentile for two-tailed test
## Compute confidence interval
CI_lower <- coef_est - t_value * se_est
CI_upper <- coef_est + t_value * se_est
## Display the Confidence Interval of the slope
c(CI_lower, CI_upper)
## OR SIMPLY USE conint()
confint(pglsModel)


#### plot in log scale
#### including the scaling hypotheses and (non-significant) fitted lines onto the graph

iso_intercept <- mean(log10(AoA_data$AoA)) - 0 * mean(log10(AoA_data$weight))
intercept_W <- mean(log10(AoA_data$AoA)) - (-1/3) * mean(log10(AoA_data$weight))

ggplot(AoA_data, aes(x = log10(weight), y = log10(AoA), label = species)) + 
  geom_point(size=7, color=colors_species_2) + 
  labs(y = 'log10(angle-of-attack)', x='log10(body mass)') +
  geom_errorbar(aes(ymax = log10(AoA + AoA_SE), ymin = log10(AoA - AoA_SE)), width = 0.025) +
  geom_errorbarh(aes(xmax = log10(weight + weight_SE), xmin = log10(weight - weight_SE)), height = 0.005) + 
  theme_classic() + 
  theme(legend.position = "none") +
  ## isometry:
  geom_abline(intercept = iso_intercept, slope = (0), color = '#ba4c36', linetype = "dashed", lwd=1.5) + # isometry
  geom_abline(intercept = 1.7144032, slope = -0.0857708, color = 'black', linetype = "solid", lwd=1.5) + # pgls
  geom_abline(intercept = intercept_W, slope = (-1/3), color = '#A0A0A0', linetype = "dotted", lwd=1.5)# required for weight support if everything else scales isometrically


# Plot in log-log scale while showing real values using Base R
plot(AoA_data$weight, AoA_data$AoA, 
     log='xy', pch=16, cex=3, bty='n', col=colors_species_2, 
     ylab='Angle-of-attack [deg]', xlab='body mass [mg]', axes=F, ylim=c(20, 75))
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

## recompute the PGLS with the log 10 values:
x_log <- log10(AoA_data$weight)
y_log <- log10(AoA_data$AoA)
pglsModel <- gls(y_log ~ x_log, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)
abline(1.7144032, -0.0857708 , col='black', lwd=2) # PGLS
abline(1.7144032,(0), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(1.7144032, (-1/3), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically

## add error bars on the log scale plot showing real values
# Convert standard errors to log space (calculate upper and lower limits)
log_weight_SE_upper <- log10(AoA_data$weight + AoA_data$weight_SE)
log_weight_SE_lower <- log10(AoA_data$weight - AoA_data$weight_SE)
log_AoA_SE_upper <- log10(AoA_data$AoA + AoA_data$AoA_SE)
log_AoA_SE_lower <- log10(AoA_data$AoA - AoA_data$AoA_SE)
# Make sure lower limits are valid (avoid log of zero or negative values)
log_weight_SE_lower[AoA_data$weight - AoA_data$weight_SE <= 0] <- NA
log_stroke_amp_SE_lower[AoA_data$AoA - AoA_data$AoA_SE <= 0] <- NA
# Add horizontal error bars (Body Mass SE)
arrows(x0 = 10^log_weight_SE_lower, x1 = 10^log_weight_SE_upper,
       y0 = AoA_data$AoA, y1 = AoA_data$AoA,
       angle=90, code=3, length=0.05, col="black")

# Add vertical error bars
arrows(x0 = AoA_data$weight, x1 = AoA_data$weight,
       y0 = 10^log_AoA_SE_lower, y1 = 10^log_AoA_SE_upper,
       angle=90, code=3, length=0.05, col="black")






##### Plot in linear scale (power function)

## instead of linear trends, we should be using power-law relationships for both the PGLS and the weight-support trends.
## In these cases, the relationship is not linear but follows a specific scaling rule :
## f = k*m^(−1/3)
## equivalent to f = k / m^(1/3)


#### find the k value
# Transform the data
x_now <- log(AoA_data$weight)
y_now <- log(AoA_data$AoA)
# Fit pgls model to the transformed data
model <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
# Extract the intercept (log(k)) and slope
intercept <- coef(model)[1]
slope_estimated <- coef(model)[2]
# Calculate k by taking the exponential of the intercept
k_estimated <- exp(intercept)

# Define the scaling function for weight-support (k * m^(-1/6)) (needed in ggplot below)
weight_support_function <- function(m) {
  k <- k_estimated  # Adjust this value to your specific constant 'k'
  return(k / m^(1/3))
}

#### Define the scaling function for PGLS (k * m^(slope)) with the estimated slope
pgls_function <- function(m) {
  k <- k_estimated # Adjust this value to your specific constant 'k'
  slope <- slope_estimated  # Adjust the slope based on your PGLS model
  return(k * m^slope)
}

iso_intercept <- mean(AoA_data$AoA) - 0 * mean(AoA_data$weight)

ggplot(AoA_data, aes(x = weight, y = AoA, label = species)) + 
  geom_point(size=7, color=colors_species_2) + 
  labs(y = 'angle-of-attack [deg]', x = 'body mass [mg]') +
  geom_errorbar(aes(ymax = AoA + AoA_SE, ymin = AoA - AoA_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 1) + 
  ylim(0,100) +
  theme_classic() + 
  theme(legend.position = "none") +
  # Plotting isometry line
  geom_abline(intercept = iso_intercept, slope = (0), color = '#ba4c36', linetype = "dashed", lwd=1.5) + # isometry
  # Plotting the weight-support scaling as a power-law curve
  geom_function(fun = weight_support_function, color = '#A0A0A0', linetype = "dotted", lwd = 1.5) + 
  # Plotting the PGLS scaling as a power-law curve
  geom_function(fun = pgls_function, color = 'black', linetype = "solid", lwd = 1.5)





####### Phylogenetic Reduced Major Axis (phyloRMA):
x_now = log(AoA_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log(AoA_data$AoA); names(y_now) = phy_filmed_only$tip.label
rma_result <- phyl.RMA(x_now, y_now, phy_filmed_only, fixed=T) ; rma_result

### test is it deviate from isometry based on function-embedded test
phyl.RMA(x_now, y_now, phy_filmed_only, fixed=T, h0=0.1)









############ ___ deviation amplitude vs mass ############

dev_amp_data <- as.data.frame(matrix(NA, nrow(full_dataset), 5))
colnames(dev_amp_data) = c('species', 'dev_amp', 'dev_amp_SE', 'weight', 'weight_SE')
dev_amp_data$species = full_dataset$genus_species
dev_amp_data$dev_amp = full_dataset$deviation_amplitude
dev_amp_data$weight = full_dataset$fresh_weight
for (i in 1:length(levels(full_dataset$genus_species))){
  dev_amp_data$dev_amp_SE[i] = se(wb_parameters$deviation_amplitude[which(wb_parameters$genus_species == dev_amp_data$species[i])])
  dev_amp_data$weight_SE[i] = se(weight_data$fresh_weight[which(weight_data$genus_species == dev_amp_data$species[i])])
}

ggplot(dev_amp_data, aes(x = weight, y = dev_amp, label = species)) + #geom_text(hjust=0, vjust=0) + 
  geom_point(size=7, color=colors_species) + labs(y = 'deviation amplitude [deg]', x='body mass [mg]') +
  geom_errorbar(aes(ymax = dev_amp + dev_amp_SE, ymin = dev_amp - dev_amp_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 1) + ylim(0,50) +
  theme_classic() + theme(legend.position = "none")
lmodel2(log10(dev_amp_data$dev_amp) ~ log10(dev_amp_data$weight))


### re-order the rows in 'dev_amp_data' to match the order of the tips in the phylogeny 'phy_filmed_only'
dev_amp_data <- dev_amp_data[match(phy_filmed_only$tip.label, dev_amp_data$species), ]

## PGLS using nlme on LOG TRANSFORMED DATA
## (following https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log(dev_amp_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log(dev_amp_data$dev_amp_SE); names(y_now) = phy_filmed_only$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)

## Obtain the adjusted R-squared (not given through gls() function so we compute it ourselves)  
residuals <- residuals(pglsModel) # extract residuals
RSS <- sum(residuals^2) # compute RSS (Residual Sum of Squares)
## Compute the Total Sum of Square (TSS)
fitted_values <- pglsModel$fitted # extract fitted values
response_values <- fitted_values + residuals # Reconstruct the response variable
TSS <- sum((response_values - mean(response_values))^2) # Compute Total Sum of Squares (TSS)
# Get number of observations (n) and predictors (k)
n <- length(residuals) # Get number of observations (n)
k <- length(coef(pglsModel)) - 1  # Get number of predictors (k)
## Compute Adjusted R-squared
R2_adj <- 1 - ((RSS / (n - k - 1)) / (TSS / (n - 1))) ; R2_adj
summary(lm(log(dev_amp_data$dev_amp) ~ log(dev_amp_data$weight)))

## Obtain IC 95% of the pgls regression slope 
coef_est <- coef(summary(pglsModel))["x_now", "Value"] # extract coefficient estimates
se_est <- coef(summary(pglsModel))["x_now", "Std.Error"]  # extract standard errors
df <- summary(pglsModel)$dims$N - summary(pglsModel)$dims$p  # get residual degrees of freedom (df) (nb of observation - nb of estimated parameters)
## Compute t-value for 95% CI (alpha = 0.05, two-tailed)
t_value <- qt(0.975, df)  # 97.5th percentile for two-tailed test
## Compute confidence interval
CI_lower <- coef_est - t_value * se_est
CI_upper <- coef_est + t_value * se_est
## Display the Confidence Interval of the slope
c(CI_lower, CI_upper)
## OR SIMPLY USE conint()
confint(pglsModel)



############ ___ rotation amplitude vs mass ############

rot_amp_data <- as.data.frame(matrix(NA, nrow(full_dataset), 5))
colnames(rot_amp_data) = c('species', 'rot_amp', 'rot_amp_SE', 'weight', 'weight_SE')
rot_amp_data$species = full_dataset$genus_species
rot_amp_data$rot_amp = full_dataset$rotation_amplitude
rot_amp_data$weight = full_dataset$fresh_weight
for (i in 1:length(levels(full_dataset$genus_species))){
  rot_amp_data$rot_amp_SE[i] = se(wb_parameters$rotation_amplitude[which(wb_parameters$genus_species == rot_amp_data$species[i])])
  rot_amp_data$weight_SE[i] = se(weight_data$fresh_weight[which(weight_data$genus_species == rot_amp_data$species[i])])
}

ggplot(rot_amp_data, aes(x = weight, y = rot_amp, label = species)) + #geom_text(hjust=0, vjust=0) + 
  geom_point(size=7, color=colors_species) + labs(y = 'rotation amplitude [deg]', x='body mass [mg]') +
  geom_errorbar(aes(ymax = rot_amp + rot_amp_SE, ymin = rot_amp - rot_amp_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 1) + ylim(0,150) +
  theme_classic() + theme(legend.position = "none")
lmodel2(log10(rot_amp_data$rot_amp) ~ log10(rot_amp_data$weight))


### re-order the rows in 'rot_amp_data' to match the order of the tips in the phylogeny 'phy_filmed_only'
rot_amp_data <- rot_amp_data[match(phy_filmed_only$tip.label, rot_amp_data$species), ]

## PGLS using nlme on LOG TRANSFORMED DATA
## (following https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log(rot_amp_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log(rot_amp_data$rot_amp); names(y_now) = phy_filmed_only$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)

## Obtain the adjusted R-squared (not given through gls() function so we compute it ourselves)  
residuals <- residuals(pglsModel) # extract residuals
RSS <- sum(residuals^2) # compute RSS (Residual Sum of Squares)
## Compute the Total Sum of Square (TSS)
fitted_values <- pglsModel$fitted # extract fitted values
response_values <- fitted_values + residuals # Reconstruct the response variable
TSS <- sum((response_values - mean(response_values))^2) # Compute Total Sum of Squares (TSS)
# Get number of observations (n) and predictors (k)
n <- length(residuals) # Get number of observations (n)
k <- length(coef(pglsModel)) - 1  # Get number of predictors (k)
## Compute Adjusted R-squared
R2_adj <- 1 - ((RSS / (n - k - 1)) / (TSS / (n - 1))) ; R2_adj
summary(lm(log(rot_amp_data$rot_amp) ~ log(rot_amp_data$weight)))

## Obtain IC 95% of the pgls regression slope 
coef_est <- coef(summary(pglsModel))["x_now", "Value"] # extract coefficient estimates
se_est <- coef(summary(pglsModel))["x_now", "Std.Error"]  # extract standard errors
df <- summary(pglsModel)$dims$N - summary(pglsModel)$dims$p  # get residual degrees of freedom (df) (nb of observation - nb of estimated parameters)
## Compute t-value for 95% CI (alpha = 0.05, two-tailed)
t_value <- qt(0.975, df)  # 97.5th percentile for two-tailed test
## Compute confidence interval
CI_lower <- coef_est - t_value * se_est
CI_upper <- coef_est + t_value * se_est
## Display the Confidence Interval of the slope
c(CI_lower, CI_upper)
## OR SIMPLY USE conint()
confint(pglsModel)




############ ___ peak stroke rate vs mass ############

peak_stroke_rate_data <- as.data.frame(matrix(NA, nrow(full_dataset), 5))
colnames(peak_stroke_rate_data) = c('species', 'peak_stroke_rate', 'peak_stroke_rate_SE', 'weight', 'weight_SE')
peak_stroke_rate_data$species = full_dataset$genus_species
peak_stroke_rate_data$peak_stroke_rate = full_dataset$peak_ang_velo_fs
peak_stroke_rate_data$weight = full_dataset$fresh_weight
for (i in 1:length(levels(full_dataset$genus_species))){
  peak_stroke_rate_data$peak_stroke_rate_SE[i] = se(wb_parameters$peak_ang_velo_fs[which(wb_parameters$genus_species == peak_stroke_rate_data$species[i])])
  peak_stroke_rate_data$weight_SE[i] = se(weight_data$fresh_weight[which(weight_data$genus_species == peak_stroke_rate_data$species[i])])
}

ggplot(peak_stroke_rate_data, aes(x = weight, y = peak_stroke_rate, label = species)) + #geom_text(hjust=0, vjust=0) + 
  geom_point(size=7, color=colors_species) + labs(y = 'peak angular velocity [rad/sec]', x='body mass [mg]') +
  geom_errorbar(aes(ymax = peak_stroke_rate + peak_stroke_rate_SE, ymin = peak_stroke_rate - peak_stroke_rate_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 1) + ylim(0,140000) +
  theme_classic() + theme(legend.position = "none")
lmodel2(log10(peak_stroke_rate_data$peak_stroke_rate) ~ log10(peak_stroke_rate_data$weight))



### re-order the rows in 'peak_stroke_rate_data' to match the order of the tips in the phylogeny 'phy_filmed_only'
peak_stroke_rate_data <- peak_stroke_rate_data[match(phy_filmed_only$tip.label, peak_stroke_rate_data$species), ]

## PGLS using nlme on LOG TRANSFORMED DATA
## (following https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log(peak_stroke_rate_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log(peak_stroke_rate_data$peak_stroke_rate); names(y_now) = phy_filmed_only$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)

## Obtain the adjusted R-squared (not given through gls() function so we compute it ourselves)  
residuals <- residuals(pglsModel) # extract residuals
RSS <- sum(residuals^2) # compute RSS (Residual Sum of Squares)
## Compute the Total Sum of Square (TSS)
fitted_values <- pglsModel$fitted # extract fitted values
response_values <- fitted_values + residuals # Reconstruct the response variable
TSS <- sum((response_values - mean(response_values))^2) # Compute Total Sum of Squares (TSS)
# Get number of observations (n) and predictors (k)
n <- length(residuals) # Get number of observations (n)
k <- length(coef(pglsModel)) - 1  # Get number of predictors (k)
## Compute Adjusted R-squared
R2_adj <- 1 - ((RSS / (n - k - 1)) / (TSS / (n - 1))) ; R2_adj
summary(lm(log(peak_stroke_rate_data$peak_stroke_rate) ~ log(peak_stroke_rate_data$weight)))

## Obtain IC 95% of the pgls regression slope 
coef_est <- coef(summary(pglsModel))["x_now", "Value"] # extract coefficient estimates
se_est <- coef(summary(pglsModel))["x_now", "Std.Error"]  # extract standard errors
df <- summary(pglsModel)$dims$N - summary(pglsModel)$dims$p  # get residual degrees of freedom (df) (nb of observation - nb of estimated parameters)
## Compute t-value for 95% CI (alpha = 0.05, two-tailed)
t_value <- qt(0.975, df)  # 97.5th percentile for two-tailed test
## Compute confidence interval
CI_lower <- coef_est - t_value * se_est
CI_upper <- coef_est + t_value * se_est
## Display the Confidence Interval of the slope
c(CI_lower, CI_upper)
## OR SIMPLY USE conint()
confint(pglsModel)




############ ___ peak rotation rate vs mass ############

peak_rotation_rate_data <- as.data.frame(matrix(NA, nrow(full_dataset), 5))
colnames(peak_rotation_rate_data) = c('species', 'peak_rotation_rate', 'peak_rotation_rate_SE', 'weight', 'weight_SE')
peak_rotation_rate_data$species = full_dataset$genus_species
peak_rotation_rate_data$peak_rotation_rate = full_dataset$peak_rot_velo
peak_rotation_rate_data$weight = full_dataset$fresh_weight
for (i in 1:length(levels(full_dataset$genus_species))){
  peak_rotation_rate_data$peak_rotation_rate_SE[i] = se(wb_parameters$peak_rot_velo[which(wb_parameters$genus_species == peak_rotation_rate_data$species[i])])
  peak_rotation_rate_data$weight_SE[i] = se(weight_data$fresh_weight[which(weight_data$genus_species == peak_rotation_rate_data$species[i])])
}

ggplot(peak_rotation_rate_data, aes(x = weight, y = peak_rotation_rate, label = species)) + #geom_text(hjust=0, vjust=0) + 
  geom_point(size=7, color=colors_species) + labs(y = 'peak rotation velocity [rad/sec]', x='body mass [mg]') +
  geom_errorbar(aes(ymax = peak_rotation_rate + peak_rotation_rate_SE, ymin = peak_rotation_rate - peak_rotation_rate_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 1) + ylim(0,140000) +
  theme_classic() + theme(legend.position = "none")
lmodel2(log10(peak_rotation_rate_data$peak_rotation_rate) ~ log10(peak_rotation_rate_data$weight))


### re-order the rows in 'peak_rotation_rate_data' to match the order of the tips in the phylogeny 'phy_filmed_only'
peak_rotation_rate_data <- peak_rotation_rate_data[match(phy_filmed_only$tip.label, peak_rotation_rate_data$species), ]

## PGLS using nlme on LOG TRANSFORMED DATA
## (following https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log(peak_rotation_rate_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log(peak_rotation_rate_data$peak_rotation_rate); names(y_now) = phy_filmed_only$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)

## Obtain the adjusted R-squared (not given through gls() function so we compute it ourselves)  
residuals <- residuals(pglsModel) # extract residuals
RSS <- sum(residuals^2) # compute RSS (Residual Sum of Squares)
## Compute the Total Sum of Square (TSS)
fitted_values <- pglsModel$fitted # extract fitted values
response_values <- fitted_values + residuals # Reconstruct the response variable
TSS <- sum((response_values - mean(response_values))^2) # Compute Total Sum of Squares (TSS)
# Get number of observations (n) and predictors (k)
n <- length(residuals) # Get number of observations (n)
k <- length(coef(pglsModel)) - 1  # Get number of predictors (k)
## Compute Adjusted R-squared
R2_adj <- 1 - ((RSS / (n - k - 1)) / (TSS / (n - 1))) ; R2_adj
summary(lm(log(peak_rotation_rate_data$peak_rotation_rate) ~ log(peak_rotation_rate_data$weight)))

## Obtain IC 95% of the pgls regression slope 
coef_est <- coef(summary(pglsModel))["x_now", "Value"] # extract coefficient estimates
se_est <- coef(summary(pglsModel))["x_now", "Std.Error"]  # extract standard errors
df <- summary(pglsModel)$dims$N - summary(pglsModel)$dims$p  # get residual degrees of freedom (df) (nb of observation - nb of estimated parameters)
## Compute t-value for 95% CI (alpha = 0.05, two-tailed)
t_value <- qt(0.975, df)  # 97.5th percentile for two-tailed test
## Compute confidence interval
CI_lower <- coef_est - t_value * se_est
CI_upper <- coef_est + t_value * se_est
## Display the Confidence Interval of the slope
c(CI_lower, CI_upper)
## OR SIMPLY USE conint()
confint(pglsModel)





############ ___ peak deviation rate vs mass ############

peak_deviation_rate_data <- as.data.frame(matrix(NA, nrow(full_dataset), 5))
colnames(peak_deviation_rate_data) = c('species', 'peak_deviation_rate', 'peak_deviation_rate_SE', 'weight', 'weight_SE')
peak_deviation_rate_data$species = full_dataset$genus_species
peak_deviation_rate_data$peak_deviation_rate = full_dataset$peak_dev_velo
peak_deviation_rate_data$weight = full_dataset$fresh_weight
for (i in 1:length(levels(full_dataset$genus_species))){
  peak_deviation_rate_data$peak_deviation_rate_SE[i] = se(wb_parameters$peak_dev_velo[which(wb_parameters$genus_species == peak_deviation_rate_data$species[i])])
  peak_deviation_rate_data$weight_SE[i] = se(weight_data$fresh_weight[which(weight_data$genus_species == peak_deviation_rate_data$species[i])])
}

ggplot(peak_deviation_rate_data, aes(x = weight, y = peak_deviation_rate, label = species)) + #geom_text(hjust=0, vjust=0) + 
  geom_point(size=7, color=colors_species) + labs(y = 'peak deviation velocity [rad/sec]', x='body mass [mg]') +
  geom_errorbar(aes(ymax = peak_deviation_rate + peak_deviation_rate_SE, ymin = peak_deviation_rate - peak_deviation_rate_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 1) + ylim(0,40000) +
  theme_classic() + theme(legend.position = "none")
lmodel2(log10(peak_deviation_rate_data$peak_deviation_rate) ~ log10(peak_deviation_rate_data$weight))


### re-order the rows in 'peak_deviation_rate_data' to match the order of the tips in the phylogeny 'phy_filmed_only'
peak_deviation_rate_data <- peak_deviation_rate_data[match(phy_filmed_only$tip.label, peak_deviation_rate_data$species), ]

## PGLS using nlme on LOG TRANSFORMED DATA
## (following https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log(peak_deviation_rate_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log(peak_deviation_rate_data$peak_deviation_rate); names(y_now) = phy_filmed_only$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)

## Obtain the adjusted R-squared (not given through gls() function so we compute it ourselves)  
residuals <- residuals(pglsModel) # extract residuals
RSS <- sum(residuals^2) # compute RSS (Residual Sum of Squares)
## Compute the Total Sum of Square (TSS)
fitted_values <- pglsModel$fitted # extract fitted values
response_values <- fitted_values + residuals # Reconstruct the response variable
TSS <- sum((response_values - mean(response_values))^2) # Compute Total Sum of Squares (TSS)
# Get number of observations (n) and predictors (k)
n <- length(residuals) # Get number of observations (n)
k <- length(coef(pglsModel)) - 1  # Get number of predictors (k)
## Compute Adjusted R-squared
R2_adj <- 1 - ((RSS / (n - k - 1)) / (TSS / (n - 1))) ; R2_adj
summary(lm(log(peak_deviation_rate_data$peak_deviation_rate) ~ log(peak_deviation_rate_data$weight)))

## Obtain IC 95% of the pgls regression slope 
coef_est <- coef(summary(pglsModel))["x_now", "Value"] # extract coefficient estimates
se_est <- coef(summary(pglsModel))["x_now", "Std.Error"]  # extract standard errors
df <- summary(pglsModel)$dims$N - summary(pglsModel)$dims$p  # get residual degrees of freedom (df) (nb of observation - nb of estimated parameters)
## Compute t-value for 95% CI (alpha = 0.05, two-tailed)
t_value <- qt(0.975, df)  # 97.5th percentile for two-tailed test
## Compute confidence interval
CI_lower <- coef_est - t_value * se_est
CI_upper <- coef_est + t_value * se_est
## Display the Confidence Interval of the slope
c(CI_lower, CI_upper)
## OR SIMPLY USE conint()
confint(pglsModel)







############ ___ peak angular acceleration vs mass ############

peak_ang_accel_data <- as.data.frame(matrix(NA, nrow(full_dataset), 5))
colnames(peak_ang_accel_data) = c('species', 'peak_ang_accel', 'peak_ang_accel_SE', 'weight', 'weight_SE')
peak_ang_accel_data$species = full_dataset$genus_species
peak_ang_accel_data$peak_ang_accel = full_dataset$peak_ang_accel
peak_ang_accel_data$weight = full_dataset$fresh_weight
for (i in 1:length(levels(full_dataset$genus_species))){
  peak_ang_accel_data$peak_ang_accel_SE[i] = se(wb_parameters$peak_ang_accel[which(wb_parameters$genus_species == peak_ang_accel_data$species[i])])
  peak_ang_accel_data$weight_SE[i] = se(weight_data$fresh_weight[which(weight_data$genus_species == peak_ang_accel_data$species[i])])
}

ggplot(peak_ang_accel_data, aes(x = weight, y = peak_ang_accel, label = species)) + #geom_text(hjust=0, vjust=0) + 
  geom_point(size=7, color=colors_species) + labs(y = 'peak angular acceleration [rad/sec^2]', x='body mass [mg]') +
  geom_errorbar(aes(ymax = peak_ang_accel + peak_ang_accel_SE, ymin = peak_ang_accel - peak_ang_accel_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 1) + ylim(0,2e+08) +
  theme_classic() + theme(legend.position = "none")
lmodel2(log10(peak_ang_accel_data$peak_ang_accel) ~ log10(peak_ang_accel_data$weight))






############ ___ body speed vs mass ############

body_speed_data <- as.data.frame(matrix(NA, nrow(full_dataset), 5))
colnames(body_speed_data) = c('species', 'body_speed', 'body_speed_SE', 'weight', 'weight_SE')
body_speed_data$species = full_dataset$genus_species
body_speed_data$body_speed = full_dataset$mean_body_speed
body_speed_data$weight = full_dataset$fresh_weight
for (i in 1:length(levels(full_dataset$genus_species))){
  body_speed_data$body_speed_SE[i] = se(wb_parameters$mean_body_speed[which(wb_parameters$genus_species == body_speed_data$species[i])])
  body_speed_data$weight_SE[i] = se(weight_data$fresh_weight[which(weight_data$genus_species == body_speed_data$species[i])])
}

ggplot(body_speed_data, aes(x = weight, y = body_speed, label = species)) + #geom_text(hjust=0, vjust=0) + 
  geom_point(size=7, color=colors_species) + labs(y = 'body speed [m/sec]', x='body mass [mg]') +
  geom_errorbar(aes(ymax = body_speed + body_speed_SE, ymin = body_speed - body_speed_SE), width = 2) +
  geom_errorbarh(aes(xmax = weight + weight_SE, xmin = weight - weight_SE), height = 1) + ylim(0,1) +
  theme_classic() + theme(legend.position = "none")
lmodel2(log10(body_speed_data$body_speed) ~ log10(body_speed_data$weight))


lmodel2(log10(full_dataset$abs_wing_speed) ~ log10(full_dataset$fresh_weight) + log10(full_dataset$mean_body_speed))
summary(lm(log10(full_dataset$abs_wing_speed) ~ log10(full_dataset$fresh_weight) + log10(full_dataset$mean_body_speed)))


### re-order the rows in 'body_speed_data' to match the order of the tips in the phylogeny 'phy_filmed_only'
body_speed_data <- body_speed_data[match(phy_filmed_only$tip.label, body_speed_data$species), ]

## PGLS using nlme on LOG TRANSFORMED DATA
## (following https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log(body_speed_data$weight); names(x_now) = phy_filmed_only$tip.label
y_now = log(body_speed_data$body_speed); names(y_now) = phy_filmed_only$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)

summary(lm(log10(body_speed_data$body_speed) ~ log10(body_speed_data$weight)))





############ [3] morphological analysis ############

############ __ [i] scaling between morphology and body mass [individual level] ############

## load the enlarged morphological dataset including the Leiden museum specimens and several measurement per species
morpho_dataset_ind_level <- read.csv(paste0(local_path,'/hoverflies_weight_and_wing_data_ind_level.csv'), sep=',')
morpho_dataset_ind_level$genus_species = as.factor(morpho_dataset_ind_level$genus_species)
morpho_dataset_ind_level$sex = as.factor(morpho_dataset_ind_level$sex)
museum <- which(morpho_dataset_ind_level$origin=='museum')
filmed <- which(morpho_dataset_ind_level$origin=='field')

## check correlation between weight and all wing morphology variables:
chart.Correlation(cbind(morpho_dataset_ind_level$fresh_weight, morpho_dataset_ind_level[,10:ncol(morpho_dataset_ind_level)]), method='pearson', histogram=T, pch=16)

### set species colors
### --> *** there are 20 species from Leiden museum that have not been filmed (morphological analysis only)
###         We have to figure out a way to distinguish them from the other in the colo?r code.

levels(morpho_dataset_ind_level$genus_species)
colors_species = c('grey75', # Ceriana_conopsoides
                   'grey75', # Dorylomorpha_sp
                   '#D1BA15', # Episyrphus_viridaureus
                   'grey75', # Eristalinus_aeneus
                   '#2B3A21', # Eristalis_tenax
                   '#78C155', # Eupeodes_nielseni
                   'grey75', # Helophilus_pendulus
                   'grey75', # Leucozona_lucorum
                   'grey75', # Leucozona_nigripila
                   'grey75', # Melangyna_guttata
                   'grey75', # Melangyna_lasiophthalma
                   'grey75', # Melanostoma_dubium
                   '#A63D40', # Melanostoma_mellinum
                   'grey75', # Microdon_analis
                   'grey75', # Myolepta_dubia
                   'grey75', # Neoascia_annexa
                   'grey75', # Neoascia_geniculata
                   'grey75', # Pipizella_viduata
                   'grey75', # Pipunculus_campestris
                   'grey75', # Platycheirus_albimanus
                   '#7BCEEF', # Platycheirus_clypeatus
                   'grey75', # Psilota_atra
                   'grey75', # Senaspis_elliotti
                   '#674ACE', # Sphaerophoria_scripta
                   '#287C34', # Syrphus_nigrilinearus
                   '#F07E15', # Tropidia_scita
                   'grey75', # Volucella_hyalinipennis
                   'grey75') # Volucella_inanis

colo_species <- colors_species[match(morpho_dataset_ind_level$genus_species, levels(morpho_dataset_ind_level$genus_species))]

# plot legend only:
plot(0, col='white',axes=F, xlab='', ylab='');legend('center', pch=16, pt.cex=1.5, c(levels(morpho_dataset_ind_level$genus_species)), col=colors_species, bty='n', cex=0.6)

### create a color variable where the filmed specimens appear in white just to facilitate there superimposition
colo_species_filmed_invisible = colo_species
colo_species_filmed_invisible[which(colo_species != "grey75")] <- "white"

### create a pch variable for differentiating Male, Female and Unknown
pch_sex <- rep(NA,nrow(morpho_dataset_ind_level))
pch_sex[which(morpho_dataset_ind_level$sex=='M')] <- 15
pch_sex[which(morpho_dataset_ind_level$sex=='F')] <- 16
pch_sex[which(is.na(morpho_dataset_ind_level$sex)==T)] <- 18


### define function to plot detailed log scale while keeping real value:
log10Tck <- function(side, type){
  lim <- switch(side, 
                x = par('usr')[1:2],
                y = par('usr')[3:4],
                stop("side argument must be 'x' or 'y'"))
  at <- floor(lim[1]) : ceil(lim[2])
  return(switch(type, 
                minor = outer(1:9, 10^(min(at):max(at))),
                major = 10^at,
                stop("type argument must be 'major' or 'minor'")
  ))
} ## from: https://stackoverflow.com/questions/47890742/logarithmic-scale-plot-in-r


### !!! --> might be needed to plot as log10(variable), and then add the log axis to ease the slope testing process





############ ___ *** double check S2 and S2* calculation ############ 


## the dimensionlss second-moment-of-area, S2* is define as
## S2' = S2 / ( S * b^2 )
## but in Ty Hedrick MATLAB code, it's computed as
## S2' = sqrt(S2) / (S * b^2)
## which is actually the r2*

## we compute our-self the S2*
morpho_dataset_ind_level$nd_2ndMoment_custom = morpho_dataset_ind_level$secondMoment / (morpho_dataset_ind_level$wing_area_cm2 * (morpho_dataset_ind_level$length_cm)^2)
plot(morpho_dataset_ind_level$nd_2ndMoment_custom , morpho_dataset_ind_level$non_dimensional_2ndMoment , pch=16, cex = 2, col=colo_species, xlab="S2' = S2 / ( S * b^2 )", ylab="WingImageProcessor S2', supposedly equals to sqrt(S2) / (S * b^2)")
## this shows that S2' from Ty, when squared, is equal to S2 / ( S * b^2 )

plot(morpho_dataset_ind_level$nd_2ndMoment_custom , (morpho_dataset_ind_level$non_dimensional_2ndMoment)^2 , pch=16, cex = 2, col=colo_species, xlab="S2' = S2 / ( S * b^2 )", ylab="(WingImageProcessor S2')^2")

## --> WE SHALL USE ( WingImageProcessor S2* values )^2 IN OUR ANALYSIS !!
morpho_dataset_ind_level$non_dimensional_2ndMoment = (morpho_dataset_ind_level$non_dimensional_2ndMoment)^2



############  ________ S2 vs mass: ############ 

plot(log10(morpho_dataset_ind_level$secondMoment)~ log10(morpho_dataset_ind_level$fresh_weight), pch=16, col=colo_species_filmed_invisible, cex=2, bty='n', xlab='log(body mass)',ylab='log(second moment of area)')
points(log10(morpho_dataset_ind_level$secondMoment)[filmed] ~ log10(morpho_dataset_ind_level$fresh_weight)[filmed], pch=21, col='black', bg=colo_species[filmed], cex=2)

## log scale but with actual values:
plot(morpho_dataset_ind_level$secondMoment ~ morpho_dataset_ind_level$fresh_weight, pch=16, col=colo_species_filmed_invisible, cex=2, bty='n', xlab='body mass [mg]',ylab='second moment of area [cm^4]', log='xy', axes=F)
points(morpho_dataset_ind_level$secondMoment[filmed] ~ morpho_dataset_ind_level$fresh_weight[filmed], pch=21, col='black', bg=colo_species[filmed], cex=2)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

lmodel2(log10(morpho_dataset_ind_level$secondMoment) ~ log10(morpho_dataset_ind_level$fresh_weight)) ## check the MA line in the output

## because S2 ~ body mass ^(4/3)", the slope should be 4/3 if S2 scales isometrically with mass
abline(-3.093464 , 1.134511  , col='black', lwd=2)          # observed OLS
iso_intercept <- mean(log10(morpho_dataset_ind_level$secondMoment)) - 4/3 * mean(log10(morpho_dataset_ind_level$fresh_weight))
abline(iso_intercept ,(4/3), col='grey35', lwd=2, lty='dashed') # isometric
legend('topleft', c('observed scaling (1.24)', 'isometric scaling (1.33)'), col=c('black', 'grey35'), bty='n',lty=c('solid','dashed'),cex=0.7,pt.cex=1.2)
## because 1.24 < 1.333 (i.e. 4/3), we have a negative allometry.
## -> this means that larger hoverflies get relatively smaller S2 for their size,
## test whether the two slopes differ significantly
sma(log10(morpho_dataset_ind_level$secondMoment) ~ log10(morpho_dataset_ind_level$fresh_weight), slope.test = 4/3)
## SIGNIFICANT: NEGATIVE ALLOMETRY

## plot for Table summing up scaling expectation:
plot(log10(morpho_dataset_ind_level$secondMoment) ~ log10(morpho_dataset_ind_level$fresh_weight), pch=16, col='grey80', cex=1.5, bty='L', xlab='log(body mass)',ylab='log(second moment of area)', asp=F)
abline(-3.093464 , 1.134511 , col='black', lwd=4.5)  # observed OLS
abline(a = iso_intercept, b = 4/3, col = '#ba4c36', lwd=4.5, lty='dashed') # isometric
intercept_W <- mean(log10(morpho_dataset_ind_level$secondMoment)) - 1 * mean(log10(morpho_dataset_ind_level$fresh_weight))
abline(a = intercept_W, b = 1, col = '#0e9594', lwd=4.5, lty='solid') # weight support



############  ________ wing length vs mass: ############ 

## body mass scales with (wing length)^3 <=> wing length ~ body mass ^1/3 (because volume is length ^3)
## for wing length to scale isometrically with size: wing length ~ body mass ^(1/3)

## first plot the museum specimen only
plot(log10(morpho_dataset_ind_level$length_cm) ~ log10(morpho_dataset_ind_level$fresh_weight), pch=16, col=colo_species_filmed_invisible, cex=2, bty='n', xlab='log(body mass)',ylab='log(wing length)')
## superimpose the filmed specimens:
points(log10(morpho_dataset_ind_level$length_cm)[filmed] ~ log10(morpho_dataset_ind_level$fresh_weight)[filmed], pch=21, col='black', bg=colo_species[filmed], cex=2)

## log scale but with actual values:
plot(morpho_dataset_ind_level$length_cm ~ morpho_dataset_ind_level$fresh_weight, pch=16, col=colo_species_filmed_invisible, cex=2, bty='n', xlab='body mass [mg]',ylab='wing length [cm]', log='xy', axes=F)
points(morpho_dataset_ind_level$length_cm[filmed] ~ morpho_dataset_ind_level$fresh_weight[filmed], pch=21, col='black', bg=colo_species[filmed], cex=2)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)


lmodel2(log10(morpho_dataset_ind_level$length_cm) ~ log10(morpho_dataset_ind_level$fresh_weight))
## check both the OLS and the MA line (are they consistent?)
## (see "Warton 2007. Bivariate line fitting methods for allometry" about RMA)

## because "wing length ~ body mass ^(1/3)", the slope should be 1/3 if wingspan scales isometrically with mass
abline(-0.5091325, 0.2872195  , col='black', lwd=2)         # observed (MA row of the lmodel2() output)
## calculate the isometric intercept to make the line pass by the mean value of the y variable
## this can be done as intercept = mean_y - slope * mean_x
iso_intercept <- mean(log10(morpho_dataset_ind_level$length_cm)) - 1/3 * mean(log10(morpho_dataset_ind_level$fresh_weight))
abline(iso_intercept,(1/3), col='grey35', lwd=2, lty='dashed') # isometric
legend('topleft', c('observed scaling (0.28)', 'isometric scaling (0.33)'), col=c('black', 'grey35'), bty='n',lty=c('solid','dashed'),cex=0.7,pt.cex=1.2)
## because 0.28 < 0.33 (i.e. 1/3), we have a negative allometry.
## -> this means that larger hoverflies get relatively smaller wing length for their size
## test whether the two slopes differ significantly:
sma(log10(morpho_dataset_ind_level$length_cm) ~ log10(morpho_dataset_ind_level$fresh_weight), slope.test = 1/3)
### SIGNIFICANT: NEGATIVE ALLOMETRY
## moreover, the confidence interval given by 'lmodel2' exclude 1/3


## plot for Table summing up scaling expectation:
plot(log10(morpho_dataset_ind_level$length_cm) ~ log10(morpho_dataset_ind_level$fresh_weight), pch=16, col='grey80', cex=1.5, bty='L', xlab='log(body mass)',ylab='log(wing length)', asp=F)
abline(-0.5091325, 0.2872195  , col='black', lwd=4.5)  # observed
abline(a = iso_intercept, b = 1/3, col = '#ba4c36', lwd=4.5, lty='dashed') # isometric
intercept_W <- mean(log10(morpho_dataset_ind_level$length_cm)) - 2/9 * mean(log10(morpho_dataset_ind_level$fresh_weight))
abline(a = intercept_W, b = 2/9, col = '#0e9594', lwd=4.5, lty='solid') # weight support




############  ________ wing area vs mass: ############ 

plot(log10(morpho_dataset_ind_level$wing_area_cm2) ~ log10(morpho_dataset_ind_level$fresh_weight), pch=16, col=colo_species_filmed_invisible, cex=2, bty='n', xlab='log(body mass)',ylab='log(wing area)')
points(log10(morpho_dataset_ind_level$wing_area_cm2)[filmed] ~ log10(morpho_dataset_ind_level$fresh_weight)[filmed], pch=21, col='black', bg=colo_species[filmed], cex=2)

## log scale but with actual values:
plot(morpho_dataset_ind_level$wing_area_cm2 ~ morpho_dataset_ind_level$fresh_weight, pch=16, col=colo_species_filmed_invisible, cex=2, bty='n', xlab='body mass [mg]',ylab='wing area [cm^2]', log='xy', axes=F)
points(morpho_dataset_ind_level$wing_area_cm2[filmed] ~ morpho_dataset_ind_level$fresh_weight[filmed], pch=21, col='black', bg=colo_species[filmed], cex=2)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

lmodel2(log10(morpho_dataset_ind_level$wing_area_cm2) ~ log10(morpho_dataset_ind_level$fresh_weight))
## check both the OLS and the MA line (are they consistent?)
## (see "Warton 2007. Bivariate line fitting methods for allometry" about RMA)

## because wing area ~ body mass ^(2/3)", the slope should be 2/3 if wing area scales isometrically with mass
abline(-1.701080 ,0.6344530 , col='black', lwd=2)          # observed
iso_intercept <- mean(log10(morpho_dataset_ind_level$wing_area_cm2)) - 2/3 * mean(log10(morpho_dataset_ind_level$fresh_weight))
abline(iso_intercept,(2/3), col='grey35', lwd=2, lty='dashed') # isometric 
legend('topleft', c('observed scaling (0.63)', 'isometric scaling (0.66)'), col=c('black', 'grey35'), bty='n',lty=c('solid','dashed'),cex=0.7,pt.cex=1.2)
## because 0.63 < 0.6666 (i.e. 2/3), we have a negative allometry.
## test whether the two slopes differ significantly
sma(log10(morpho_dataset_ind_level$wing_area_cm2) ~ log10(morpho_dataset_ind_level$fresh_weight), slope.test = 2/3)
### NON SIGNIFICANT --> area scale isometrically


## plot for Table summing up scaling expectation:
plot(log10(morpho_dataset_ind_level$wing_area_cm2) ~ log10(morpho_dataset_ind_level$fresh_weight), pch=16, col='grey80', cex=1.5, bty='L', xlab='log(body mass)',ylab='log(wing area)', asp=F)
abline(-1.701080 ,0.6344530  , col='black', lwd=4.5)  # observed
abline(a = iso_intercept, b = 2/3, col = '#ba4c36', lwd=4.5, lty='dashed') # isometric
intercept_W <- mean(log10(morpho_dataset_ind_level$wing_area_cm2)) - 1/3 * mean(log10(morpho_dataset_ind_level$fresh_weight))
abline(a = intercept_W, b = 1/3, col = '#0e9594', lwd=4.5, lty='solid') # weight support






############  ________ aspect ratio vs mass: ############ 

plot(log10(morpho_dataset_ind_level$aspect_ratio) ~ log10(morpho_dataset_ind_level$fresh_weight), pch=16, col=colo_species_filmed_invisible, cex=2, bty='n', xlab='log(body mass)',ylab='log(aspect ratio)')
points(log10(morpho_dataset_ind_level$aspect_ratio)[filmed] ~ log10(morpho_dataset_ind_level$fresh_weight)[filmed], pch=21, col='black', bg=colo_species[filmed], cex=2)

## log scale but with actual values:
plot(morpho_dataset_ind_level$aspect_ratio ~ morpho_dataset_ind_level$fresh_weight, pch=16, col=colo_species_filmed_invisible, cex=2, bty='n', xlab='body mass [mg]',ylab='aspect ratio', log='xy', axes=F)
points(morpho_dataset_ind_level$aspect_ratio[filmed] ~ morpho_dataset_ind_level$fresh_weight[filmed], pch=21, col='black', bg=colo_species[filmed], cex=2)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

lmodel2(log10(morpho_dataset_ind_level$aspect_ratio) ~ log10(morpho_dataset_ind_level$fresh_weight))
## check both the OLS and the MA line (are they consistent?)
## (see "Warton 2007. Bivariate line fitting methods for allometry" about RMA)

## because aspect ratio ~ body mass ^(0)", the slope should be 0 if wing aspect ratio scales isometrically with mass
abline(0.9620381 , -0.04360696 , col='black', lwd=2)  # observed
iso_intercept <- mean(log10(morpho_dataset_ind_level$aspect_ratio)) - 0 * mean(log10(morpho_dataset_ind_level$fresh_weight))
abline(iso_intercept, 0, col='grey35', lwd=2, lty='dashed')   # isometric
legend('topright', c('observed scaling (-0.04)', 'isometric scaling (0)'), col=c('black', 'grey35'), bty='n',lty=c('solid','dashed'),cex=0.7,pt.cex=1.2)
## because -0.04360696 < 0, we have a negative allometry.
## test whether the two slopes differ significantly
sma(log10(morpho_dataset_ind_level$aspect_ratio) ~ log10(morpho_dataset_ind_level$fresh_weight), slope.test = 0)
### SIGNIFICANT: NEGATIVE ALLOMETRY
## moreover, the confidence interval given by 'lmodel2' exclude 0





###########  ________ wing chord vs mass: ############ 

## body mass scales with (wing chord)^3 <=> wing chord ~ body mass ^1/3 (because volume is wing chord ^3)
## for wing chord to scale isometrically with size: wing chord ~ body mass ^(1/3)

plot(log10(morpho_dataset_ind_level$average_chord)~ log10(morpho_dataset_ind_level$fresh_weight), pch=16, col=colo_species_filmed_invisible, cex=2, bty='n', xlab='log(body mass)',ylab='log(wing chord)')
points(log10(morpho_dataset_ind_level$average_chord)[filmed] ~ log10(morpho_dataset_ind_level$fresh_weight)[filmed], pch=21, col='black', bg=colo_species[filmed], cex=2)

## log scale but with actual values:
plot(morpho_dataset_ind_level$average_chord ~ morpho_dataset_ind_level$fresh_weight, pch=16, col=colo_species_filmed_invisible, cex=2, bty='n', xlab='body mass [mg]',ylab='wing chord [cm]', log='xy', axes=F)
points(morpho_dataset_ind_level$average_chord[filmed] ~ morpho_dataset_ind_level$fresh_weight[filmed], pch=21, col='black', bg=colo_species[filmed], cex=2)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

lmodel2(log10(morpho_dataset_ind_level$average_chord) ~ log10(morpho_dataset_ind_level$fresh_weight)) ## check the MA line in the output

## because "wing chord ~ body mass ^(1/3)", the slope should be 1/3 if chord scales isometrically with mass
abline(-1.170615 , 0.3315192  , col='black', lwd=2)           # observed
iso_intercept <- mean(log10(morpho_dataset_ind_level$average_chord)) - 1/3 * mean(log10(morpho_dataset_ind_level$fresh_weight))
abline(iso_intercept ,(1/3), col='grey35', lwd=2, lty='dashed')  # isometric
legend('topleft', c('observed scaling (0.33)', 'isometric scaling (0.33)'), col=c('black', 'grey35'), bty='n',lty=c('solid','dashed'),cex=0.7,pt.cex=1.2)
sma(log10(morpho_dataset_ind_level$average_chord) ~ log10(morpho_dataset_ind_level$fresh_weight), slope.test = 1/3)

## NON SIGNIFICANT: isometric scaling


## plot for Table summing up scaling expectation:
plot(log10(morpho_dataset_ind_level$average_chord) ~ log10(morpho_dataset_ind_level$fresh_weight), pch=16, col='grey80', cex=1.5, bty='L', xlab='log(body mass)',ylab='log(wing chord)', asp=F)
abline(-1.170615 , 0.3315192  , col='black', lwd=4.5)  # observed
abline(a = iso_intercept, b = 1/3, col = '#ba4c36', lwd=4.5, lty='dashed') # isometric
intercept_W <- mean(log10(morpho_dataset_ind_level$average_chord)) - 2/3 * mean(log10(morpho_dataset_ind_level$fresh_weight))
abline(a = intercept_W, b = 2/3, col = '#0e9594', lwd=4.5, lty='solid') # weight support





############  ________ non-dim S2 vs mass: ############ 

plot(log10(morpho_dataset_ind_level$non_dimensional_2ndMoment)~ log10(morpho_dataset_ind_level$fresh_weight), pch=16, col=colo_species_filmed_invisible, cex=2, bty='n', xlab='log(body mass)',ylab='log(second moment of area)')
points(log10(morpho_dataset_ind_level$non_dimensional_2ndMoment)[filmed] ~ log10(morpho_dataset_ind_level$fresh_weight)[filmed], pch=21, col='black', bg=colo_species[filmed], cex=2)

## log scale but with actual values:
plot(morpho_dataset_ind_level$non_dimensional_2ndMoment ~ morpho_dataset_ind_level$fresh_weight, pch=16, col=colo_species_filmed_invisible, cex=2, bty='n', xlab='body mass [mg]',ylab='dimensionless S2 [-]', log='xy', axes=F)
points(morpho_dataset_ind_level$non_dimensional_2ndMoment[filmed] ~ morpho_dataset_ind_level$fresh_weight[filmed], pch=21, col='black', bg=colo_species[filmed], cex=2)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

lmodel2(log10(morpho_dataset_ind_level$non_dimensional_2ndMoment) ~ log10(morpho_dataset_ind_level$fresh_weight)) ## check the MA line in the output

## because nd_S2 ~ body mass ^(0)", the slope should be 0 if S2 scales isometrically with mass
abline(-0.2089225 , -0.02107210 , col='black', lwd=2)    # observed
iso_intercept <- mean(log10(morpho_dataset_ind_level$non_dimensional_2ndMoment)) - 0 * mean(log10(morpho_dataset_ind_level$fresh_weight))
abline(iso_intercept,(0), col='grey35', lwd=2, lty='dashed') # isometric
legend('top', c('observed scaling (-0.021)', 'isometric scaling (0)'), col=c('black', 'grey35'), bty='n',lty=c('solid','dashed'),cex=0.7,pt.cex=1.2)
## because -0.021 < 0, we have a negative allometry.
## -> this means that larger hoverflies get relatively smaller non dimensional S2 for their size
## test whether the two slopes differ significantly 
sma(log10(morpho_dataset_ind_level$non_dimensional_2ndMoment) ~ log10(morpho_dataset_ind_level$fresh_weight), slope.test = 0)
### SIGNIFICANT
## moreover, the confidence interval given by 'lmodel2' exclude 0



## plot for Table summing up scaling expectation:
plot(log10(morpho_dataset_ind_level$non_dimensional_2ndMoment) ~ log10(morpho_dataset_ind_level$fresh_weight), pch=16, col='grey80', cex=1.5, bty='L', xlab='log(body mass)',ylab='log(dimensionless S2)', asp=F)
abline(-0.2089225 , -0.02107210  , col='black', lwd=4.5)  # observed
abline(a = iso_intercept, b = 0, col = '#ba4c36', lwd=4.5, lty='dashed') # isometric
intercept_W <- mean(log10(morpho_dataset_ind_level$non_dimensional_2ndMoment)) - (-1/3) * mean(log10(morpho_dataset_ind_level$fresh_weight))
abline(a = intercept_W, b = -1/3, col = '#0e9594', lwd=4.5, lty='solid') # weight support





############ [4] effect of sex on morphology ############

## check the number of M and F in the different species
as.data.frame.matrix(table(morpho_dataset_ind_level$genus_species, morpho_dataset_ind_level$sex))

### look for difference between sex within species

############  ________ body mass ############ 
anova(lm(morpho_dataset_ind_level$fresh_weight ~ morpho_dataset_ind_level$sex + morpho_dataset_ind_level$genus_species))
## Mixed Model with species as a random effect is the preferred way to control for the species effect:
model <- lmer(fresh_weight ~ sex + (1 | genus_species), data = morpho_dataset_ind_level) ; summary(model)

### significantly different between sex
### investigate the direction of that difference:
# Compute mean fresh weight per species and sex
mean_weights <- aggregate(fresh_weight ~ genus_species + sex, data = morpho_dataset_ind_level, FUN = mean)
# Reshape data: Spread into separate columns for Male and Female weights
library(tidyr)
mean_weights_wide <- spread(mean_weights, key = sex, value = fresh_weight)
# Compute the difference (Female - Male)
mean_weights_wide$diff_female_male <- mean_weights_wide$F - mean_weights_wide$M
# Add a column indicating which sex has the higher mean weight
mean_weights_wide$largest <- ifelse(mean_weights_wide$diff_female_male > 0, "female", "male")
# remove NAs
mean_weights_wide <- na.omit(mean_weights_wide)

# Count the number of species where females show larger values (X) and male show larger values (Y)
X <- sum(mean_weights_wide$largest == "female", na.rm = TRUE)
Y <- sum(mean_weights_wide$largest == "male", na.rm = TRUE)
cat("Body mass was higher in females for", X, "species and higher in males for", Y, "species.\n")



############  ________ wingspan ############ 
anova(lm(morpho_dataset_ind_level$length_cm ~ morpho_dataset_ind_level$sex + morpho_dataset_ind_level$genus_species))
## non significant
## Mixed Model:
model <- lmer(length_cm ~ sex + (1 | genus_species), data = morpho_dataset_ind_level) ; summary(model)


############  ________ wing chord ############ 
anova(lm(morpho_dataset_ind_level$average_chord ~ morpho_dataset_ind_level$sex + morpho_dataset_ind_level$genus_species))
## Mixed Model:
model <- lmer(average_chord ~ sex + (1 | genus_species), data = morpho_dataset_ind_level) ; summary(model)


### significantly different between sex
### investigate the direction of that difference:
mean_chord <- aggregate(average_chord ~ genus_species + sex, data = morpho_dataset_ind_level, FUN = mean)
# Reshape data: Spread into separate columns for Male and Female weights
mean_chord_wide <- spread(mean_chord, key = sex, value = average_chord)
# Compute the difference (Female - Male)
mean_chord_wide$diff_female_male <- mean_chord_wide$F - mean_chord_wide$M
# Add a column indicating which sex has the higher mean weight
mean_chord_wide$largest <- ifelse(mean_chord_wide$diff_female_male > 0, "female", "male")
# remove NAs
mean_chord_wide <- na.omit(mean_chord_wide)

# Count the number of species where females show larger values (X) and male show larger values (Y)
X <- sum(mean_chord_wide$largest == "female", na.rm = TRUE)
Y <- sum(mean_chord_wide$largest == "male", na.rm = TRUE)
cat("Mean chord was higher in females for", X, "species and higher in males for", Y, "species.\n")




############  ________ wing area ############ 
anova(lm(morpho_dataset_ind_level$wing_area_cm2 ~ morpho_dataset_ind_level$sex + morpho_dataset_ind_level$genus_species))
## Mixed Model:
model <- lmer(wing_area_cm2 ~ sex + (1 | genus_species), data = morpho_dataset_ind_level) ; summary(model)

### significantly different between sex
### investigate the direction of that difference:
mean_area <- aggregate(wing_area_cm2 ~ genus_species + sex, data = morpho_dataset_ind_level, FUN = mean)
# Reshape data: Spread into separate columns for Male and Female weights
mean_area_wide <- spread(mean_area, key = sex, value = wing_area_cm2)
# Compute the difference (Female - Male)
mean_area_wide$diff_female_male <- mean_area_wide$F - mean_area_wide$M
# Add a column indicating which sex has the higher mean weight
mean_area_wide$largest <- ifelse(mean_area_wide$diff_female_male > 0, "female", "male")
# remove NAs
mean_area_wide <- na.omit(mean_area_wide)

# Count the number of species where females show larger values (X) and male show larger values (Y)
X <- sum(mean_area_wide$largest == "female", na.rm = TRUE)
Y <- sum(mean_area_wide$largest == "male", na.rm = TRUE)
cat("Wing area was larger in females for", X, "species and higher in males for", Y, "species.\n")



############  ________ S2 ############ 
anova(lm(morpho_dataset_ind_level$secondMoment ~ morpho_dataset_ind_level$sex + morpho_dataset_ind_level$genus_species))
## Mixed Model:
model <- lmer(secondMoment ~ sex + (1 | genus_species), data = morpho_dataset_ind_level) ; summary(model)


### significantly different between sex
### investigate the direction of that difference:
mean_S2 <- aggregate(secondMoment ~ genus_species + sex, data = morpho_dataset_ind_level, FUN = mean)
# Reshape data: Spread into separate columns for Male and Female weights
mean_S2_wide <- spread(mean_S2, key = sex, value = secondMoment)
# Compute the difference (Female - Male)
mean_S2_wide$diff_female_male <- mean_S2_wide$F - mean_S2_wide$M
# Add a column indicating which sex has the higher mean weight
mean_S2_wide$largest <- ifelse(mean_S2_wide$diff_female_male > 0, "female", "male")
# remove NAs
mean_S2_wide <- na.omit(mean_S2_wide)

# Count the number of species where females show larger values (X) and male show larger values (Y)
X <- sum(mean_S2_wide$largest == "female", na.rm = TRUE)
Y <- sum(mean_S2_wide$largest == "male", na.rm = TRUE)
cat("S2 was larger in females for", X, "species and higher in males for", Y, "species.\n")


############  ________ non-dim S2 ############ 
anova(lm(morpho_dataset_ind_level$non_dimensional_2ndMoment ~ morpho_dataset_ind_level$sex + morpho_dataset_ind_level$genus_species))
## non significant
## Mixed Model:
model <- lmer(non_dimensional_2ndMoment ~ sex + (1 | genus_species), data = morpho_dataset_ind_level) ; summary(model)




############ [5] phylogenetic analyses ############


########## check number of measurements per species

## nb of species
length(levels(morpho_dataset_ind_level$genus_species))
## number of individuals:
nrow(morpho_dataset_ind_level)
## number of individuals per species
compile_inds_per_species = as.data.frame(matrix(NA, length(levels(morpho_dataset_ind_level$genus_species)), 2))
colnames(compile_inds_per_species) = c('species', 'nb_of_inds')
compile_inds_per_species$species = as.factor(levels(morpho_dataset_ind_level$genus_species))
for (i in 1:nrow(compile_inds_per_species)){
  species_now = compile_inds_per_species$species[i]
  subset_now = morpho_dataset_ind_level[which(morpho_dataset_ind_level$genus_species==species_now),]
  subset_now[] <- lapply(subset_now, function(x) if(is.factor(x)) factor(x) else x) # drop levels
  compile_inds_per_species$nb_of_inds[i] = length(subset_now$individual_ID)
}
## mean nb of individuals per species = 4.21 +/- 1.75
mean(compile_inds_per_species$nb_of_inds)
sd(compile_inds_per_species$nb_of_inds)



################ ____ create species-level dataset ################

morpho_dataset_species_level <- as.data.frame(matrix(NA, length(levels(morpho_dataset_ind_level$genus_species)), length(5:ncol(morpho_dataset_ind_level))+1))
colnames(morpho_dataset_species_level) = c('genus_species', colnames(morpho_dataset_ind_level)[5:ncol(morpho_dataset_ind_level)])
morpho_dataset_species_level$genus_species = as.factor(levels(morpho_dataset_ind_level$genus_species))

## fill the categorical variables:
for (i in 1:nrow(morpho_dataset_species_level)){
  morpho_dataset_species_level$origin[i] = as.character(morpho_dataset_ind_level$origin[which(morpho_dataset_ind_level$genus_species==morpho_dataset_species_level$genus_species[i])][1])
  morpho_dataset_species_level$filmed[i] = as.character(morpho_dataset_ind_level$filmed[which(morpho_dataset_ind_level$genus_species==morpho_dataset_species_level$genus_species[i])][1])
}
## fill the morphological parameters:
for (i in 4:ncol(morpho_dataset_species_level)){
  morpho_dataset_species_level[,i] = tapply(morpho_dataset_ind_level[,which(colnames(morpho_dataset_ind_level)==colnames(morpho_dataset_species_level)[i])], morpho_dataset_ind_level$genus_species, mean)
}

## add normalised wingspan: R* = R / m^0.33
morpho_dataset_species_level$normalised_R <- rep(NA,nrow(morpho_dataset_species_level))
for (i in 1:nrow(morpho_dataset_species_level)) {
  morpho_dataset_species_level$normalised_R[i] = (morpho_dataset_species_level$length_cm[i]/100) / ( (morpho_dataset_species_level$fresh_weight[i]/1e+06)^0.33 )
}
# View(morpho_dataset_species_level[,c(1,15, 6)])

### save "morpho_dataset_species_level"
# write.csv(morpho_dataset_species_level, paste0(local_path, '/hoverflies_weight_and_wing_data_species_level.csv'), row.names = F)




### re-order the rows in "morpho_dataset_species_level" such as matching the order of the tips in the phylogeny
morpho_dataset_species_level_tip_order = morpho_dataset_species_level
rownames(morpho_dataset_species_level_tip_order) = phy$tip.label
for (i in 1:nrow(morpho_dataset_species_level_tip_order)){
  morpho_dataset_species_level_tip_order[i,] = morpho_dataset_species_level[which(morpho_dataset_species_level$genus_species==phy$tip.label[i]),]
}

### re-order the rows in "full_dataset" such as matching the order of the tips in the phylogeny
full_dataset_tip_order = full_dataset
rownames(full_dataset_tip_order) = phy_filmed_only$tip.label
for (i in 1:nrow(full_dataset_tip_order)){
  full_dataset_tip_order[i,] = full_dataset[which(full_dataset$genus_species==phy_filmed_only$tip.label[i]),]
}

museum_tip_order <- which(morpho_dataset_species_level_tip_order$origin=='museum')
filmed_tip_order <- which(morpho_dataset_species_level_tip_order$origin=='field')




############ __ [i] phylogenetic signal ############


######### ____ [i.a] phylogenetic signal on morphological trait ######### 

data_phylosig_test <- morpho_dataset_species_level_tip_order[,c(6:ncol(morpho_dataset_species_level_tip_order))]

### we use the function 'phylosig' from Phytools to compute phylogenetic signal of continuous traits.
### and we store the K of Blomberg and the p-value in a table
phylosig_result = as.data.frame(matrix(NA, length(colnames(data_phylosig_test)), 2))
colnames(phylosig_result) = c('Blomberg K', 'p-value')
rownames(phylosig_result) = colnames(data_phylosig_test)

### filling results in "phylosig_result"
## ----> method Blomberg's K
for (i in 1:nrow(phylosig_result)){
  p = phylosig(phy, data_phylosig_test[,which(colnames(data_phylosig_test)==rownames(phylosig_result)[i])], method="K",test=T)
  phylosig_result[i,1] = p$K
  phylosig_result[i,2] = p$P
}  ; colnames(phylosig_result)[1] = c('Blomberg K')

## ----> method Pagel's lambda (might be better for phylo signal below 1 according to Revell and Harmon 2022. Phylogenetic Comparative Methods in R; Princeton University)
# for (i in 1:nrow(phylosig_result)){
#   p = phylosig(phy, data_phylosig_test[,which(colnames(data_phylosig_test)==rownames(phylosig_result)[i])], method="lambda",test=T)
#   phylosig_result[i,1] = p$lambda
#   phylosig_result[i,2] = p$P
# } ; colnames(phylosig_result)[1] = c('Pagel lambda')


## order by decreasing values of phylo signal
phylosig_result_ordered = phylosig_result[order(phylosig_result[,1],decreasing = T),]

## --> All morphological parameters show significant phylogenetic signal,
##     excepted the 2nd and 3rd moment of wing area.
## The signal is  strong (lambda almost 1), trait variation perfectly matches the phylogenetic structure.



######### ____ [i.b] phylogenetic signal on flight trait ######### 

data_phylosig_test <- full_dataset_tip_order[,c(2:20)]

### we use the function 'phylosig' from Phytools to compute phylogenetic signal of continuous traits.
### and we store the K of Blomberg and the p-value in a table
phylosig_result = as.data.frame(matrix(NA, length(colnames(data_phylosig_test)), 2))
colnames(phylosig_result) = c('Blomberg K', 'p-value')
rownames(phylosig_result) = colnames(data_phylosig_test)

### filling results in "phylosig_result"
## ----> method Blomberg's K
for (i in 1:nrow(phylosig_result)){
  p = phylosig(phy_filmed_only, data_phylosig_test[,which(colnames(data_phylosig_test)==rownames(phylosig_result)[i])], method="K",test=T)
  phylosig_result[i,1] = p$K
  phylosig_result[i,2] = p$P
}  ; colnames(phylosig_result)[1] = c('Blomberg K')

## ----> method Pagel's lambda (might be better for phylo signal below 1 according to Revell and Harmon 2022. Phylogenetic Comparative Methods in R; Princeton University)
# for (i in 1:nrow(phylosig_result)){
#   p = phylosig(phy_filmed_only, data_phylosig_test[,which(colnames(data_phylosig_test)==rownames(phylosig_result)[i])], method="lambda",test=T)
#   phylosig_result[i,1] = p$lambda
#   phylosig_result[i,2] = p$P
# } ; colnames(phylosig_result)[1] = c('Pagel lambda')


## order by decreasing values of phylo signal
phylosig_result_ordered = phylosig_result[order(phylosig_result[,1],decreasing = T),]

## --> No phylogenetic  signal detected on the flight parameter (could be due to small sample size...)




############ __ [ii] scaling between morphology and body mass [phylogenetic level] ############


##### set colors for species
tip_cols <- colors_species[match(morpho_dataset_species_level_tip_order$genus_species, levels(morpho_dataset_species_level_tip_order$genus_species))]

# plot legend only:
plot(0, col='white',axes=F, xlab='', ylab='');legend('center', pch=16, pt.cex=1.5, c(levels(morpho_dataset_species_level_tip_order$genus_species)), col=colors_species, bty='n', cex=0.6)

### create a color variable where the filmed specimens appear in white just to facilitate there superimposition
colo_species_filmed_invisible_tip_order = tip_cols
colo_species_filmed_invisible_tip_order[which(tip_cols != "grey75")] <- "white"



############  ________ *** plot body mass at the tip of the phylogeny ############ 

plot(phy, cex=0.6, edge.width = 0.5,  use.edge.length = F)
phy_2 <- phy
phy_2$edge.length <- NULL # this is done to ignore edge length on the figure

plotTree.barplot(phy_2, setNames(morpho_dataset_species_level_tip_order$fresh_weight, rownames(morpho_dataset_species_level_tip_order)),
                 args.plotTree=list(ftype="off", lwd=0.9),
                 args.barplot=list(xlab="body mass [mg]", space=1, border='grey30'))



##### compute standard error for body mass measurements
dataset_ind_lvl <- read.csv(paste0(local_path,'/hoverflies_weight_and_wing_data_ind_level.csv'), sep=',')

### compute SD using the tapply function:
species_sd <- tapply(dataset_ind_lvl$fresh_weight, dataset_ind_lvl$genus_species, sd, na.rm = T)
## re-order to match "morpho_dataset_species_level_tip_order"
species_sd_tip_order <- species_sd[morpho_dataset_species_level_tip_order$genus_species]
morpho_dataset_species_level_tip_order$fresh_weight_SD <- species_sd_tip_order # add it to the main dataset


#### plot including the SD bars
fresh_weight <- setNames(morpho_dataset_species_level_tip_order$fresh_weight, 
                         rownames(morpho_dataset_species_level_tip_order))
fresh_weight_sd <- setNames(morpho_dataset_species_level_tip_order$fresh_weight_SD, 
                            rownames(morpho_dataset_species_level_tip_order)) 

## Plot tree + barplot
plotTree.barplot(phy_2, fresh_weight,
                 args.plotTree = list(ftype = "i", fsize = 0.8, lwd = 0.9),  # Add species names (ftype = "off" for removing species names)
                 args.barplot = list(xlab = "body mass [mg]", space = 1, border = 'grey30'))

## Get barplot positions for adding error bars
bar_x <- fresh_weight
bar_y <- 1:length(fresh_weight)  # Vertical positions of bars
## Add error bars
arrows(x0 = bar_x - fresh_weight_sd, x1 = bar_x + fresh_weight_sd, 
       y0 = bar_y, y1 = bar_y, code = 3, angle = 90, length = 0.02, col = "black", lwd = 1)
#### -----> TO BE ADJUSTED IN ILLUSTRATOR DIRECTLY
dev.off()





############  ________ wing length vs mass: ############ 

## body mass scales with (wing length)^3 <=> wing length ~ body mass ^1/3 (because volume is length ^3)
## for wing length to scale isometrically with size: wing length ~ body mass ^(1/3)

## first plot the museum specimen only
plot(log10(morpho_dataset_species_level_tip_order$length_cm) ~ log10(morpho_dataset_species_level_tip_order$fresh_weight), pch=16, col=colo_species_filmed_invisible_tip_order, cex=2.7, bty='n', xlab='log(body mass)',ylab='log(wing length)')
## superimpose the filmed specimens:
points(log10(morpho_dataset_species_level_tip_order$length_cm)[filmed_tip_order] ~ log10(morpho_dataset_species_level_tip_order$fresh_weight)[filmed_tip_order], pch=21, col='black', bg=tip_cols[filmed_tip_order], cex=2.7)

## log scale but with actual values:
plot(morpho_dataset_species_level_tip_order$length_cm ~ morpho_dataset_species_level_tip_order$fresh_weight, pch=16, col=colo_species_filmed_invisible_tip_order, cex=2.7, bty='n', xlab='body mass [mg]',ylab='wing length [cm]', log='xy', axes=F)
points(morpho_dataset_species_level_tip_order$length_cm[filmed_tip_order] ~ morpho_dataset_species_level_tip_order$fresh_weight[filmed_tip_order], pch=21, col='black', bg=tip_cols[filmed_tip_order], cex=2.7)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)




############  _____________ PGLS ############ 

## PGLS using nlme (for some reason, the pgls using caper package doesn't work with hoverfly data)
## (we follow  https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(morpho_dataset_species_level_tip_order$fresh_weight); names(x_now) = phy$tip.label
y_now = log10(morpho_dataset_species_level_tip_order$length_cm); names(y_now) = phy$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy), method = "ML")
summary(pglsModel)

## Obtain the adjusted R-squared (not given through gls() function so we compute it ourselves)  
residuals <- residuals(pglsModel) # extract residuals
RSS <- sum(residuals^2) # compute RSS (Residual Sum of Squares)
## Compute the Total Sum of Square (TSS)
fitted_values <- pglsModel$fitted # extract fitted values
response_values <- fitted_values + residuals # Reconstruct the response variable
TSS <- sum((response_values - mean(response_values))^2) # Compute Total Sum of Squares (TSS)
# Get number of observations (n) and predictors (k)
n <- length(residuals) # Get number of observations (n)
k <- length(coef(pglsModel)) - 1  # Get number of predictors (k)
## Compute Adjusted R-squared
R2_adj <- 1 - ((RSS / (n - k - 1)) / (TSS / (n - 1))) ; R2_adj
# summary(lm(log10(morpho_dataset_species_level_tip_order$length_cm) ~ log10(morpho_dataset_species_level_tip_order$fresh_weight)))

## Obtain IC 95% of the pgls regression slope 
coef_est <- coef(summary(pglsModel))["x_now", "Value"] # extract coefficient estimates
se_est <- coef(summary(pglsModel))["x_now", "Std.Error"]  # extract standard errors
df <- summary(pglsModel)$dims$N - summary(pglsModel)$dims$p  # get residual degrees of freedom (df) (nb of observation - nb of estimated parameters)
## Compute t-value for 95% CI (alpha = 0.05, two-tailed)
t_value <- qt(0.975, df)  # 97.5th percentile for two-tailed test
## Compute confidence interval
CI_lower <- coef_est - t_value * se_est
CI_upper <- coef_est + t_value * se_est
## Display the Confidence Interval of the slope
c(CI_lower, CI_upper)



#### pgls (caper package) ---------> DOES NOT WORK !
# d = as.data.frame(matrix(NA, nrow(morpho_dataset_species_level_tip_order), 3 ))
# d[,1] = as.factor(phy$tip.label)
# d[,2] = log10(morpho_dataset_species_level_tip_order$fresh_weight)
# d[,3] = log10(morpho_dataset_species_level_tip_order$length_cm)
# colnames(d) = c('species', 'weight', 'wing_length')
# ## create a 'comparative data' object including the phylogeny, the data matching the tips, and the covariance matrix (with vcv=T)
# comp_data = comparative.data(phy, d, species, vcv=T)
# p = pgls(wing_length ~ weight, comp_data, lambda = 'ML') ; p ; anova(p)
# pgls_model <- summary(p) ; pgls_model
# ## get C.I 95% for the pgls slope:
# slope <- pgls_model$coefficients[2, 1]  # Coefficient (slope of wing_length)
# slope_se <- pgls_model$coefficients[2, 2]  # Standard error of the slope
# n <- nrow(morpho_dataset) # get sample size
# ## Compute the critical t-value for a 95% confidence interval
# t_crit <- qt(0.975, df = n - 2)  # Degrees of freedom = n - 2 for simple regression
# ## Calculate the 95% confidence interval for the slope
# lower_bound <- slope - t_crit * slope_se
# upper_bound <- slope + t_crit * slope_se
# ## Display the confidence interval for the slope
# c(lower_bound, upper_bound)


##### plot in log scale while keeping real values
plot(morpho_dataset_species_level_tip_order$length_cm ~ morpho_dataset_species_level_tip_order$fresh_weight, pch=16, col=colo_species_filmed_invisible_tip_order, cex=2.7, bty='n', xlab='body mass [mg]',ylab='wing length [cm]', log='xy', axes=F)
points(morpho_dataset_species_level_tip_order$length_cm[filmed_tip_order] ~ morpho_dataset_species_level_tip_order$fresh_weight[filmed_tip_order], pch=21, col='black', bg=tip_cols[filmed_tip_order], cex=2.7)
## calculate the isometric intercept to make the line pass by the mean value of the y variable
## this can be done as intercept = mean_y - slope * mean_x
iso_intercept <- mean(log10(morpho_dataset_species_level_tip_order$length_cm)) - 1/3 * mean(log10(morpho_dataset_species_level_tip_order$fresh_weight))
abline(iso_intercept, (1/3), col = '#ba4c36', lwd=4.5, lty='dashed') # isometric
abline(-0.4792984, 0.2552412, col='black', lwd=4.5) # PGLS regression
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)
intercept_W <- mean(log10(morpho_dataset_species_level_tip_order$length_cm)) - 2/9 * mean(log10(morpho_dataset_species_level_tip_order$fresh_weight))
abline(a = intercept_W, b = 2/9, col = 'grey30', lwd=6, lty='dotted') # required for weight support if everything else scales isometrically




############  _____________ phyRMA ############ 


####### Phylogenetic Reduced Major Axis (phyloRMA):
x_now = log10(morpho_dataset_species_level_tip_order$fresh_weight); names(x_now) = phy$tip.label
y_now = log10(morpho_dataset_species_level_tip_order$length_cm); names(y_now) = phy$tip.label
rma_result <- phyl.RMA(x_now, y_now, phy, fixed=T) ; rma_result

### test is it deviate from isometry based on function-embedded test
phyl.RMA(x_now, y_now, phy, fixed=T, h0=1/3)


##### Estimate the Confidence Interval using phylogenetic bootstrap residual 
# # Extract residuals
# predicted_y <- rma_result$RMA.beta[1] + rma_result$RMA.beta[2] * x_now
# residuals <- y_now - predicted_y
# 
# # Bootstrap function
# set.seed(42)  # For reproducibility
# num_boot <- 1000  # Number of bootstrap iterations
# boot_slopes <- numeric(num_boot)
# 
# for (i in 1:num_boot) {
#   # Resample residuals with replacement
#   boot_residuals <- sample(residuals, replace = TRUE)
#   # Construct new y-values by adding bootstrapped residuals to predicted values
#   y_boot <- predicted_y + boot_residuals
#   # Fit new RMA model
#   boot_rma <- phyl.RMA(x_now, y_boot, phy, fixed = TRUE)
#   # Store the bootstrapped slope
#   boot_slopes[i] <- boot_rma$RMA.beta[2]
# }
# 
# # Step 4: Compute 95% Confidence Interval
# ci_lower <- quantile(boot_slopes, 0.025)
# ci_upper <- quantile(boot_slopes, 0.975)
# ## Display the Confidence Interval of the slope
# c(ci_lower, ci_upper)


##### Estimate the Confidence Interval using the Standard Error (CI = slope - 1.96 * slope_se)
slope <- rma_result$RMA.beta[2] # Extract slope estimate
slope_var <- rma_result$V[2,2] # Extract variance of the slope from the VCV matrix
slope_se <- sqrt(slope_var) # Compute standard error (SE)
# Approximate 95% confidence interval using Normal distribution
ci_lower <- slope - 1.96 * slope_se
ci_upper <- slope + 1.96 * slope_se
## Display the Confidence Interval of the slope
c(ci_lower, ci_upper)
### this CI looks very wide... !



###################################### OLD 

# ## superimpose phylogeny:
# # pheno = data.frame(log10(morpho_dataset_species_level_tip_order$fresh_weight), log10(morpho_dataset_species_level_tip_order$length_cm)) ; row.names(pheno) <- phy$tip.label
# # plot(pheno, pch=16, lwd=3, col=tip_cols, cex=2, bty='n', xlab='log(body mass)',ylab='log(wing length)')
# # phylomorphospace(phy, pheno, lwd=1, node.size=c(0.7,0), label="off", fsize=1.2, add=T)
# # points(pheno, pch=16, lwd=3, col=tip_cols, cex=1.1)
# 
# ## classic Reduced Major Axis regression (RMA):
# lmodel2(log10(morpho_dataset_species_level_tip_order$length_cm) ~ log10(morpho_dataset_species_level_tip_order$fresh_weight))
# 
# ## phylogenetic RMA:
# x_now = log10(morpho_dataset_species_level_tip_order$fresh_weight); names(x_now) = phy$tip.label
# y_now = log10(morpho_dataset_species_level_tip_order$length_cm); names(y_now) = phy$tip.label
# phyl.RMA(x_now, y_now, phy, h0=1/3)
# 
# abline(-0.5041353, 0.2844828  , col='black', lwd=2)         # observed (MA row of the lmodel2() output)
# # abline(-0.5154928, 0.2972948   , col='darkred', lwd=2)    # phylogenetic RMA (phyl.RMA output)
# ## calculate the intercept such as the line pass by the mean value of the y variable
# ## this can be done as intercept = mean_y - slope * mean_x
# iso_intercept <- mean(log10(morpho_dataset_species_level_tip_order$length_cm)) - 1/3 * mean(log10(morpho_dataset_species_level_tip_order$fresh_weight))
# abline(iso_intercept,(1/3), col='grey35', lwd=2, lty='dashed') # isometric
# legend('topleft', c('observed scaling (0.28)', 'isometric scaling (0.33)'), col=c('black', 'grey35'), bty='n',lty=c('solid','dashed'),cex=0.7,pt.cex=1.2)
# 
# ## test whether the two slopes differ significantly:
# sma(log10(morpho_dataset_species_level_tip_order$length_cm) ~ log10(morpho_dataset_species_level_tip_order$fresh_weight), slope.test = 1/3)
# ## NON SIGNIFICANT (but CI almost exclude 0.33)




############  ________ wing chord vs mass: ############ 

plot(log10(morpho_dataset_species_level_tip_order$average_chord) ~ log10(morpho_dataset_species_level_tip_order$fresh_weight), pch=16, col=colo_species_filmed_invisible_tip_order, cex=2, bty='n', xlab='log(body mass)',ylab='log(wing chord)')
points(log10(morpho_dataset_species_level_tip_order$average_chord)[filmed_tip_order] ~ log10(morpho_dataset_species_level_tip_order$fresh_weight)[filmed_tip_order], pch=21, col='black', bg=tip_cols[filmed_tip_order], cex=2.7)

## log scale but with actual values:
plot(morpho_dataset_species_level_tip_order$average_chord ~ morpho_dataset_species_level_tip_order$fresh_weight, pch=16, col=colo_species_filmed_invisible_tip_order, cex=2, bty='n', xlab='body mass [mg]',ylab='wing chord [cm]', log='xy', axes=F)
points(morpho_dataset_species_level_tip_order$average_chord[filmed_tip_order] ~ morpho_dataset_species_level_tip_order$fresh_weight[filmed_tip_order], pch=21, col='black', bg=tip_cols[filmed_tip_order], cex=2.7)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)


############  _____________ PGLS ############ 


## PGLS using nlme (for some reason, the pgls using caper package doesn't work with hoverfly data)
## (we follow  https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(morpho_dataset_species_level_tip_order$fresh_weight); names(x_now) = phy$tip.label
y_now = log10(morpho_dataset_species_level_tip_order$average_chord); names(y_now) = phy$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy), method = "ML")
summary(pglsModel)

## Obtain the adjusted R-squared (not given through gls() function so we compute it ourselves)  
residuals <- residuals(pglsModel) # extract residuals
RSS <- sum(residuals^2) # compute RSS (Residual Sum of Squares)
## Compute the Total Sum of Square (TSS)
fitted_values <- pglsModel$fitted # extract fitted values
response_values <- fitted_values + residuals # Reconstruct the response variable
TSS <- sum((response_values - mean(response_values))^2) # Compute Total Sum of Squares (TSS)
# Get number of observations (n) and predictors (k)
n <- length(residuals) # Get number of observations (n)
k <- length(coef(pglsModel)) - 1  # Get number of predictors (k)
## Compute Adjusted R-squared
R2_adj <- 1 - ((RSS / (n - k - 1)) / (TSS / (n - 1))) ; R2_adj

## Obtain IC 95% of the pgls regression slope 
coef_est <- coef(summary(pglsModel))["x_now", "Value"] # extract coefficient estimates
se_est <- coef(summary(pglsModel))["x_now", "Std.Error"]  # extract standard errors
df <- summary(pglsModel)$dims$N - summary(pglsModel)$dims$p  # get residual degrees of freedom (df) (nb of observation - nb of estimated parameters)
## Compute t-value for 95% CI (alpha = 0.05, two-tailed)
t_value <- qt(0.975, df)  # 97.5th percentile for two-tailed test
## Compute confidence interval
CI_lower <- coef_est - t_value * se_est
CI_upper <- coef_est + t_value * se_est
## Display the Confidence Interval of the slope
c(CI_lower, CI_upper)


##### plot in log scale while keeping real values
plot(morpho_dataset_species_level_tip_order$average_chord ~ morpho_dataset_species_level_tip_order$fresh_weight, pch=16, col=colo_species_filmed_invisible_tip_order, cex=2.7, bty='n', xlab='body mass [mg]',ylab='wing chord [cm]', log='xy', axes=F)
points(morpho_dataset_species_level_tip_order$average_chord[filmed_tip_order] ~ morpho_dataset_species_level_tip_order$fresh_weight[filmed_tip_order], pch=21, col='black', bg=tip_cols[filmed_tip_order], cex=2.7)
## calculate the isometric intercept to make the line pass by the mean value of the y variable
## this can be done as intercept = mean_y - slope * mean_x
iso_intercept <- mean(log10(morpho_dataset_species_level_tip_order$average_chord)) - 1/3 * mean(log10(morpho_dataset_species_level_tip_order$fresh_weight))
abline(iso_intercept, (1/3), col = '#ba4c36', lwd=4.5, lty='dashed') # isometric
abline(-1.1449591, 0.2944177, col='black', lwd=4.5) # PGLS regression
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)
intercept_W <- mean(log10(morpho_dataset_species_level_tip_order$average_chord)) - 0 * mean(log10(morpho_dataset_species_level_tip_order$fresh_weight))
abline(a = intercept_W, b = 0, col = 'grey30', lwd=6, lty='dotted') # required for weight support if everything else scales isometrically




################## OLD

# ## classic Reduced Major Axis regression (RMA):
# lmodel2(log10(morpho_dataset_species_level_tip_order$average_chord) ~ log10(morpho_dataset_species_level_tip_order$fresh_weight)) ## check the MA line in the output
# 
# ## phylogenetic RMA:
# x_now = log10(morpho_dataset_species_level_tip_order$fresh_weight); names(x_now) = phy$tip.label
# y_now = log10(morpho_dataset_species_level_tip_order$average_chord); names(y_now) = phy$tip.label
# phyl.RMA(x_now, y_now, phy, h0=1/3)
# 
# ## because "wing chord ~ body mass ^(1/3)", the slope should be 1/3 if chord scales isometrically with mass
# abline(-1.176880 , 0.3343728  , col='black', lwd=2)          # observed (MA line)
# # abline(-1.1885144, 0.3349356   , col='darkred', lwd=2)         # phylogenetic RMA (phyl.RMA output)
# iso_intercept <- mean(log10(morpho_dataset_species_level_tip_order$average_chord)) - 1/3 * mean(log10(morpho_dataset_species_level_tip_order$fresh_weight))
# abline(iso_intercept ,(1/3), col='grey35', lwd=2, lty='dashed') # isometric
# legend('topleft', c('observed scaling (1.33)', 'isometric scaling (0.33)'), col=c('black', 'grey35'), bty='n',lty=c('solid','dashed'),cex=0.7, pt.cex=1.2)
# ## test whether the two slopes differ significantly
# sma(log10(morpho_dataset_species_level_tip_order$average_chord) ~ log10(morpho_dataset_species_level_tip_order$fresh_weight), slope.test = 1/3)
# ## NON SIGNIFICANT



############  _____________ phyRMA ############ 

####### Phylogenetic Reduced Major Axis (phyloRMA):
x_now = log10(morpho_dataset_species_level_tip_order$fresh_weight); names(x_now) = phy$tip.label
y_now = log10(morpho_dataset_species_level_tip_order$average_chord); names(y_now) = phy$tip.label
rma_result <- phyl.RMA(x_now, y_now, phy, fixed=T) ; rma_result

### test is it deviate from isometry based on function-embedded test
phyl.RMA(x_now, y_now, phy, fixed=T, h0=1/3)







############  ________ S2 vs mass: ############ 


plot(log10(morpho_dataset_species_level_tip_order$secondMoment) ~ log10(morpho_dataset_species_level_tip_order$fresh_weight), pch=16, col=colo_species_filmed_invisible_tip_order, cex=2.7, bty='n', xlab='log(body mass)',ylab='log(second moment of area)')
points(log10(morpho_dataset_species_level_tip_order$secondMoment)[filmed_tip_order] ~ log10(morpho_dataset_species_level_tip_order$fresh_weight)[filmed_tip_order], pch=21, col='black', bg=tip_cols[filmed_tip_order], cex=2.7)

## log scale but with actual values:
plot(morpho_dataset_species_level_tip_order$secondMoment ~ morpho_dataset_species_level_tip_order$fresh_weight, pch=16, col=colo_species_filmed_invisible_tip_order, cex=2.7, bty='n', xlab='body mass [mg]',ylab='second moment of area [cm^4]', log='xy', axes=F)
points(morpho_dataset_species_level_tip_order$secondMoment[filmed_tip_order] ~ morpho_dataset_species_level_tip_order$fresh_weight[filmed_tip_order], pch=21, col='black', bg=tip_cols[filmed_tip_order], cex=2.7)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)


############  _____________ PGLS ############ 


## PGLS using nlme (for some reason, the pgls using caper package doesn't work with hoverfly data)
## (we follow  https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(morpho_dataset_species_level_tip_order$fresh_weight); names(x_now) = phy$tip.label
y_now = log10(morpho_dataset_species_level_tip_order$secondMoment); names(y_now) = phy$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy), method = "ML")
summary(pglsModel)

## Obtain the adjusted R-squared (not given through gls() function so we compute it ourselves)  
residuals <- residuals(pglsModel) # extract residuals
RSS <- sum(residuals^2) # compute RSS (Residual Sum of Squares)
## Compute the Total Sum of Square (TSS)
fitted_values <- pglsModel$fitted # extract fitted values
response_values <- fitted_values + residuals # Reconstruct the response variable
TSS <- sum((response_values - mean(response_values))^2) # Compute Total Sum of Squares (TSS)
# Get number of observations (n) and predictors (k)
n <- length(residuals) # Get number of observations (n)
k <- length(coef(pglsModel)) - 1  # Get number of predictors (k)
## Compute Adjusted R-squared
R2_adj <- 1 - ((RSS / (n - k - 1)) / (TSS / (n - 1))) ; R2_adj

## Obtain IC 95% of the pgls regression slope 
coef_est <- coef(summary(pglsModel))["x_now", "Value"] # extract coefficient estimates
se_est <- coef(summary(pglsModel))["x_now", "Std.Error"]  # extract standard errors
df <- summary(pglsModel)$dims$N - summary(pglsModel)$dims$p  # get residual degrees of freedom (df) (nb of observation - nb of estimated parameters)
## Compute t-value for 95% CI (alpha = 0.05, two-tailed)
t_value <- qt(0.975, df)  # 97.5th percentile for two-tailed test
## Compute confidence interval
CI_lower <- coef_est - t_value * se_est
CI_upper <- coef_est + t_value * se_est
## Display the Confidence Interval of the slope
c(CI_lower, CI_upper)


##### plot in log scale while keeping real values
plot(morpho_dataset_species_level_tip_order$secondMoment ~ morpho_dataset_species_level_tip_order$fresh_weight, pch=16, col=colo_species_filmed_invisible_tip_order, cex=2.7, bty='n', xlab='body mass [mg]',ylab='second moment of area [cm^4]', log='xy', axes=F)
points(morpho_dataset_species_level_tip_order$secondMoment[filmed_tip_order] ~ morpho_dataset_species_level_tip_order$fresh_weight[filmed_tip_order], pch=21, col='black', bg=tip_cols[filmed_tip_order], cex=2.7)
## calculate the isometric intercept to make the line pass by the mean value of the y variable
## this can be done as intercept = mean_y - slope * mean_x
iso_intercept <- mean(log10(morpho_dataset_species_level_tip_order$secondMoment)) - 4/3 * mean(log10(morpho_dataset_species_level_tip_order$fresh_weight))
abline(iso_intercept, (4/3), col = '#ba4c36', lwd=4.5, lty='dashed') # isometric
abline(-2.949733, 1.008725, col='black', lwd=4.5) # PGLS regression
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)
intercept_W <- mean(log10(morpho_dataset_species_level_tip_order$secondMoment)) - 1 * mean(log10(morpho_dataset_species_level_tip_order$fresh_weight))
abline(a = intercept_W, b = 1, col = 'grey30', lwd=6, lty='dotted') # required for weight support if everything else scales isometrically



####### Phylogenetic Reduced Major Axis (phyloRMA):
x_now = log10(morpho_dataset_species_level_tip_order$fresh_weight); names(x_now) = phy$tip.label
y_now = log10(morpho_dataset_species_level_tip_order$secondMoment); names(y_now) = phy$tip.label
rma_result <- phyl.RMA(x_now, y_now, phy, fixed=T) ; rma_result

### test is it deviate from isometry based on function-embedded test
phyl.RMA(x_now, y_now, phy, fixed=T, h0=4/3)








############  ________ non-dim S2 vs mass: ############ 

plot(log10(morpho_dataset_species_level_tip_order$non_dimensional_2ndMoment) ~ log10(morpho_dataset_species_level_tip_order$fresh_weight), pch=16, col=colo_species_filmed_invisible_tip_order, cex=2.7, bty='n', xlab='log(body mass)',ylab='log(second moment of area)')
points(log10(morpho_dataset_species_level_tip_order$non_dimensional_2ndMoment)[filmed_tip_order] ~ log10(morpho_dataset_species_level_tip_order$fresh_weight)[filmed_tip_order], pch=21, col='black', bg=tip_cols[filmed_tip_order], cex=2.7)

## log scale but with actual values:
plot(morpho_dataset_species_level_tip_order$non_dimensional_2ndMoment ~ morpho_dataset_species_level_tip_order$fresh_weight, pch=16, col=colo_species_filmed_invisible_tip_order, cex=2.7, bty='n', xlab='body mass [mg]',ylab='dimensionless S2 [-]', log='xy', axes=F)
points(morpho_dataset_species_level_tip_order$non_dimensional_2ndMoment[filmed_tip_order] ~ morpho_dataset_species_level_tip_order$fresh_weight[filmed_tip_order], pch=21, col='black', bg=tip_cols[filmed_tip_order], cex=2.7)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)



############  _____________ PGLS ############ 


## PGLS using nlme (for some reason, the pgls using caper package doesn't work with hoverfly data)
## (we follow  https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(morpho_dataset_species_level_tip_order$fresh_weight); names(x_now) = phy$tip.label
y_now = log10(morpho_dataset_species_level_tip_order$non_dimensional_2ndMoment); names(y_now) = phy$tip.label
pglsModel <- gls(y_now ~ x_now, correlation = corBrownian(phy = phy), method = "ML")
summary(pglsModel)

## Obtain the adjusted R-squared (not given through gls() function so we compute it ourselves)  
residuals <- residuals(pglsModel) # extract residuals
RSS <- sum(residuals^2) # compute RSS (Residual Sum of Squares)
## Compute the Total Sum of Square (TSS)
fitted_values <- pglsModel$fitted # extract fitted values
response_values <- fitted_values + residuals # Reconstruct the response variable
TSS <- sum((response_values - mean(response_values))^2) # Compute Total Sum of Squares (TSS)
# Get number of observations (n) and predictors (k)
n <- length(residuals) # Get number of observations (n)
k <- length(coef(pglsModel)) - 1  # Get number of predictors (k)
## Compute Adjusted R-squared
R2_adj <- 1 - ((RSS / (n - k - 1)) / (TSS / (n - 1))) ; R2_adj

## Obtain IC 95% of the pgls regression slope 
coef_est <- coef(summary(pglsModel))["x_now", "Value"] # extract coefficient estimates
se_est <- coef(summary(pglsModel))["x_now", "Std.Error"]  # extract standard errors
df <- summary(pglsModel)$dims$N - summary(pglsModel)$dims$p  # get residual degrees of freedom (df) (nb of observation - nb of estimated parameters)
## Compute t-value for 95% CI (alpha = 0.05, two-tailed)
t_value <- qt(0.975, df)  # 97.5th percentile for two-tailed test
## Compute confidence interval
CI_lower <- coef_est - t_value * se_est
CI_upper <- coef_est + t_value * se_est
## Display the Confidence Interval of the slope
c(CI_lower, CI_upper)


##### plot in log scale while keeping real values
plot(morpho_dataset_species_level_tip_order$non_dimensional_2ndMoment ~ morpho_dataset_species_level_tip_order$fresh_weight, pch=16, col=colo_species_filmed_invisible_tip_order, cex=2.7, bty='n', xlab='body mass [mg]',ylab='dimensionless S2 [-]', log='xy', axes=F)
points(morpho_dataset_species_level_tip_order$non_dimensional_2ndMoment[filmed_tip_order] ~ morpho_dataset_species_level_tip_order$fresh_weight[filmed_tip_order], pch=21, col='black', bg=tip_cols[filmed_tip_order], cex=2.7)
## calculate the isometric intercept to make the line pass by the mean value of the y variable
## this can be done as intercept = mean_y - slope * mean_x
iso_intercept <- mean(log10(morpho_dataset_species_level_tip_order$non_dimensional_2ndMoment)) - 0 * mean(log10(morpho_dataset_species_level_tip_order$fresh_weight))
abline(iso_intercept, (0), col = '#ba4c36', lwd=4.5, lty='dashed') # isometric
abline(-0.4135408, -0.0347662, col='black', lwd=4.5) # PGLS regression
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
# axis(2, at=log10Tck('x','minor'), tcl= 0.2)
axis(2)
intercept_W <- mean(log10(morpho_dataset_species_level_tip_order$non_dimensional_2ndMoment)) - (-1/3) * mean(log10(morpho_dataset_species_level_tip_order$fresh_weight))
abline(a = intercept_W, b = (-1/3), col = 'grey30', lwd=6, lty='dotted') # required for weight support if everything else scales isometrically


################### OLD

# ## classic Reduced Major Axis regression (RMA):
# lmodel2(log10(morpho_dataset_species_level_tip_order$non_dimensional_2ndMoment) ~ log10(morpho_dataset_species_level_tip_order$fresh_weight)) ## check the MA line in the output
# 
# ## phylogenetic RMA: TEST AGAINST slope 0 IS NOT ALLOWED
# # x_now = log10(morpho_dataset_species_level_tip_order$fresh_weight); names(x_now) = phy$tip.label
# # y_now = log10(morpho_dataset_species_level_tip_order$non_dimensional_2ndMoment); names(y_now) = phy$tip.label
# # phyl.RMA(x_now, y_now, phy, h0=0)
# 
# ## because nd_S2 ~ body mass ^(0)", the slope should be 0 if S2 scales isometrically with mass
# abline(-0.2090921 , -0.01983373  , col='black', lwd=2)          # observed (MA line)
# # abline(-0.19565829, -0.02652135   , col='darkred', lwd=2)         # phylogenetic RMA (phyl.RMA output)
# iso_intercept <- mean(log10(morpho_dataset_species_level_tip_order$non_dimensional_2ndMoment)) - 0 * mean(log10(morpho_dataset_species_level_tip_order$fresh_weight))
# abline(iso_intercept ,0, col='grey35', lwd=2, lty='dashed') # isometric
# legend('topright', c('observed scaling (1.33)', 'isometric scaling (-0.02)'), col=c('black', 'grey35'), bty='n',lty=c('solid','dashed'),cex=0.7, pt.cex=1.2)
# ## test whether the two slopes differ significantly
# sma(log10(morpho_dataset_species_level_tip_order$secondMoment) ~ log10(morpho_dataset_species_level_tip_order$fresh_weight), slope.test = 0)
# ## SIGNIFICANT




############  _____________ phyRMA ############ 

####### Phylogenetic Reduced Major Axis (phyloRMA):
x_now = log10(morpho_dataset_species_level_tip_order$fresh_weight); names(x_now) = phy$tip.label
y_now = log10(morpho_dataset_species_level_tip_order$non_dimensional_2ndMoment); names(y_now) = phy$tip.label
rma_result <- phyl.RMA(x_now, y_now, phy, fixed=T) ; rma_result

### test is it deviate from isometry based on function-embedded test
phyl.RMA(x_now, y_now, phy, fixed=T, h0=0.1)





