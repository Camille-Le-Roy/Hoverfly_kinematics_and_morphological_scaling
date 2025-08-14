





rm(list=ls())



library(scico)
library(colorspace)
library(tagcloud)
library(RColorBrewer)
library(lmodel2)
library(pracma)
library(nlme)
library(phytools)
library(ggplot2)

local_path = dirname(rstudioapi::getSourceEditorContext()$path)


######################## [1] preamble ######################## 

## Thomas Engels performed CFD simulation using the wingbeat kinematics of hoverflies (averaged across all species)
## 3 series of simulations were performed:
## 
## Series 1: wingbeat frequency and wing length are kept the same, species-specific wing shape 
## --> wing shape effect only on force production
## Series 2: wingbeat frequency kept the same, species-specific wing shape and size 
## --> wing shape and size effect on force production
## Series 3: species-specific wingbeat frequency 
## --> wingbeat frequency, wing shape and size effect on force production
## (comparing Series 2 vs Series 3 allow assessing the effect of Reynolds number)




######################## [2] formatting data ######################## 

full_dataset <- read.csv(paste0(substr(local_path, 1, as.numeric(gregexpr(pattern ='/04_CFD Thomas Engels',local_path)[[1]])), '/full_dataset_species_level.csv'))
## data generated in the code "BBhov_analyzing_WBkinematics_vs_morphology.R"


############ ___ (i) CFD data from Series 1 ############ 

#### load CFD data from Series 1
CFD_data_Series1 <- array(NA, dim = c(400, 9, nrow(full_dataset)))
dimnames(CFD_data_Series1)[[3]] <- full_dataset$genus_species
col_names <-  colnames(read.csv(paste0(local_path,'/data received from T. Engels/time data series 1/time_data_Episyrphus_viridaureus.csv'), sep=';'))
dimnames(CFD_data_Series1)[[2]] <- c(col_names, 'total_force')
path_to_files = paste0(local_path,'/data received from T. Engels/time data series 1/')
for (i in 1:dim(CFD_data_Series1)[3]){
  file_name_now <- list.files(paste0(path_to_files), pattern=".csv")[i]
  data_now <- read.csv(paste0(path_to_files,file_name_now), sep = ';')
  CFD_data_Series1[,1:ncol(data_now),i] <- as.matrix(data_now)
}
## compute total force as sqrt(Fz^2 + Fx^2)
for (i in 1:dim(CFD_data_Series1)[3]){
  for (j in 1:dim(CFD_data_Series1)[1]){
    # CFD_data_Series1[j,'total_force', i] = sqrt( CFD_data_Series1[j,'fz', i]^2 + CFD_data_Series1[j,'fx', i]^2 )
    CFD_data_Series1[j,'total_force', i] = CFD_data_Series1[j,'fz', i] }} # just take Fz
# View(CFD_data_Series1[,,1])

## the output from CFD simulations performed by Thomas is dimensionless.
## to obtain the dimensional forces, we multiply F with rho_air * R^4 * f^2
CFD_data_Series1_dimensional <- CFD_data_Series1[,c(1,9),] # we just select column 'time' and 'total_force'
dimnames(CFD_data_Series1_dimensional)[[3]] <- full_dataset$genus_species
for (i in 1:dim(CFD_data_Series1_dimensional)[3]){
  for (j in 1:dim(CFD_data_Series1_dimensional)[1]){
    F_now = CFD_data_Series1_dimensional[j,'total_force', i]
    CFD_data_Series1_dimensional[j,'total_force', i] = F_now * 1.225 * (mean(full_dataset$length_cm)/100)^4 * (mean(full_dataset$WB_freq))^2
  }} # we use the mean wing length and the mean WB frequency because those stay constant in series 1
# View(CFD_data_Series1_dimensional[,,1])



############ ___ (ii) CFD data from Series 2 ############ 

#### load CFD data from Series 2
CFD_data_Series2 <- array(NA, dim = c(400, 9, nrow(full_dataset)))
dimnames(CFD_data_Series2)[[3]] <- full_dataset$genus_species
col_names <-  colnames(read.csv(paste0(local_path,'/data received from T. Engels/time data series 2/time_data_Episyrphus_viridaureus_series2_series2.csv'), sep=';'))
dimnames(CFD_data_Series2)[[2]] <- c(col_names, 'total_force')

path_to_files = paste0(local_path,'/data received from T. Engels/time data series 2/')
for (i in 1:dim(CFD_data_Series2)[3]){
  file_name_now <- list.files(paste0(path_to_files), pattern=".csv")[i]
  data_now <- read.csv(paste0(path_to_files,file_name_now), sep = ';')
  CFD_data_Series2[,1:ncol(data_now),i] <- as.matrix(data_now)
}
## compute total force as sqrt(Fz^2 + Fx^2)
for (i in 1:dim(CFD_data_Series2)[3]){
  for (j in 1:dim(CFD_data_Series2)[1]){
    # CFD_data_Series2[j,'total_force', i] = sqrt( CFD_data_Series2[j,'fz', i]^2 + CFD_data_Series2[j,'fx', i]^2 )
    CFD_data_Series2[j,'total_force', i] = CFD_data_Series2[j,'fz', i] }} # just take Fz
# View(CFD_data_Series2[,,1])

## the output from CFD simulations performed by Thomas is dimensionless.
## to obtain the dimensional forces, we multiply F with rho_air * R^4 * f^2
CFD_data_Series2_dimensional <- CFD_data_Series1[,c(1,9),] # we just select column 'time' and 'total_force'
dimnames(CFD_data_Series2_dimensional)[[3]] <- full_dataset$genus_species
for (i in 1:dim(CFD_data_Series2_dimensional)[3]){
  for (j in 1:dim(CFD_data_Series2_dimensional)[1]){
    F_now = CFD_data_Series2_dimensional[j,'total_force', i]
    CFD_data_Series2_dimensional[j,'total_force', i] = F_now * 1.225 * (full_dataset$length_cm[i]/100)^4 * (mean(full_dataset$WB_freq))^2
  }} # we use the species specific wing length and the mean WB frequency (Series 2 configuration)
# View(CFD_data_Series2_dimensional[,,1])



############ ___ (iii) CFD data from Series 3 ############ 

## the Series 3 has actually not been performed
## we do not have time series but we have the average force per species.
## load the output data from series 3, where variation in WB frequency was included:
dimensional_data_series3 <- read.csv(paste0(local_path,'/data received from T. Engels/data_series3_dimensional_values.csv'))
dimensional_data_series3$total_force_N <- rep(NA, nrow(dimensional_data_series3))
for (i in 1:nrow(dimensional_data_series3)){
  dimensional_data_series3$total_force_N[i] = sqrt( dimensional_data_series3$Fz_N[i]^2 + dimensional_data_series3$Fx_N[i]^2 )
}

## we can however compute the time series for Series 3 ourselves
## by doing: F from Series 2 (using mean freq) * (f^2 of the current species / f^2 average)
CFD_data_Series3_dimensional <- CFD_data_Series2_dimensional
for (i in 1:dim(CFD_data_Series3_dimensional)[3]){
  for (j in 1:dim(CFD_data_Series3_dimensional)[1]){
    F_now = CFD_data_Series3_dimensional[j,'total_force', i]
    CFD_data_Series3_dimensional[j,'total_force', i] = F_now * (full_dataset$WB_freq[i]^2 / (mean(full_dataset$WB_freq))^2 )
  }}




########### ___ (vi) extract average forces for each Series ########### 

for (i in 1:nrow(full_dataset)){
  ## Series 1
  full_dataset$force_nondim_Series1[i] = mean(CFD_data_Series1[,'total_force',i])
  full_dataset$force_dimens_Series1[i] = mean(CFD_data_Series1_dimensional[,'total_force',i])
  
  ## Series 2
  full_dataset$force_nondim_Series2[i] = mean(CFD_data_Series2[,'total_force',i])
  full_dataset$force_dimens_Series2[i] = mean(CFD_data_Series2_dimensional[,'total_force',i])
  
  ## Series 3
  # full_dataset$force_dimens_Series3[i] = dimensional_data_series3$total_force_N[i]
  # full_dataset$force_dimens_Series3[i] = mean(CFD_data_Series3_dimensional[,'total_force',i])
  full_dataset$force_dimens_Series3[i] = trapz(CFD_data_Series3_dimensional[,'time',i], CFD_data_Series3_dimensional[,'total_force',i])
}

## ** CDF simulations don't have a constant sampling frequency.
## Because of this we cannot simply take the mean of all the temporal dynamic of forces.
## We have to do a Trapezoidal Integration, i.e. compute the area under the curve, using trapz() from package pracma




########### ___ (v) color coding ########### 

full_dataset$genus_species <- as.factor(full_dataset$genus_species)
levels(full_dataset$genus_species)

colors_species = c('#D1BA15', # Episyrphus_viridaureus  ## darker (greener) color are heavier species here
                   '#2B3A21', # Eristalis_tenax
                   '#78C155', # Eupeodes_nielseni
                   '#A63D40', # Melanostoma_mellinum
                   '#7BCEEF', # Platycheirus_clypeatus
                   '#674ACE', # Sphaerophoria_scripta
                   '#287C34', # Syrphus_nigrilinearus
                   '#F07E15') # Tropidia_scita


####### set color according to weight:
# weight_color <- scico(25, palette = 'bamako')[seq(1,22,by=3)] #choose gradient among: scico_palette_show()
# weight_color = weight_color[length(weight_color):1] # reverse the gradient
weight_color <- brewer.pal(9,'BuPu')[4:9] # (i) choose gradient among: display.brewer.all()

# plot(rep(1,length(weight_color)), pch=15, cex=15,col=weight_color, bty='n', axes=F) # (ii) visually check the gradient
### set a gradient
variable_col_gradient = full_dataset$fresh_weight # (iii) associate the gradient to the values of a quantitative variable
gradient_now <- smoothPalette(variable_col_gradient, weight_color)
# plot(full_dataset$fresh_weight, as.factor(full_dataset$genus_species), col=gradient_now, pch=16, cex=4, xlab='body mass [mg]')



########## ____ *** save full_dataset ########## 

# write.csv(full_dataset, paste0(local_path, '/full_dataset_incl_mean_forces.csv'))



######################## [3] temporal dynamic of forces ######################## 

############ ___ (i) dimensionless force ############ 
## Series 1 dimensionless
p=1
plot((CFD_data_Series1[,'time', p]-2)*100, CFD_data_Series1[,'total_force', p], typ='l', lwd=2, ylim=c(0,4), col=colors_species[p], xlab=c('wingbeat cycle [%]'), ylab=c('dimensionless force'), bty='n', main='series 1')
for (p in 2:dim(CFD_data_Series1)[3]){points((CFD_data_Series1[,'time', p]-2)*100, CFD_data_Series1[,'total_force', p], typ='l', lwd=2, col=colors_species[p])}
## Series 2 dimensionless
p=1
plot((CFD_data_Series2[,'time', p]-2)*100, CFD_data_Series2[,'total_force', p], typ='l', lwd=2, ylim=c(0,4), col=colors_species[p], xlab=c('wingbeat cycle [%]'), ylab=c('dimensionless force'), bty='n', main='series 2')
for (p in 2:dim(CFD_data_Series2)[3]){points((CFD_data_Series2[,'time', p]-2)*100, CFD_data_Series2[,'total_force', p], typ='l', lwd=2, col=colors_species[p])}

####### save it for Xiaolei Mou:
### save Series 1 dimensionless
for (i in 1:dim(CFD_data_Series1)[3]) {
  slice_df <- as.data.frame(CFD_data_Series1[,,i])
  file_name <- paste0(dimnames(CFD_data_Series1)[[3]][i], ".csv")
  saving_location <- paste0(local_path, '/CFD data processed (for Xiaolei Mou)/CFD_Series_1_dimensionless/')
  write.csv(slice_df, file = paste0(saving_location,file_name), row.names = F)
}
### save Series 2 dimensionless
for (i in 1:dim(CFD_data_Series2)[3]) {
  slice_df <- as.data.frame(CFD_data_Series2[,,i])
  file_name <- paste0(dimnames(CFD_data_Series2)[[3]][i], ".csv")
  saving_location <- paste0(local_path, '/CFD data processed (for Xiaolei Mou)/CFD_Series_2_dimensionless/')
  write.csv(slice_df, file = paste0(saving_location,file_name), row.names = F)
}

############ ___ (ii) dimensional force ############ 
## Series 1 dimensional
p=1
plot((CFD_data_Series1_dimensional[,'time', p]-2)*100, CFD_data_Series1_dimensional[,'total_force', p]*2*1000, typ='l', lwd=2, ylim=c(0,2), col=colors_species[p], xlab=c('wingbeat cycle [%]'), ylab=c('force [mN]'), bty='n', main='series 1')
for (p in 2:dim(CFD_data_Series1_dimensional)[3]){points((CFD_data_Series1_dimensional[,'time', p]-2)*100, CFD_data_Series1_dimensional[,'total_force', p]*2*1000, typ='l', lwd=2, col=colors_species[p])}

## Series 2 dimensional
p=1
plot((CFD_data_Series2_dimensional[,'time', p]-2)*100, CFD_data_Series2_dimensional[,'total_force', p]*2*1000, typ='l', lwd=2, ylim=c(0,3), col=colors_species[p], xlab=c('wingbeat cycle [%]'), ylab=c('force [mN]'), bty='n', main='series 2')
for (p in 2:dim(CFD_data_Series2_dimensional)[3]){points((CFD_data_Series2_dimensional[,'time', p]-2)*100, CFD_data_Series2_dimensional[,'total_force', p]*2*1000, typ='l', lwd=2, col=colors_species[p])}
## Series 2 dimensional normalized by m*g (should fall on top of each other)
p=1
plot((CFD_data_Series2_dimensional[,'time', p]-2)*100, (CFD_data_Series2_dimensional[,'total_force', p]*2)/((full_dataset$fresh_weight[p]*10^-6)*9.81), typ='l', lwd=2, ylim=c(0,6), col=colors_species[p], xlab=c('wingbeat cycle [%]'), ylab=c('weight-normalized force F / mg [-]'), bty='n', main='series 2')
for (p in 2:dim(CFD_data_Series2_dimensional)[3]){
  points((CFD_data_Series2_dimensional[,'time', p]-2)*100, (CFD_data_Series2_dimensional[,'total_force', p]*2)/((full_dataset$fresh_weight[p]*10^-6)*9.81), typ='l', lwd=2, col=colors_species[p])
} ## curve are not falling on top of each other


###### save it for Xiaolei Mou:
### save Series 1 dimensional
for (i in 1:dim(CFD_data_Series1_dimensional)[3]) {
  slice_df <- as.data.frame(CFD_data_Series1_dimensional[,,i])
  file_name <- paste0(dimnames(CFD_data_Series1_dimensional)[[3]][i], ".csv")
  saving_location <- paste0(local_path, '/CFD data processed (for Xiaolei Mou)/CFD_Series_1_dimensional/')
  write.csv(slice_df, file = paste0(saving_location,file_name), row.names = F)
}
### save Series 2 dimensional
for (i in 1:dim(CFD_data_Series2_dimensional)[3]) {
  slice_df <- as.data.frame(CFD_data_Series2_dimensional[,,i])
  file_name <- paste0(dimnames(CFD_data_Series2_dimensional)[[3]][i], ".csv")
  saving_location <- paste0(local_path, '/CFD data processed (for Xiaolei Mou)/CFD_Series_2_dimensional/')
  write.csv(slice_df, file = paste0(saving_location,file_name), row.names = F)
}


## Series 3 dimensional
## (when considering species-specific WB freq, it should scale with mass)
p=1
plot((CFD_data_Series3_dimensional[,'time', p]-2)*100, CFD_data_Series3_dimensional[,'total_force', p]*2*1000, typ='l', lwd=2, ylim=c(0,3), col=colors_species[p], xlab=c('wingbeat cycle [%]'), ylab=c('force [mN]'), bty='n', main='series 3')
for (p in 2:dim(CFD_data_Series3_dimensional)[3]){points((CFD_data_Series3_dimensional[,'time', p]-2)*100, CFD_data_Series3_dimensional[,'total_force', p]*2*1000, typ='l', lwd=2, col=colors_species[p])}

## Series 3 dimensional normalized by m*g (should fall on top of each other)
p=1
plot((CFD_data_Series3_dimensional[,'time', p]-2)*100, (CFD_data_Series3_dimensional[,'total_force', p]*2)/((full_dataset$fresh_weight[p]*10^-6)*9.81), typ='l', lwd=2, ylim=c(0,5), col=colors_species[p], xlab=c('wingbeat cycle [%]'), ylab=c('weight-normalized force F / mg [-]'), bty='n', main='series 3')
for (p in 2:dim(CFD_data_Series3_dimensional)[3]){
  points((CFD_data_Series3_dimensional[,'time', p]-2)*100, (CFD_data_Series3_dimensional[,'total_force', p]*2)/((full_dataset$fresh_weight[p]*10^-6)*9.81), typ='l', lwd=2, col=colors_species[p])
  # abline(mean(CFD_data_Series3_dimensional[,'total_force', p]/((full_dataset$fresh_weight[p]*10^-6)*9.81)),0, col=colors_species[p], lty='longdash')
  # abline((full_dataset$force_dimens_Series3[p]*2)/((full_dataset$fresh_weight[p]*10^-6)*9.81), 0 , col=colors_species[p], lty='longdash')
}

###### save it for Xiaolei Mou:
### save Series 3 dimensional
for (i in 1:dim(CFD_data_Series3_dimensional)[3]) {
  slice_df <- as.data.frame(CFD_data_Series3_dimensional[,,i])
  file_name <- paste0(dimnames(CFD_data_Series3_dimensional)[[3]][i], ".csv")
  saving_location <- paste0(local_path, '/CFD data processed (for Xiaolei Mou)/CFD_Series_3_dimensional/')
  write.csv(slice_df, file = paste0(saving_location,file_name), row.names = F)
}


############ for supplementary, we want the 4 panels:
## mean WB freq force, mean WB freq normalized force
## mean sp-specific freq force, sp-specific freq normalized force
## because things are not following well on top of each other, we try

## Series 2 (mean WB freq)
## force normalized by mean force of each species (F*=F/F_mean)
p=1
plot((CFD_data_Series2_dimensional[,'time', p]-2)*100, (CFD_data_Series2_dimensional[,'total_force', p]*2)/mean(CFD_data_Series2_dimensional[,'total_force', p]), typ='l', ylim=c(0,5), lwd=2, col=colors_species[p], xlab=c('wingbeat cycle [%]'), ylab=c('normalised force [-] (F*=F/F_mean)'), bty='n', main='series 2 (mean WB freq)')
for (p in 2:dim(CFD_data_Series2_dimensional)[3]){points((CFD_data_Series2_dimensional[,'time', p]-2)*100, (CFD_data_Series2_dimensional[,'total_force', p]*2)/mean(CFD_data_Series2_dimensional[,'total_force', p]), typ='l', lwd=2, col=colors_species[p])}
## force normalized by WB-average force (F*=F/F_mean) (F*=F/F_mean) AND normalised by m*g
p=1
plot((CFD_data_Series2_dimensional[,'time', p]-2)*100, ((CFD_data_Series2_dimensional[,'total_force', p]*2)/mean(CFD_data_Series2_dimensional[,'total_force', p]))/((full_dataset$fresh_weight[p]*10^-6)*9.81), ylim=c(0,60000), typ='l', lwd=2, col=colors_species[p], xlab=c('wingbeat cycle [%]'), ylab=c('normalised force [-] (F*/mg)'), bty='n', main='series 2 (mean WB freq) divided by mg')
for (p in 2:dim(CFD_data_Series2_dimensional)[3]){
  points((CFD_data_Series2_dimensional[,'time', p]-2)*100, ((CFD_data_Series2_dimensional[,'total_force', p]*2)/mean(CFD_data_Series2_dimensional[,'total_force', p]))/((full_dataset$fresh_weight[p]*10^-6)*9.81), typ='l', lwd=2, col=colors_species[p])
} 


## Series 3 (species-specific WB freq)
## force normalized by WB-average force (F*=F/F_mean)
p=1
plot((CFD_data_Series3_dimensional[,'time', p]-2)*100, (CFD_data_Series3_dimensional[,'total_force', p]*2*1000)/mean(CFD_data_Series3_dimensional[,'total_force', p]), typ='l', lwd=2, col=colors_species[p], xlab=c('wingbeat cycle [%]'), ylab=c('normalised force [-] (F*=F/F_mean)'), bty='n', main='series 3 (species-specific WB freq)')
for (p in 2:dim(CFD_data_Series3_dimensional)[3]){points((CFD_data_Series3_dimensional[,'time', p]-2)*100, (CFD_data_Series3_dimensional[,'total_force', p]*2*1000)/mean(CFD_data_Series3_dimensional[,'total_force', p]), typ='l', lwd=2, col=colors_species[p])}
## force normalized by mean force of each species (F*=F/F_mean) AND normalised by m*g
p=1
plot((CFD_data_Series3_dimensional[,'time', p]-2)*100, ((CFD_data_Series3_dimensional[,'total_force', p]*2)/mean(CFD_data_Series3_dimensional[,'total_force', p]))/((full_dataset$fresh_weight[p]*10^-6)*9.81), typ='l', lwd=2, col=colors_species[p], ylim=c(0,60000), xlab=c('wingbeat cycle [%]'), ylab=c('normalised force [-] (F*/mg)'), bty='n', main='series 3 (species-specific WB freq) divided by mg')
for (p in 2:dim(CFD_data_Series3_dimensional)[3]){points((CFD_data_Series3_dimensional[,'time', p]-2)*100, ((CFD_data_Series3_dimensional[,'total_force', p]*2)/mean(CFD_data_Series3_dimensional[,'total_force', p]))/((full_dataset$fresh_weight[p]*10^-6)*9.81), typ='l', lwd=2, col=colors_species[p])}








######################## [4] average forces vs morphology ######################## 

############ ___ (i) force vs body mass ############ 

### aerodynamic force vs body mass
plot(full_dataset$fresh_weight, (full_dataset$force_dimens_Series3*1000)*2, pch=21, cex=3, ylim=c(0,2), bty='n', bg=colors_species, xlab='body mass [mg]', ylab='mean force [mN]')
## *2 because 2 wings

## fit a line
lmodel2((full_dataset$force_dimens_Series3*1000)*2 ~ full_dataset$fresh_weight)
abline(0.05437399 , 0.009214414       , col='black', lwd=2)   # observed
abline(0, 0.01, col='grey55', lwd=2, lty='dashed')   # weight support



### get the observed fitting line vs theoretical weight support line:
## for this we transformed the variable in weight [N] vs force [N]  
# plot((full_dataset$fresh_weight*10^-6)*9.81, full_dataset$force_dimens_Series3, pch=21, cex=3, ylim=c(0,1.5e-3), bty='n', bg=colors_species, xlab='weight [N]', ylab='force [mN]')
# Y = full_dataset$force_dimens_Series3
# X = (full_dataset$fresh_weight*10^-6)*9.81
# lmodel2(Y ~ X)
# abline(4.268134e-05 , 0.7259081 , col='black', lwd=2)   # observed
# abline(0 , 1, col='grey55', lwd=2, lty='dashed')


############ ___ (i) weight-normalized force vs body mass ############ 

### if force generate the required weight support, F/mg vs body mass should result in a horizontal relationship
## forces from Series 3 from me:
plot(full_dataset$fresh_weight, (full_dataset$force_dimens_Series3)*2/((full_dataset$fresh_weight*10^-6)*9.81), ylim=c(0,4), pch=21, cex=3, bty='n', bg=colors_species, xlab='body mass [mg]', ylab='mean weight-normalized force F/mg [-]')
## forces from Series 3 from Thomas:
plot(full_dataset$fresh_weight, (dimensional_data_series3$total_force_N)*2/((full_dataset$fresh_weight*10^-6)*9.81), ylim=c(0,3), pch=21, cex=3, bty='n', bg=colors_species, xlab='body mass [mg]', ylab='mean weight-normalized force F/mg [-]')
### (we multiply the force by 2 because 2 wings)


## fit a line
X = (full_dataset$force_dimens_Series3)*2/((full_dataset$fresh_weight*10^-6)*9.81)
lmodel2(X ~ full_dataset$fresh_weight )
abline(1.274011, -0.003680456, col='black', lwd=2)   # observed






### to show variation in WB freq on the aerodynamic results figure
plot(full_dataset$fresh_weight, full_dataset$WB_freq, pch=21, col='black', bg=colors_species, cex=3, ylim=c(0,250), bty='n', xlab='body mass [mg]', ylab='wingbeat frequency[Hz]')
abline(h = mean(full_dataset$WB_freq), col = "black", lty = 2)







############ ___ (iii) force vs wing parameters ############ 

### aerodynamic force vs wing parameters:
### mean frequency (Series 2)
plot(full_dataset$secondMoment ~ full_dataset$force_dimens_Series2, pch=21, cex=3, bty='n', bg=colors_species, ylab='S2 [cm^2]', xlab='force [N]')
plot(full_dataset$length_cm ~ full_dataset$force_dimens_Series2, pch=21, cex=3, bty='n', bg=colors_species, ylab='wing length [cm]', xlab='force [N]')
plot(full_dataset$average_chord ~ full_dataset$force_dimens_Series2, pch=21, cex=3, bty='n', bg=colors_species, ylab='average chord [cm]', xlab='force [N]')
plot(full_dataset$non_dimensional_2ndMoment ~ full_dataset$force_dimens_Series2, pch=21, cex=3, bty='n', bg=colors_species, ylab='dimensionless S2 [-]', xlab='force [N]')

### species-specific frequency (Series 3)
plot(full_dataset$force_dimens_Series3*2*1000, full_dataset$secondMoment, pch=21, cex=3, bty='n', bg=colors_species, ylab='S2 [cm^2]', xlab='mean vertical force [mN]')
plot(full_dataset$force_dimens_Series3*2*1000, full_dataset$length_cm, pch=21, cex=3, bty='n', bg=colors_species, ylab='wing length [cm]', xlab='mean vertical force [mN]')
plot(full_dataset$force_dimens_Series3*2*1000, full_dataset$average_chord, pch=21, cex=3, bty='n', bg=colors_species, ylab='average chord [cm]', xlab='mean vertical force [mN]')
plot(full_dataset$force_dimens_Series3*2*1000, full_dataset$non_dimensional_2ndMoment, pch=21, cex=3, bty='n', bg=colors_species, ylab='dimensionless S2 [-]', xlab='mean vertical force [mN]')


###### plots force vs wing parameters for paper figures:

# define function to plot detailed log scale while keeping real value:
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



## load the phylogeny so we can perform PGLS regression between wing parameters and aerodynamic force
## tree matching the small flight dataset:
load(paste0(substr(local_path, 1, as.numeric(gregexpr(pattern ='/04_CFD Thomas Engels',local_path)[[1]])), 'phylogeny from Daniel Wong et al. 2023/hoverfly_tree_(flight_dataset).RData'))
# plot(phy_filmed_only, cex=0.6, use.edge.length = F)

## reorder full_dataset to match phylogeny
full_dataset_tip_order = full_dataset
rownames(full_dataset_tip_order) = as.character(phy_filmed_only$tip.label)
for (i in 1:nrow(full_dataset_tip_order)){full_dataset_tip_order[i,] = full_dataset[which(full_dataset$genus_species==phy_filmed_only$tip.label[i]),]}

## also reorder color vector:
colors_species_tip_order = c('#F07E15', # Tropidia_scita
                             '#2B3A21', # Eristalis_tenax
                             '#A63D40', # Melanostoma_mellinum
                             '#7BCEEF', # Platycheirus_clypeatus
                             '#D1BA15', # Episyrphus_viridaureus
                             '#674ACE', # Sphaerophoria_scripta
                             '#287C34', # Syrphus_nigrilinearus
                             '#78C155') # Eupeodes_nielseni
colors_species_tip_order = colors_species_tip_order[length(colors_species_tip_order):1] # for some reason it's inversed




############  ________ S2 vs forces: ############ 

######## normal regression
plot(log10(full_dataset$secondMoment) ~ log10(full_dataset$force_dimens_Series3*2*1000), pch=21, cex=3, bty='n', bg=colors_species, ylab='log(S2)', xlab='log(mean vertical force)')
plot((full_dataset$force_dimens_Series3*2*1000), full_dataset$secondMoment, log='xy', pch=21, cex=3, bty='n', bg=colors_species, ylab='S2 [cm^4]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

lmodel2(log10(full_dataset$secondMoment) ~ log10(full_dataset$force_dimens_Series3*2*1000)) 
# check "Warton 2007. Bivariate line fitting methods for allometry" about SMA

abline(-0.8019922, 1.135573  , col='black', lwd=2) # observed (SMA row of the lmodel2() output)
abline(-0.8019922,(4/3), col='grey55', lwd=2, lty='dashed') # isometric

m <- lm(log10(full_dataset$secondMoment) ~ log10(full_dataset$force_dimens_Series3*2*1000)) ; summary(m)
confint(m)

######## phylogenetic regression

## PGLS using nlme (for some reason, the pgls using caper package doesn't work with hoverfly data)
## (we follow  https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(full_dataset_tip_order$force_dimens_Series3*2*1000); names(x_now) = phy_filmed_only$tip.label
y_now = log10(full_dataset_tip_order$secondMoment); names(y_now) = phy_filmed_only$tip.label
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
summary(lm(y_now ~ x_now))

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


#### plot including the different fitting lines
plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$secondMoment, log='xy', pch=21, cex=3, bty='n', bg=colors_species_tip_order, ylab='S2 [cm^2]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)
abline(-0.8395809, 1.1133493  , col='black', lwd=2) # PGLS
abline(-0.8441588, 1.063914  , col='black', lwd=2) # OLS
abline(-0.8395809,(4/3), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(-0.8395809,1, col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically





############  ________ wing length vs forces: ############ 

plot(log10(full_dataset$length_cm) ~ log10(full_dataset$force_dimens_Series3*2*1000), pch=21, cex=3, bty='n', bg=colors_species, ylab='log(wingspan)', xlab='log(mean vertical force)')
plot((full_dataset$force_dimens_Series3*2*1000), full_dataset$length_cm, log='xy', pch=21, cex=3, bty='n', bg=colors_species, ylab='wingspan [cm]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

lmodel2(log10(full_dataset$length_cm) ~ log10(full_dataset$force_dimens_Series3*2*1000))
# check "Warton 2007. Bivariate line fitting methods for allometry" about SMA

abline(0.06874105, 0.2899559  , col='black', lwd=2) # observed (SMA row of the lmodel2() output)
abline(0.06874105,(1/3), col='grey55', lwd=2, lty='dashed') # isometric


m <- lm(log10(full_dataset$length_cm) ~ log10(full_dataset$force_dimens_Series3*2*1000)) ; summary(m)
confint(m)




######## phylogenetic regression

## PGLS using nlme (for some reason, the pgls using caper package doesn't work with hoverfly data)
## (we follow  https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(full_dataset_tip_order$force_dimens_Series3*2*1000); names(x_now) = phy_filmed_only$tip.label
y_now = log10(full_dataset_tip_order$length_cm); names(y_now) = phy_filmed_only$tip.label
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
summary(lm(y_now ~ x_now))

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


#### plot including the different fitting lines
plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$length_cm, log='xy', pch=21, cex=3, bty='n', bg=colors_species_tip_order, ylab='wingspan [cm]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)
abline(0.04239217, 0.27634836  , col='black', lwd=2) # PGLS
abline(0.05852589, 0.2725960  , col='black', lwd=2) # OLS
abline(0.04239217,(1/3), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(0.04239217, (2/9), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically




############  ________ wing chord vs forces: ############ 

plot(log10(full_dataset$average_chord) ~ log10(full_dataset$force_dimens_Series3*2*1000), pch=21, cex=3, bty='n', bg=colors_species, ylab='log(wing chord)', xlab='log(mean vertical force)')
plot((full_dataset$force_dimens_Series3*2*1000), full_dataset$average_chord, log='xy', pch=21, cex=3, bty='n', bg=colors_species, ylab='wing chord [cm]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

lmodel2(log10(full_dataset$average_chord) ~ log10(full_dataset$force_dimens_Series3*2*1000))
# check "Warton 2007. Bivariate line fitting methods for allometry" about SMA

abline(-0.5095045, 0.2976487  , col='black', lwd=2) # observed (SMA row of the lmodel2() output)
abline(-0.5095045,(1/3), col='grey55', lwd=2, lty='dashed') # isometric


m <- lm(log10(full_dataset$average_chord) ~ log10(full_dataset$force_dimens_Series3*2*1000)) ; summary(m)
confint(m)



######## phylogenetic regression

## PGLS using nlme (for some reason, the pgls using caper package doesn't work with hoverfly data)
## (we follow  https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(full_dataset_tip_order$force_dimens_Series3*2*1000); names(x_now) = phy_filmed_only$tip.label
y_now = log10(full_dataset_tip_order$average_chord); names(y_now) = phy_filmed_only$tip.label
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
summary(lm(y_now ~ x_now))

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


#### plot including the different fitting lines
plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$average_chord, log='xy', pch=21, cex=3, bty='n', bg=colors_species_tip_order, ylab='average chord [cm]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)
abline(-0.52522, 0.27095  , col='black', lwd=2) # PGLS
abline(-0.5252156, 0.2709487  , col='black', lwd=2) # OLS
abline(-0.52522,(1/3), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(-0.52522, (0), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically




############  ________ non-dim S2 vs forces: ############ 

plot(log10(full_dataset$non_dimensional_2ndMoment) ~ log10(full_dataset$force_dimens_Series3*2*1000), pch=21, cex=3, bty='n', bg=colors_species, ylab='log(dimensionless S2)', xlab='log(mean vertical force)')
plot((full_dataset$force_dimens_Series3*2*1000), full_dataset$non_dimensional_2ndMoment, log='xy', pch=21, cex=3, bty='n', bg=colors_species, ylab='dimensionless S2 [-]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

lmodel2(log10(full_dataset$non_dimensional_2ndMoment) ~ log10(full_dataset$force_dimens_Series3*2*1000)) ## check the SMA line in the output
# check "Warton 2007. Bivariate line fitting methods for allometry" about SMA

abline(-0.2562559, -0.02168560  , col='black', lwd=2) # observed (SMA row of the lmodel2() output)
abline(mean(log10(full_dataset$non_dimensional_2ndMoment)),(0), col='grey55', lwd=2, lty='dashed') # isometric


m <- lm(log10(full_dataset$non_dimensional_2ndMoment) ~ log10(full_dataset$force_dimens_Series3*2*1000)) ; summary(m)
confint(m)


full_dataset$non_dimensional_2ndMoment[which(full_dataset$genus_species=='Platycheirus_clypeatus')]





######## phylogenetic regression

## PGLS using nlme (for some reason, the pgls using caper package doesn't work with hoverfly data)
## (we follow  https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(full_dataset_tip_order$force_dimens_Series3*2*1000); names(x_now) = phy_filmed_only$tip.label
y_now = log10(full_dataset_tip_order$non_dimensional_2ndMoment); names(y_now) = phy_filmed_only$tip.label
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
summary(lm(y_now ~ x_now))

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


#### plot including the different fitting lines
plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$non_dimensional_2ndMoment, log='xy', pch=21, cex=3, bty='n', bg=colors_species_tip_order, ylab='non-dim S2 [-]', xlab='mean vertical force [mN]', axes=F, ylim=c(0.25,0.4))
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2)
# abline(-0.5034996, -0.0230802  , col='black', lwd=2) # PGLS
abline(-0.5047312, -0.03014837  , col='black', lwd=2) # OLS
abline(-0.5047312,(0), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(-0.5047312, (-0.33), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically


##### adjust Y-axis range

# Convert variables to log scale
x_log <- log10(full_dataset_tip_order$force_dimens_Series3 * 2 * 1000)
y_log <- log10(full_dataset_tip_order$non_dimensional_2ndMoment)

lmodel2(y_log ~ x_log)

# Plot with log-log scale
plot((full_dataset_tip_order$force_dimens_Series3 * 2 * 1000), 
     full_dataset_tip_order$non_dimensional_2ndMoment, 
     log='xy', pch=21, cex=3, bty='n', bg=colors_species_tip_order, 
     ylab='non-dim S2 [-]', xlab='mean vertical force [mN]', axes=F, ylim=c(0.5, 0.6))
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

## recompute the PGLS with the log 10 values:
pglsModel <- gls(y_log ~ x_log, correlation = corBrownian(phy = phy_filmed_only), method = "ML")
summary(pglsModel)
abline(-0.25174978, -0.01154012 , col='black', lwd=2) # PGLS
abline(-0.2523656, -0.01507419  , col='black', lwd=2) # OLS
abline(-0.25174978,(0), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(-0.25174978, (-0.33), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically




full_dataset$non_dimensional_2ndMoment[which(full_dataset$genus_species=='Eristalis_tenax')]









############ ___ (iv) force vs kinematic parameters ############ 


############  ________ WB frequency vs forces: ############ 


plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$WB_freq, log='xy', pch=21, cex=3, bty='n', bg=colors_species, ylab='wingbeat frequency [Hz]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)


m <- lm(log10(full_dataset$WB_freq) ~ log10(full_dataset$force_dimens_Series3*2*1000)) ; summary(m)
confint(m)

######## phylogenetic regression

## PGLS using nlme (for some reason, the pgls using caper package doesn't work with hoverfly data)
## (we follow  https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(full_dataset_tip_order$force_dimens_Series3*2*1000); names(x_now) = phy_filmed_only$tip.label
y_now = log10(full_dataset_tip_order$WB_freq); names(y_now) = phy_filmed_only$tip.label
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
summary(lm(y_now ~ x_now))

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


#### plot including the different fitting lines
plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$WB_freq, log='xy', pch=21, cex=3, bty='n', bg=colors_species_tip_order, ylab='wingbeat frequency [Hz]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)
abline(2.252515, -0.042664  , col='black', lwd=2) # PGLS
abline(2.21887, -0.03715  , col='black', lwd=2) # OLS
abline(2.252515,(0), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(2.252515, (-0.17), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically

### adjust Y-axis range
plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$WB_freq, log='xy', pch=21, cex=3, bty='n', bg=colors_species_tip_order, ylab='wingbeat frequency [Hz]', xlab='mean vertical force [mN]',
     ylim=c(30, 400), axes=F) # we exclude zero because it would not work in log scale
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)
abline(2.252515, -0.042664  , col='black', lwd=2) # PGLS
abline(2.21887, -0.03715  , col='black', lwd=2) # OLS
abline(2.252515,(0), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(2.252515, (-0.17), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically






############  ________ stroke amplitude vs forces: ############ 


plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$stroke_amplitude, log='xy', pch=21, cex=3, bty='n', bg=colors_species, ylab='stroke amplitude [deg]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)


######## phylogenetic regression

## PGLS using nlme (for some reason, the pgls using caper package doesn't work with hoverfly data)
## (we follow  https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(full_dataset_tip_order$force_dimens_Series3*2*1000); names(x_now) = phy_filmed_only$tip.label
y_now = log10(full_dataset_tip_order$stroke_amplitude); names(y_now) = phy_filmed_only$tip.label
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
summary(lm(y_now ~ x_now))

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


#### plot including the different fitting lines
plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$stroke_amplitude, log='xy', pch=21, cex=3, bty='n', bg=colors_species_tip_order, ylab='stroke amplitude [deg]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)
abline(2.252515, -0.042664  , col='black', lwd=2) # PGLS
abline(2.252515,(0), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(2.252515, (-0.17), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically






############  ________ angular speed (fxA) vs forces: ############ 

full_dataset_tip_order$angular_speed_fA <- full_dataset_tip_order$WB_freq * full_dataset_tip_order$stroke_amplitude

plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$angular_speed_fA, log='xy', pch=21, cex=3, bty='n', bg=colors_species, ylab='angular speed [deg/sec]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)


######## phylogenetic regression

## PGLS using nlme (for some reason, the pgls using caper package doesn't work with hoverfly data)
## (we follow  https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(full_dataset_tip_order$force_dimens_Series3*2*1000); names(x_now) = phy_filmed_only$tip.label
y_now = log10(full_dataset_tip_order$angular_speed_fA); names(y_now) = phy_filmed_only$tip.label
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
summary(lm(y_now ~ x_now))

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


#### plot including the different fitting lines
plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$angular_speed_fA, log='xy', pch=21, cex=3, bty='n', bg=colors_species_tip_order, ylab='angular speed [deg/sec]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)
abline(4.223616, -0.063396  , col='black', lwd=2) # PGLS
abline(4.223616,(0), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(4.223616, (-0.17), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically







############  ________ AoA vs forces: ############ 


plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$AoA_peak_ang_velo_bs, log='xy', pch=21, cex=3, bty='n', bg=colors_species, ylab='angle-of-attack [deg]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)


######## phylogenetic regression

## PGLS using nlme (for some reason, the pgls using caper package doesn't work with hoverfly data)
## (we follow  https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(full_dataset_tip_order$force_dimens_Series3*2*1000); names(x_now) = phy_filmed_only$tip.label
y_now = log10(full_dataset_tip_order$AoA_peak_ang_velo_bs); names(y_now) = phy_filmed_only$tip.label
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
summary(lm(y_now ~ x_now))

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


#### plot including the different fitting lines
plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$AoA_peak_ang_velo_bs, log='xy', pch=21, cex=3, bty='n', bg=colors_species_tip_order, ylab='angle-of-attack [deg]', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)
abline(1.53186, -0.07424  , col='black', lwd=2) # PGLS
abline(1.53186,(0), col='#ba4c36', lwd=2, lty='dashed') # isometric
abline(1.53186, (-0.17), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically





############  ________ (frequency)^2 x S2 vs force ############ 

full_dataset_tip_order$freq_times_S2 <- (full_dataset_tip_order$WB_freq)^2 * full_dataset_tip_order$secondMoment

plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$freq_times_S2, log='xy', pch=21, cex=3, bty='n', bg=colors_species_tip_order, ylab='WB frequency * S2', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)

Y = log10(full_dataset_tip_order$freq_times_S2)
X = log10(full_dataset_tip_order$force_dimens_Series3*2*1000)
m <- lm(Y ~ X) ; summary(m)
confint(m)


######## phylogenetic regression

## PGLS using nlme (for some reason, the pgls using caper package doesn't work with hoverfly data)
## (we follow  https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/)
x_now = log10(full_dataset_tip_order$force_dimens_Series3*2*1000); names(x_now) = phy_filmed_only$tip.label
y_now = log10(full_dataset_tip_order$freq_times_S2); names(y_now) = phy_filmed_only$tip.label
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
summary(lm(y_now ~ x_now))


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


#### plot including the different fitting lines
plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$freq_times_S2, log='xy', pch=21, cex=3, bty='n', bg=colors_species_tip_order, ylab='(WB frequency)^2 * S2', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)
abline(3.594381, 0.991612  , col='black', lwd=2) # PGLS
abline(3.593575, 0.989607  , col='black', lwd=2) # OLS
abline(3.594381,(1.33), col='#ba4c36', lwd=2, lty='dashed') # isometric
# abline(1.37740, (-0.17), col='#A0A0A0', lwd=2, lty='dotted')# required for weight support if everything else scales isometrically





############  ________ (angular speed)^2 x S2 vs force ############ 

full_dataset_tip_order$omega_times_S2 <- (full_dataset_tip_order$angular_speed_fA)^2 * full_dataset_tip_order$secondMoment

plot((full_dataset_tip_order$force_dimens_Series3*2*1000), full_dataset_tip_order$freq_times_S2, log='xy', pch=21, cex=3, bty='n', bg=colors_species_tip_order, ylab='(angular speed)^2 * S2', xlab='mean vertical force [mN]', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)


## here w=f*A, with f the species-specific frequency, and A the mean amplitude (as used in CFD).
## Because A is thus constant, the results should look exactly the same as when using (frequency)^2 x S2.

#### give data used for FIgure 7 to Florian
data <- cbind(as.character(full_dataset_tip_order$genus_species),
              full_dataset_tip_order$force_dimens_Series2*2*1000,
              full_dataset_tip_order$force_dimens_Series3*2*1000,
              full_dataset_tip_order$WB_freq,
              full_dataset_tip_order$stroke_amplitude,
              full_dataset_tip_order$length_cm,
              full_dataset_tip_order$average_chord,
              full_dataset_tip_order$secondMoment,
              full_dataset_tip_order$non_dimensional_2ndMoment)
colnames(data) = c('genus_species', 'Fv_meanFreq_mN', 'Fv_variableFreq_mN', 'frequency_Hz', 'amplitude_deg', 'wingspan_cm', 'mean_chord_cm', 'S2', 'dimensionaless_S2')
# write.csv(data, paste0(local_path, '/data_fig_7.csv'))


######################## [5] effect of including WB freq variation on weight support estimation ######################## 


#### show difference between Series 3 and Series 2
plot(full_dataset$fresh_weight, (full_dataset$force_dimens_Series2*2)/((full_dataset$fresh_weight*10^-6)*9.81), ylim=c(0,3), pch=22, cex=3, bty='n', bg=colors_species, xlab='body mass [mg]', ylab='mean weight-normalized force F/mg [-]')
points(full_dataset$fresh_weight, (full_dataset$force_dimens_Series3*2)/((full_dataset$fresh_weight*10^-6)*9.81), pch=21, cex=3, bg=colors_species)
## fit the line for series 3 (circles):
abline(lm(((full_dataset$force_dimens_Series3*2)/((full_dataset$fresh_weight*10^-6)*9.81)) ~ full_dataset$fresh_weight ), col='black', lwd=2)
## fit the line for series 2 (squares):
abline(lm(((full_dataset$force_dimens_Series2*2)/((full_dataset$fresh_weight*10^-6)*9.81)) ~ full_dataset$fresh_weight ), col='black', lwd=2, lty='longdash')
## show the weight support line
abline(1, 0, col='grey55', lwd=2, lty='dashed')


## compute the deviation from weight support (as abs(F/mg - 1)*100) for both Series 2 and Series 3 
data_dev_from_weight_supp <- as.data.frame(matrix(NA, 16, 2))
colnames(data_dev_from_weight_supp) = c('Series', 'values')
data_dev_from_weight_supp$Series = c(rep('Series_2',8), rep('Series_3',8))

for (i in which(data_dev_from_weight_supp$Series=='Series_2')) {
  dev_from_Fmg_now = ((full_dataset$force_dimens_Series2*2)/((full_dataset$fresh_weight*10^-6)*9.81))[i] - 1
  data_dev_from_weight_supp$values[i] = abs(dev_from_Fmg_now)*100
}
for (i in which(data_dev_from_weight_supp$Series=='Series_3')) {
  dev_from_Fmg_now = ((full_dataset$force_dimens_Series3*2)/((full_dataset$fresh_weight*10^-6)*9.81))[i-8] - 1
  data_dev_from_weight_supp$values[i] = abs(dev_from_Fmg_now)*100
}

ggplot(data_dev_from_weight_supp, aes(x = Series, y = values, shape = Series)) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Series)) +
  geom_jitter(position = position_jitter(0.1), cex = 4) +
  labs(y = 'deviation from weight support: abs(F/mg -1)', x=' ') +  ylim(0,150) +
  scale_shape_manual(values =  c(22,21) ) + # specify pch
  theme_classic() + theme(legend.position = "none")

### --> does not show what we were hoping for!


#### compare variation in F/mg
weight_supp_for_boxplot <- data_dev_from_weight_supp
for (i in 1:8) { weight_supp_for_boxplot$values[i] = (full_dataset$force_dimens_Series2[i]*2)/((full_dataset$fresh_weight[i]*10^-6)*9.81) }
for (i in 9:16) { weight_supp_for_boxplot$values[i] = (full_dataset$force_dimens_Series3[i-8]*2)/((full_dataset$fresh_weight[i-8]*10^-6)*9.81) }

ggplot(weight_supp_for_boxplot, aes(x = Series, y = values, shape = Series)) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Series)) +
  geom_jitter(position = position_jitter(0.1), cex = 4) +
  labs(y = 'mean weight-normalized force F/mg [-]', x=' ') +  ylim(0,3) +
  scale_shape_manual(values =  c(22,21) ) + # specify pch
  theme_classic() + theme(legend.position = "none")





#### Differences are better seen by plotting force vs mass of series 2 and series 3

plot(full_dataset$fresh_weight, (full_dataset$force_dimens_Series2*2*1000), ylim=c(0,2), pch=22, cex=3, bty='n', bg=colors_species, xlab='body mass [mg]', ylab='mean vertical force [mN]')
points(full_dataset$fresh_weight, (full_dataset$force_dimens_Series3*2*1000), pch=21, cex=3, bg=colors_species)
abline(0, 0.01, col='grey55', lwd=2, lty='dashed')# weight support
abline(lm((full_dataset$force_dimens_Series3*2*1000) ~ full_dataset$fresh_weight), col='black', lwd=2)# series 3 fitting
abline(lm((full_dataset$force_dimens_Series2*2*1000) ~ full_dataset$fresh_weight), col='black', lwd=2, lty='longdash')# series 2 fitting




######################## [8] sum up the relative contribution of R c and S2 ######################## 

## we analytically estimated the contribution of R, c and S2* to the S2 and hence to force production.
## we show the observed contribution and the theoretical contribution under isometry

contrib_wing_param <- as.data.frame(matrix(NA, 4*2, 3))
colnames(contrib_wing_param) = c('parameters', 'type', 'values')
contrib_wing_param$parameters = c('wing_length', 'wing_chord', 'nd_S2', 'R_c_nd_S2')
contrib_wing_param$parameters <- factor(contrib_wing_param$parameters, levels = c('wing_length', 'wing_chord', 'nd_S2', 'R_c_nd_S2'))
contrib_wing_param$type <- c(rep('observed',4), rep('isometry',4))
contrib_wing_param$values[which(contrib_wing_param$type=='observed')] = c(73, 16, 11, 100)
contrib_wing_param$values[which(contrib_wing_param$type=='isometry')] = c(75, 25, 0, 100)

ggplot(contrib_wing_param, aes(x = parameters, y=values, fill = factor(type))) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  labs(title = '', x = ' ', y = 'contribution to aerodynamic force production [%]') +
  scale_fill_manual(values = c("observed" = '#0B0C0C', "isometry" = '#8C8C8C')) +
  theme_classic()





######################## [8] sum up the relative contribution of R c and S2 AND WB FREQ ########################


## we assessed the relative contribution of morphological and kinematic traits in producing weight support
## across hoverflies using the scaling exponents
## (see manuscript and excel file "Diptera_WeightSupport_AllometricScaling_FLORIAN_2.xlsx" for the methods).

## make a visualization of those relative contributions:

## contribution of wingbeat frequency and S2
relativ_contrib <- as.data.frame(matrix(NA,2,2))
colnames(relativ_contrib) = c('parameter', 'contribution')
relativ_contrib$parameter <- as.factor(c('S2', 'frequency'))
relativ_contrib$contribution <- as.numeric(c(66, 37))
### THOSE VALUES ARE COMPUTED IN "Diptera_WeightSupport_AllometricScaling_FLORIAN_2.xlsx"

## reorder levels:
relativ_contrib$parameter <- factor(relativ_contrib$parameter,  levels = c('frequency', 'S2'))

### pie chart of relative contributions
ggplot(relativ_contrib, aes(x = "", y = contribution, fill = parameter)) +
  geom_bar(width = 1, stat = "identity") +  # Make a bar chart with one "slice"
  coord_polar("y", start = 0) +             # Transform to a pie chart
  theme_void() +                            # Remove background and axes
  labs(title = "Relative Contribution of Parameters") +
  scale_fill_manual(values = c("grey50", "grey72")) +  # Shades of grey
  theme(legend.title = element_blank())     # Remove legend title




## contribution of wingspan, chord, non-dim S2, and frequency
relativ_contrib <- as.data.frame(matrix(NA,4,2))
colnames(relativ_contrib) = c('parameter', 'contribution')
relativ_contrib$parameter <- as.factor(c('wingspan', 'frequency', 'chord', 'non-dim S2'))
relativ_contrib$contribution <- as.numeric(c(44, 37, 18, 3))
### THOSE VALUES ARE COMPUTED IN "Diptera_WeightSupport_AllometricScaling_FLORIAN_2.xlsx"

## reorder levels:
relativ_contrib$parameter <- factor(relativ_contrib$parameter,  levels = c('wingspan', 'frequency', 'chord', 'non-dim S2'))


### pie chart of relative contributions
ggplot(relativ_contrib, aes(x = "", y = contribution, fill = parameter)) +
  geom_bar(width = 1, stat = "identity") +  # Make a bar chart with one "slice"
  coord_polar("y", start = 0) +             # Transform to a pie chart
  theme_void() +                            # Remove background and axes
  labs(title = "Relative Contribution of Parameters") +
  scale_fill_manual(values = c('grey35', 'grey50', 'grey72', 'grey85')) +  # Shades of grey
  theme(legend.title = element_blank())     # Remove legend title











