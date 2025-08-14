




rm(list=ls())




library(R.matlab)
library(ggplot2)
library(rgl)
library(pracma)
library(dplyr)
library(tagcloud)
library(RColorBrewer)
library(PerformanceAnalytics)
library(factoextra)
library(scico)


local_path = dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(local_path, '/AAhov_formatting_wingbeats.R')) # runs 'AAhov_formatting_wingbeats.R'



############################# [1] import the fourier-series-fitted wingbeats #############################


## load the file generated using the MATLAB code
## 'AA_batch_fit_fourier_wb.m' on the wingbeat extracted during the previous step:

## **the order of magnitude of the fourier series is for now 4, we could debate on whether the 2, 3, 4 5 or 6 are better fitting the data...
## 4th order of magnitude for all three angles (444):
fitted_wb_from_MATLAB <- readMat(paste0(local_path,'/02_fit_fourier_series_hoverflies/data/wb_one_by_one/tracked_and_fourier_fits_all_wingbeats_onebyone_444.mat'), h=T, sep=';')


## the list 'fitted_wb_from_MATLAB' contains the following elements:
##  - str/dev/rot.tracked are the tracked (raw) angle values. each column is a wingbeat.
##    the length of each tracked wingbeat is different.
##  - time.tracked contains the corresponding time columns associated with the tracked values.
##
##  - str/dev/rot.fit are the fitted ('smoothed') angle values. each column is a wingbeat.
##    they all have the same length of 101 time step.
##  - time.fit contains the corresponding time columns associated with the fitted values (all of length 101).
##
##  - str/dev/rot.fitresults is a list of list. each list in the list contains the model associated with each curve fitting. 
##    apparently it cannot be open in R. see the MATLAB object directly in MATLAB 
##  - str/dev/rot.gofs is a list of list. each list in the list contains the goodness of fit (rmse, r-squared etc)
##    associated with the model. it can be accessed using for example 'fitted_wb_from_MATLAB$dev.gofs[[1]]' for the first wb
##  - str/dev/rot.fcoeffs contains the 9 coefficients (a0,a1,a2,a3,a4 and b1,b2,b3,b4) of the model
##    it is a vector of 9 values for each wingbeats
##
##  - str/dev/rot.dot.fit contains the derivative of the fitted values, i.e. the stroke/dev/rot rate (or velocity)
##  - str/dev/rot.dotdot.fit contains the 2nd derivative of the fitted values, i.e. the stroke/dev/rot acceleration

### create a new array compiling the fitted data:
WB_array_fitted = array(NA, dim = c(nrow(fitted_wb_from_MATLAB$time.fit), 14, dim(WB_array)[3]))
dimnames(WB_array_fitted)[[2]] <-  c('wingbeat_cycle', 'stroke_fitted', 'deviation_fitted', 'rotation_fitted',
                                     'stroke_rate_fitted', 'deviation_rate_fitted', 'rotation_rate_fitted',
                                     'stroke_accel_fitted', 'deviation_accel_fitted', 'rotation_accel_fitted',
                                     'angle_of_attack','angle_of_attack_2','abs_wing_speed', 'body_pitch')
dimnames(WB_array_fitted)[[3]] <- dimnames(WB_array)[[3]]

## fill the array:
for (i in 1:dim(WB_array)[3]){
  WB_array_fitted[ , 'wingbeat_cycle', i] = fitted_wb_from_MATLAB$time.fit*100
  WB_array_fitted[ , 'stroke_fitted', i] = fitted_wb_from_MATLAB$str.fit[,i]
  WB_array_fitted[ , 'deviation_fitted', i] = fitted_wb_from_MATLAB$dev.fit[,i]
  WB_array_fitted[ , 'rotation_fitted', i] = fitted_wb_from_MATLAB$rot.fit[,i]
  WB_array_fitted[ , 'stroke_rate_fitted', i] = fitted_wb_from_MATLAB$str.dot.fit[,i]
  WB_array_fitted[ , 'deviation_rate_fitted', i] = fitted_wb_from_MATLAB$dev.dot.fit[,i]
  WB_array_fitted[ , 'rotation_rate_fitted', i] = fitted_wb_from_MATLAB$rot.dot.fit[,i]
  WB_array_fitted[ , 'stroke_accel_fitted', i] = fitted_wb_from_MATLAB$str.dotdot.fit[,i]
  WB_array_fitted[ , 'deviation_accel_fitted', i] = fitted_wb_from_MATLAB$dev.dotdot.fit[,i]
  WB_array_fitted[ , 'rotation_accel_fitted', i] = fitted_wb_from_MATLAB$rot.dotdot.fit[,i]
  WB_array_fitted[10:89 , 'body_pitch', i] = as.numeric(WB_array[ , 'body_pitch' , i])
  ## !! the body pitch value are not to be read according to 'wingbeat_cycle' (there is no match with the WB cycle)
}
# View(WB_array_fitted[,,1])


IDs_file_seq_lvl <- read.csv(paste0(local_path,'/ID_file_hoverflies_seq_level.csv'),h=T, sep=',')
IDs_file_seq_lvl$genus_species = as.factor(IDs_file_seq_lvl$genus_species)





### ---> the angle values  obtained from MATLAB code are in degrees per cycle.

### ---> we want the angle rate in deg/sec:
### we have to multiply angle rates by the wingbeat frequency, and angle accelerations by (wingbeat frequency)^2
### get frequency:
WB_freq <- IDs_file_seq_lvl$full_seq_WB_freq # this freq is more accurate because computed on the full sequence
# WB_freq <- rep(NA, dim(WB_array)[3])
# for (i in 1:length(WB_freq)){
#   time_step_now = as.numeric(WB_array_norma[2,col_time,i]) - as.numeric(WB_array_norma[1,col_time,i])
#   frame_rate_now = 1/time_step_now
#   # wingbeat frequency is computed as 1/time needed to complete a full wingbeat
#   WB_freq[i] = 1 / ((norma_wb_length[i]-2*X) * time_step_now)
# }

### convert the angle rates by multiplying by the WB_freq:
for (i in 1:dim(WB_array)[3]) {
  WB_array_fitted[,'stroke_rate_fitted',i] = WB_array_fitted[,'stroke_rate_fitted',i] * WB_freq[i]
  WB_array_fitted[,'deviation_rate_fitted',i] = WB_array_fitted[,'deviation_rate_fitted',i] * WB_freq[i]
  WB_array_fitted[,'rotation_rate_fitted',i] = WB_array_fitted[,'rotation_rate_fitted',i] * WB_freq[i]
  # correct the angle acceleration by multiplying by (wingbeat frequency)^2:
  WB_array_fitted[,'stroke_accel_fitted',i] = WB_array_fitted[,'stroke_accel_fitted',i] * (WB_freq[i])^2
  WB_array_fitted[,'deviation_accel_fitted',i] = WB_array_fitted[,'deviation_accel_fitted',i] * (WB_freq[i])^2
  WB_array_fitted[,'rotation_accel_fitted',i] = WB_array_fitted[,'rotation_accel_fitted',i] * (WB_freq[i])^2
}





############################# [2] wingbeat dynamic #############################

##### set wingbeat frequency in a color gradient:
variable_col_gradient = WB_freq

# # palette_now <- brewer.pal(9,'YlOrRd')[4:9] 
# palette_now <- brewer.pal(9,'Blues')[4:9]  #choose gradient among: display.brewer.all()
# # palette_now <- scico(9, palette = 'lapaz')[1:5] #choose gradient among: scico_palette_show()
# gradient_now <- smoothPalette(variable_col_gradient, palette_now)
# t <- colorRampPalette(c(palette_now[1], palette_now[length(palette_now)]))
# seq_now = seq(min(variable_col_gradient), max(variable_col_gradient))

##### set average weight in a color gradient:
variable_col_gradient = IDs_file_seq_lvl$mean_weight_mg
palette_now <- scico(25, palette = 'bamako')[c(1,3,6,8,10,12,13,15,17,19)]
# palette_now = palette_now[length(palette_now):1] # reverse the gradient
palette_now <- brewer.pal(9,'BuPu')[4:9]  #choose gradient among: display.brewer.all()
gradient_now <- smoothPalette(variable_col_gradient, palette_now)
t <- colorRampPalette(c(palette_now[1], palette_now[length(palette_now)]))
seq_now = seq(min(variable_col_gradient), max(variable_col_gradient))



################ ____ stroke angle ############### 

p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_fitted', p], typ='l', lwd=2, ylim=c(-100,100), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('stroke angle [deg]'), bty='n', main='stroke angle')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}

layout(matrix(1:2,ncol=2), width = c(2,1), height = c(1,1))
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_fitted', p], typ='l', lwd=2, ylim=c(-100,100), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('stroke angle [deg]'), bty='n', main='stroke angle')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
legend_image <- as.raster(matrix(t(length(seq_now)), ncol=1))
plot(c(0,2), c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'body mass [mg]', cex.main = 1)
text(x=1.4, y = seq(0,1,l=5), labels = seq(0, 100, l=5))
rasterImage(legend_image[length(seq_now):1], 0.3, 0, 0.7, 1) # inverse color direction here with "[length(seq_now):1]"
layout(matrix(1,1))

## check wingbeats one-by-one:
# for (p in dim(WB_array_fitted)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_fitted', p], typ='l', lwd=2, ylim=c(-100,100), col='dodgerblue4', xlab=c('wingbeat cycle [%]'), ylab=c('stroke angle [deg]'), bty='n', main=paste0(specimen_now,' - ',sequence_now))
#   points(as.numeric(WB_array_norma[(X+1):(norma_wb_length[p]-X),col_wingbeat_cycle, p]), as.numeric(WB_array_norma[(X+1):(norma_wb_length[p]-X),col_stroke_mean, p]), pch=16, cex=1, col='darkred')
# }






################ ____ stroke rate ########################

## stroke rate (or angular velocity) is the rate of change of the stroke angle.
## the angular velocity is usually expressed in radian/sec or in deg/sec.
## It has already been computed in MATLAB as the derivative or the fitted stroke angle

p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_rate_fitted', p], typ='l', lwd=2, ylim=c(-1e+05, 1e+05), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('stroke rate [deg/sec]'), bty='n', main='stroke rate')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_rate_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}

### ** to see the stroke rate in gradient, multiply it by *((pi/180)*WB_freq[i])

layout(matrix(1:2,ncol=2), width = c(2,1), height = c(1,1))
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_rate_fitted', p], typ='l', lwd=2, ylim=c(-1e+05, 1e+05), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('stroke rate [deg/sec]'), bty='n', main='stroke rate')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_rate_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
legend_image <- as.raster(matrix(t(length(seq_now)), ncol=1))
plot(c(0,2), c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'WB frequency [Hz]', cex.main = 1)
text(x=1.4, y = seq(0,1,l=5), labels = seq(100, 300, l=5))
rasterImage(legend_image[length(seq_now):1], 0.3, 0, 0.7, 1) # inverse color direction here with "[length(seq_now):1]"
layout(matrix(1,1))

# # check wingbeats one-by-one:
# for (p in dim(WB_array_fitted)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_rate_fitted', p], typ='l', lwd=2, ylim=c(-140000,140000), col='darkorchid4', xlab=c('wingbeat cycle [%]'), ylab=c('stroke rate [deg/sec]'), bty='n', main=paste0(specimen_now,' - ',sequence_now))
# }


##### compute and plot the average stroke rate at each frame (mean of all wingbeats)

average_stroke_rate = as.data.frame(matrix(NA,nrow(WB_array_fitted)[1], 2))
colnames(average_stroke_rate) <-  c('wingbeat_cycle', 'average_stroke_rate')
average_stroke_rate$wingbeat_cycle = as.numeric(WB_array_fitted[,'wingbeat_cycle', 1])
for (i in 1:nrow(average_stroke_rate)){
  average_stroke_rate$average_stroke_rate[i] = mean(WB_array_fitted[i,'stroke_rate_fitted',])
}

## plot stroke rate and its average
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_rate_fitted', p], typ='l', lwd=2, ylim=c(-1e+05, 1e+05), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('stroke rate [deg/sec]'), bty='n', main='stroke rate')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_rate_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
points(average_stroke_rate$wingbeat_cycle, average_stroke_rate$average_stroke_rate, typ='l', lwd=4, col='black')




################ ____ stroke acceleration ########################

## stroke acceleration has been computed in MATLAB as the 2nd derivative or the fitted stroke angle

p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_accel_fitted', p], typ='l', lwd=2, ylim=c(-2e+08, 2e+08), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('stroke acceleration [deg/sec^-2]'), bty='n', main='stroke acceleration')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_accel_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}

layout(matrix(1:2,ncol=2), width = c(2,1), height = c(1,1))
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_accel_fitted', p], typ='l', lwd=2, ylim=c(-2e+08, 2e+08), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('stroke acceleration [deg/sec^-2]'), bty='n', main='stroke acceleration')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_accel_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
legend_image <- as.raster(matrix(t(length(seq_now)), ncol=1))
plot(c(0,2), c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'WB frequency [Hz]', cex.main = 1)
text(x=1.4, y = seq(0,1,l=5), labels = seq(100, 300, l=5))
rasterImage(legend_image[length(seq_now):1], 0.3, 0, 0.7, 1) # inverse color direction here with "[length(seq_now):1]"
layout(matrix(1,1))


## check wingbeats one-by-one:
# for (p in dim(WB_array_fitted)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_accel_fitted', p], typ='l', lwd=2, ylim=c(-2e+08, 2e+08), col='deeppink4', xlab=c('wingbeat cycle [%]'), ylab=c('stroke acceleration [deg/sec^-2]'), bty='n', main=paste0(specimen_now,' - ',sequence_now))
# }






################ ____ deviation angle ###############

p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_fitted', p], typ='l', lwd=2, ylim=c(-100,100), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('deviation angle [deg]'), bty='n', main='deviation angle')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}

layout(matrix(1:2,ncol=2), width = c(2,1), height = c(1,1))
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_fitted', p], typ='l', lwd=2, ylim=c(-100,100), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('deviation angle [deg]'), bty='n', main='deviation angle')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
legend_image <- as.raster(matrix(t(length(seq_now)), ncol=1))
plot(c(0,2), c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'WB frequency [Hz]', cex.main = 1)
text(x=1.4, y = seq(0,1,l=5), labels = seq(100, 300, l=5))
rasterImage(legend_image[length(seq_now):1], 0.3, 0, 0.7, 1) # inverse color direction here with "[length(seq_now):1]"
layout(matrix(1,1))

## check wingbeats one-by-one:
# for (p in dim(WB_array_fitted)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_fitted', p], typ='l', lwd=2, ylim=c(-50,50), col='darkgreen', xlab=c('wingbeat cycle [%]'), ylab=c('deviation angle [deg]'), bty='n', main=paste0(specimen_now,' - ',sequence_now))
#   points(as.numeric(WB_array_norma[(X+1):(norma_wb_length[p]-X),col_wingbeat_cycle, p]), as.numeric(WB_array_norma[(X+1):(norma_wb_length[p]-X), col_deviation_mean, p]), pch=16, cex=1, col='darkred')
# }


################ ____ deviation rate ########################

## just like the stroke rate, the deviation rate is the rate of change of the deviation angle.
## expressed in radian/sec or in deg/sec.

p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_rate_fitted', p], typ='l', lwd=2, ylim=c(-1e+05, 1e+05), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('deviation rate [deg/sec]'), bty='n', main='deviation rate')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_rate_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}

layout(matrix(1:2,ncol=2), width = c(2,1), height = c(1,1))
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_rate_fitted', p], typ='l', lwd=2, ylim=c(-1e+05, 1e+05), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('deviation rate [deg/sec]'), bty='n', main='deviation rate')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_rate_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
legend_image <- as.raster(matrix(t(length(seq_now)), ncol=1))
plot(c(0,2), c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'WB frequency [Hz]', cex.main = 1)
text(x=1.4, y = seq(0,1,l=5), labels = seq(100, 300, l=5))
rasterImage(legend_image[length(seq_now):1], 0.3, 0, 0.7, 1) # inverse color direction here with "[length(seq_now):1]"
layout(matrix(1,1))


## check wingbeats one-by-one:
# for (p in dim(WB_array_fitted)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_rate_fitted', p], typ='l', lwd=2, ylim=c(-60000,60000), col='seagreen4', xlab=c('wingbeat cycle [%]'), ylab=c('deviation rate [deg/sec]'), bty='n', main=paste0(specimen_now,' - ',sequence_now))
# }




################ ____ deviation acceleration ########################


p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_accel_fitted', p], typ='l', lwd=2, ylim=c(-2e+08, 2e+08), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('deviation acceleration [deg/sec^-2]'), bty='n', main='deviation acceleration')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_accel_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}

layout(matrix(1:2,ncol=2), width = c(2,1), height = c(1,1))
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_accel_fitted', p], typ='l', lwd=2, ylim=c(-2e+08, 2e+08), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('deviation acceleration [deg/sec^-2]'), bty='n', main='deviation acceleration')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_accel_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
legend_image <- as.raster(matrix(t(length(seq_now)), ncol=1))
plot(c(0,2), c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'WB frequency [Hz]', cex.main = 1)
text(x=1.4, y = seq(0,1,l=5), labels = seq(100, 300, l=5))
rasterImage(legend_image[length(seq_now):1], 0.3, 0, 0.7, 1) # inverse color direction here with "[length(seq_now):1]"
layout(matrix(1,1))


## check wingbeats one-by-one:
# for (p in dim(WB_array_fitted)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_accel_fitted', p], typ='l', lwd=2, ylim=c(-5000,5000), col='turquoise4', xlab=c('wingbeat cycle [%]'), ylab=c('deviation acceleration [deg/sec^-2]'), bty='n', main=paste0(specimen_now,' - ',sequence_now))
# }






################ ____ rotation angle ###############

p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_fitted', p], typ='l', lwd=2, ylim=c(-100,100), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('rotation angle [deg]'), bty='n', main='rotation angle')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}

layout(matrix(1:2,ncol=2), width = c(2,1), height = c(1,1))
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_fitted', p], typ='l', lwd=2, ylim=c(-100,100), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('rotation angle [deg]'), bty='n', main='rotation angle')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
legend_image <- as.raster(matrix(t(length(seq_now)), ncol=1))
plot(c(0,2), c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'WB frequency [Hz]', cex.main = 1)
text(x=1.4, y = seq(0,1,l=5), labels = seq(100, 300, l=5))
rasterImage(legend_image[length(seq_now):1], 0.3, 0, 0.7, 1) # inverse color direction here with "[length(seq_now):1]"
layout(matrix(1,1))



## check wingbeats one-by-one:
# for (p in dim(WB_array_fitted)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_fitted', p], typ='l', lwd=2, ylim=c(-100,100), col='darkred', xlab=c('wingbeat cycle [%]'), ylab=c('rotation angle [deg]'), bty='n', main=paste0(specimen_now,' - ',sequence_now))
#   points(as.numeric(WB_array_norma[(X+1):(norma_wb_length[p]-X),col_wingbeat_cycle, p]), as.numeric(WB_array_norma[(X+1):(norma_wb_length[p]-X),col_rotation_mean, p]), pch=16, cex=1, col='black')
# }
# 


################ ____ rotation rate ########################

## the rotation rate is the rate of change of the rotation angle.
## expressed in radian/sec or in deg/sec.

p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_rate_fitted', p], typ='l', lwd=2, ylim=c(-1e+05, 1e+05), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('rotation rate [deg/sec]'), bty='n', main='rotation rate')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_rate_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}

layout(matrix(1:2,ncol=2), width = c(2,1), height = c(1,1))
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_rate_fitted', p], typ='l', lwd=2, ylim=c(-1e+05, 1e+05), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('rotation rate [deg/sec]'), bty='n', main='rotation rate')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_rate_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
legend_image <- as.raster(matrix(t(length(seq_now)), ncol=1))
plot(c(0,2), c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'WB frequency [Hz]', cex.main = 1)
text(x=1.4, y = seq(0,1,l=5), labels = seq(100, 300, l=5))
rasterImage(legend_image[length(seq_now):1], 0.3, 0, 0.7, 1) # inverse color direction here with "[length(seq_now):1]"
layout(matrix(1,1))


## check wingbeats one-by-one:
# for (p in dim(WB_array_fitted)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_rate_fitted', p], typ='l', lwd=2, ylim=c(-1000,1000), col='violetred4', xlab=c('wingbeat cycle [%]'), ylab=c('rotation rate [deg/sec]'), bty='n', main=paste0(specimen_now,' - ',sequence_now))
# }



################ ____ rotation acceleration ########################

p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_accel_fitted', p], typ='l', lwd=2, ylim=c(-2e+08, 2e+08), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('rotation acceleration [deg/sec^-2]'), bty='n', main='rotation acceleration')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_accel_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}

layout(matrix(1:2,ncol=2), width = c(2,1), height = c(1,1))
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_accel_fitted', p], typ='l', lwd=2, ylim=c(-2e+08, 2e+08), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('rotation acceleration [deg/sec^-2]'), bty='n', main='rotation acceleration')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_accel_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
legend_image <- as.raster(matrix(t(length(seq_now)), ncol=1))
plot(c(0,2), c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'WB frequency [Hz]', cex.main = 1)
text(x=1.4, y = seq(0,1,l=5), labels = seq(100, 300, l=5))
rasterImage(legend_image[length(seq_now):1], 0.3, 0, 0.7, 1) # inverse color direction here with "[length(seq_now):1]"
layout(matrix(1,1))



# # check wingbeats one-by-one:
# for (p in dim(WB_array_fitted)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_accel_fitted', p], typ='l', lwd=2, ylim=c(-18000,18000), col='purple4', xlab=c('wingbeat cycle [%]'), ylab=c('rotation acceleration [deg/sec^-2]'), bty='n', main=paste0(specimen_now,' - ',sequence_now))
# }







######################## ________ transform rotation angle into AoA ########################

## 0deg rotation angle is 90deg AoA,
## this correspond to the situation where the wing is perpendicular to the body axis

## when rot angle > 0, rot + AoA = 90deg, so  AoA =  90deg - rot
## when rot angle < 0, AoA - rot = 90deg, so AoA = 90deg + rot

for (i in 1:dim(WB_array_fitted)[3]){
  for (j in 1:nrow(WB_array_fitted) ){
    if (WB_array_fitted[j,'rotation_fitted',i] > 0) {WB_array_fitted[j,'angle_of_attack',i] = 90 - as.numeric(WB_array_fitted[j,'rotation_fitted',i])}
    if (WB_array_fitted[j,'rotation_fitted',i] < 0) {WB_array_fitted[j,'angle_of_attack',i] = 90 + as.numeric(WB_array_fitted[j,'rotation_fitted',i])}
  }}


######### plot angle-of-attack

p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'angle_of_attack', p], typ='l', lwd=2, ylim=c(0,100), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('angle-of-attack [deg]'), bty='n', main='angle-of-attack')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'angle_of_attack', p], typ='l', lwd=2, col=gradient_now[p])
}

layout(matrix(1:2,ncol=2), width = c(2,1), height = c(1,1))
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'angle_of_attack', p], typ='l', lwd=2, ylim=c(0,100), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('angle-of-attack [deg]'), bty='n', main='angle-of-attack')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'angle_of_attack', p], typ='l', lwd=2, col=gradient_now[p])
}
legend_image <- as.raster(matrix(t(length(seq_now)), ncol=1))
plot(c(0,2), c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'WB frequency [Hz]', cex.main = 1)
text(x=1.4, y = seq(0,1,l=5), labels = seq(100, 300, l=5))
rasterImage(legend_image[length(seq_now):1], 0.3, 0, 0.7, 1) # inverse color direction here with "[length(seq_now):1]"
layout(matrix(1,1))

## to check wingbeats one-by-one:
# for (p in dim(WB_array_fitted)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'angle_of_attack', p], typ='l', lwd=2, ylim=c(0,100), col='darkred', xlab=c('wingbeat cycle [%]'), ylab=c('angle-of-attack [deg]'), bty='n', main=paste0(specimen_now,' - ',sequence_now))
# }



################ ____ compute "true" AoA ################

## the angle of attack computed from the rotation angle is not the most accurate as it does not consider wing velocity vector.
## Ilam Govindasamy gave me the code below computing properly the AoA, using the
## stroke rate, deviation angle, rotation angle, and deviation rate
## (the code was translated from formulas in "Wang et al. (2016). A predictive quasi-steady model of aerodynamic loads on flapping wings. Journal of Fluid Mechanics, 800, 688-719")

omega_y = array(NA, dim = c(nrow(WB_array_fitted), 1, dim(WB_array_fitted)[3]))
omega_z = array(NA, dim = c(nrow(WB_array_fitted), 1, dim(WB_array_fitted)[3]))

for (i in 1:dim(WB_array_fitted)[3]){
  stroke_rate_now = WB_array_fitted[,'stroke_rate_fitted',i]
  deviation_rate_now = WB_array_fitted[,'deviation_rate_fitted',i]
  deviation_angle_now = WB_array_fitted[,'deviation_fitted',i]
  rotation_angle_now = WB_array_fitted[,'rotation_fitted',i]
  
  omega_y[,1,i] = -stroke_rate_now * cos(deviation_angle_now * pi/180) * cos(rotation_angle_now * pi/180) + deviation_rate_now * sin(rotation_angle_now * pi/180)
  omega_z[,1,i] = stroke_rate_now * cos(deviation_angle_now * pi/180) * sin(rotation_angle_now * pi/180) + deviation_rate_now * cos(rotation_angle_now * pi/180)
  
  WB_array_fitted[,'angle_of_attack_2',i] = abs(atan2(-omega_y[,1,i], omega_z[,1,i])*180/pi)
}


######### plot the true angle of attack

p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'angle_of_attack_2', p], typ='l', lwd=2, ylim=c(0,150), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('angle-of-attack [deg]'), bty='n', main='angle-of-attack')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'angle_of_attack_2', p], typ='l', lwd=2, col=gradient_now[p])
}

layout(matrix(1:2,ncol=2), width = c(2,1), height = c(1,1))
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'angle_of_attack_2', p], typ='l', lwd=2, ylim=c(0,150), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('angle-of-attack [deg]'), bty='n', main='angle-of-attack')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'angle_of_attack_2', p], typ='l', lwd=2, col=gradient_now[p])
}
legend_image <- as.raster(matrix(t(length(seq_now)), ncol=1))
plot(c(0,2), c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'WB frequency [Hz]', cex.main = 1)
text(x=1.4, y = seq(0,1,l=5), labels = seq(50, 600, l=5))
rasterImage(legend_image[length(seq_now):1], 0.3, 0, 0.7, 1) # inverse color direction here with "[length(seq_now):1]"
layout(matrix(1,1))

## to check wingbeats one-by-one:
# for (p in dim(WB_array_fitted)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'angle_of_attack_2', p], typ='l', lwd=2, ylim=c(0,150), col='darkred', xlab=c('wingbeat cycle [%]'), ylab=c('angle-of-attack [deg]'), bty='n', main=paste0(specimen_now,' - ',sequence_now))
# }



##### compute and plot the average angle of attack at each frame (mean of all wingbeats)

average_AoA = as.data.frame(matrix(NA,nrow(WB_array_fitted)[1], 2))
colnames(average_AoA) <-  c('wingbeat_cycle', 'average_AoA')
average_AoA$wingbeat_cycle = as.numeric(WB_array_fitted[,'wingbeat_cycle', 1])
for (i in 1:nrow(average_AoA)){
  average_AoA$average_AoA[i] = mean(WB_array_fitted[i,'angle_of_attack_2',])
}

## plot angle of attack and its average
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'angle_of_attack_2', p], typ='l', lwd=2, ylim=c(0,180), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('angle-of-attack [deg]'), bty='n', main='angle-of-attack')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'angle_of_attack_2', p], typ='l', lwd=2, col=gradient_now[p])
}
points(average_AoA$wingbeat_cycle, average_AoA$average_AoA, typ='l', lwd=4, col='black')







################ ____ absolute wing speed ################

## we integrate both the stroke and the deviation movement to obtain the absolute wing speed (omega).
## following pythagore: omega = sqrt((stroke rate)^2 + (deviation rate)^2) 

for (i in 1:dim(WB_array_fitted)[3]){
  for (j in 1:nrow(WB_array_fitted)){
    WB_array_fitted[j, 'abs_wing_speed', i] = sqrt( as.numeric(WB_array_fitted[j,'stroke_rate_fitted',i])^2 + as.numeric(WB_array_fitted[j,'deviation_rate_fitted',i])^2)
  }
}

### plot absolute wing speed:
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'abs_wing_speed', p], typ='l', lwd=2, ylim=c(0,1e+05), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('absolute wing speed [deg/sec]'), bty='n', main='absolute wing speed')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'abs_wing_speed', p], typ='l', lwd=2, col=gradient_now[p])
}

layout(matrix(1:2,ncol=2), width = c(2,1), height = c(1,1))
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'abs_wing_speed', p], typ='l', lwd=2, ylim=c(0,1e+05), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('absolute wing speed [deg/sec]'), bty='n', main='absolute wing speed')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'abs_wing_speed', p], typ='l', lwd=2, col=gradient_now[p])
}
legend_image <- as.raster(matrix(t(length(seq_now)), ncol=1))
plot(c(0,2), c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'WB frequency [Hz]', cex.main = 1)
text(x=1.4, y = seq(0,1,l=5), labels = seq(100, 300, l=5))
rasterImage(legend_image[length(seq_now):1], 0.3, 0, 0.7, 1) # inverse color direction here with "[length(seq_now):1]"
layout(matrix(1,1))









############################# [2.5] grand average fourier fit [family-level] ############################# 

## load the file generated using the MATLAB code
## 'AA_batch_fit_fourier_wb_averaged.m'

## **the order of magnitude of the fourier series is for now 4, we could debate on whether the 2, 3, 4 5 or 6 are better fitting the data...
## 4th order of magnitude for all three angles (444):
AVERAGE_fitted_wb_from_MATLAB <- readMat(paste0(local_path,'/02_fit_fourier_series_hoverflies/data/wb_averaged_family_level/tracked_and_fourier_fits_all_wingbeats_averaged.mat'), h=T, sep=';')


## the list 'AVERAGE_fitted_wb_from_MATLAB' contains the following elements:
##  - str/dev/rot.tracked are the tracked (raw) angle values. each column is a wingbeat.
##    the length of each tracked wingbeat is different.
##  - time.tracked contains the corresponding time columns associated with the tracked values.
##
##  - str/dev/rot.fit are the fitted ('smoothed') angle values. each column is a wingbeat.
##    they all have the same length of 101 time step.
##  - time.fit contains the corresponding time columns associated with the fitted values (all of length 101).
##
##  - str/dev/rot.fitresults is a list of list. each list in the list contains the model associated with each curve fitting. 
##    apparently it cannot be open in R. see the MATLAB object directly in MATLAB 
##  - str/dev/rot.gofs is a list of list. each list in the list contains the goodness of fit (rmse, r-squared etc)
##    associated with the model. it can be accessed using for example 'AVERAGE_fitted_wb_from_MATLAB$dev.gofs[[1]]' for the first wb
##  - str/dev/rot.fcoeffs contains the 9 coefficients (a0,a1,a2,a3,a4 and b1,b2,b3,b4) of the model
##    it is a vector of 9 values for each wingbeats
##
##  - str/dev/rot.dot.fit contains the derivative of the fitted values, i.e. the stroke/dev/rot rate (or velocity)
##  - str/dev/rot.dotdot.fit contains the 2nd derivative of the fitted values, i.e. the stroke/dev/rot acceleration

### create a dataset compiling the fitted data (we don't need an array because there is only one average wingbeat):
AVERAGE_WB = as.data.frame(matrix(NA,nrow(AVERAGE_fitted_wb_from_MATLAB$time.fit), 13))
colnames(AVERAGE_WB) <-  c('wingbeat_cycle', 'stroke_fitted', 'deviation_fitted', 'rotation_fitted',
                           'stroke_rate_fitted', 'deviation_rate_fitted', 'rotation_rate_fitted',
                           'stroke_accel_fitted', 'deviation_accel_fitted', 'rotation_accel_fitted',
                           'angle_of_attack','angle_of_attack_2','abs_wing_speed')

AVERAGE_WB[ , 'wingbeat_cycle'] = as.numeric(AVERAGE_fitted_wb_from_MATLAB$time.fit*100)
AVERAGE_WB[ , 'stroke_fitted'] = AVERAGE_fitted_wb_from_MATLAB$str.fit[,1]
AVERAGE_WB[ , 'deviation_fitted'] = AVERAGE_fitted_wb_from_MATLAB$dev.fit[,1]
AVERAGE_WB[ , 'rotation_fitted'] = AVERAGE_fitted_wb_from_MATLAB$rot.fit[,1]
AVERAGE_WB[ , 'stroke_rate_fitted'] = AVERAGE_fitted_wb_from_MATLAB$str.dot.fit[,1]
AVERAGE_WB[ , 'deviation_rate_fitted'] = AVERAGE_fitted_wb_from_MATLAB$dev.dot.fit[,1]
AVERAGE_WB[ , 'rotation_rate_fitted'] = AVERAGE_fitted_wb_from_MATLAB$rot.dot.fit[,1]
AVERAGE_WB[ , 'stroke_accel_fitted'] = AVERAGE_fitted_wb_from_MATLAB$str.dotdot.fit[,1]
AVERAGE_WB[ , 'deviation_accel_fitted'] = AVERAGE_fitted_wb_from_MATLAB$dev.dotdot.fit[,1]
AVERAGE_WB[ , 'rotation_accel_fitted'] = AVERAGE_fitted_wb_from_MATLAB$rot.dotdot.fit[,1]

# View(AVERAGE_WB[,])


#### plot the average fitting

## stroke angle
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_fitted', p], typ='l', lwd=2, ylim=c(-100,100), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('stroke angle [deg]'), bty='n', main='stroke angle')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
points(AVERAGE_WB[,'wingbeat_cycle'], AVERAGE_WB[,'stroke_fitted'], typ='l', lwd=4, col='black')

## deviation angle
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_fitted', p], typ='l', lwd=2, ylim=c(-100,100), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('deviation angle [deg]'), bty='n', main='deviation angle')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
points(AVERAGE_WB[,'wingbeat_cycle'], AVERAGE_WB[,'deviation_fitted'], typ='l', lwd=4, col='black')

## rotation angle
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_fitted', p], typ='l', lwd=2, ylim=c(-100,100), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('rotation angle [deg]'), bty='n', main='rotation angle')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
points(AVERAGE_WB[,'wingbeat_cycle'], AVERAGE_WB[,'rotation_fitted'], typ='l', lwd=4, col='black')



## for some unknown reason, the angle rate can be fitted if multiplied by 160:
AVERAGE_WB[,'stroke_rate_fitted'] <- AVERAGE_WB[,'stroke_rate_fitted']*160
AVERAGE_WB[,'deviation_rate_fitted'] <- AVERAGE_WB[,'deviation_rate_fitted']*160
AVERAGE_WB[,'rotation_rate_fitted'] <- AVERAGE_WB[,'rotation_rate_fitted']*160

## stroke rate
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_rate_fitted', p], typ='l', lwd=2, ylim=c(-1e+05,1e+05), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('stroke rate [deg/sec]'), bty='n', main='stroke rate')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_rate_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
points(AVERAGE_WB[,'wingbeat_cycle'], AVERAGE_WB[,'stroke_rate_fitted'], typ='l', lwd=4, col='black')

## deviation rate
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_rate_fitted', p], typ='l', lwd=2, ylim=c(-1e+05,1e+05), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('deviation rate [deg/sec]'), bty='n', main='deviation rate')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'deviation_rate_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
points(AVERAGE_WB[,'wingbeat_cycle'], AVERAGE_WB[,'deviation_rate_fitted'], typ='l', lwd=4, col='black')

## rotation rate
p=1
plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_rate_fitted', p], typ='l', lwd=2, ylim=c(-1e+05,1e+05), col=gradient_now[p], xlab=c('wingbeat cycle [%]'), ylab=c('rotation rate [deg/sec]'), bty='n', main='rotation rate')
for (p in 2:dim(WB_array_fitted)[3]){
  points(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'rotation_rate_fitted', p], typ='l', lwd=2, col=gradient_now[p])
}
points(AVERAGE_WB[,'wingbeat_cycle'], AVERAGE_WB[,'rotation_rate_fitted'], typ='l', lwd=4, col='black')




## save the average wingbeat kinematics of all studied hoverflies to send to Thomas Engels
# write.csv(AVERAGE_WB, paste0(local_path, '/data sent to Thomas Engels/average_wb_kinematics_hoverflies.csv'), row.names=F)


############# ___ plot superimposed all angles ############# 

plot(WB_array_fitted[,'wingbeat_cycle', p], WB_array_fitted[,'stroke_fitted', p], typ='l', col='#DD6D2F', lwd=4, ylim=c(-60,60), xlab=c('wingbeat cycle [%]'), ylab=c('angle [deg]'), bty='n')
points(AVERAGE_WB[,'wingbeat_cycle'], AVERAGE_WB[,'deviation_fitted'], typ='l', lwd=4, col='#387EB9')
points(AVERAGE_WB[,'wingbeat_cycle'], AVERAGE_WB[,'rotation_fitted'], typ='l', lwd=4, col='#4DAF49')





############################# [3] body dynamic ############################# 

differentiate_central <- function(x, dt){
  seq_now = which(is.na(x)==F)
  dx = matrix(NA,length(x),1)
  if (length(x) > 2) {
    dx[2:(length(dx)-1),] = (x[3:length(x)] - x[1:(length(x)-2)]) /2/dt
    dx[(seq_now[1]+1):(last(seq_now)-1),] = (x[(seq_now[1]+2):last(seq_now)] - x[seq_now[1]:(last(seq_now)-2)]) /2/dt
    dx[seq_now[1],] = (-3*x[seq_now[1]]+4*x[(seq_now[1]+1)]+x[(seq_now[1]+2)]) /2/dt
    dx[last(seq_now),] = (3*x[last(seq_now)] - 4*x[(last(seq_now)-1)] + x[(last(seq_now)-2)]) /2/dt 
  } else {
    dx[1,] = (x[2] - x[1]) /dt
    dx[length(x)] = (x[length(x)] - x[(length(x)-1)]) /dt
  }
  return(dx)
}

################ ____ body speed ############### 

# the instantaneous speed is the norm (or 'magnitude') of the velocity vector
# function to compute the norm of a vector
norma <- function(x){
  sqrt(x[1]^2 + x[2]^2 + x[3]^2)
}


## compute velocity in the x, y and z directions (u,v,w) using the "differentiate_central" function:
velocity_vectors = array(NA, dim = c(dim(WB_array)[1],3,dim(WB_array)[3]))# to store the computed velocity vectors
dimnames(velocity_vectors)[[2]] = c('u','v','w')
for (i in 1:dim(WB_array)[3]){
  ## create empty vectors to store the x,y and z component of velocity
  X_velocity_now = Y_velocity_now = Z_velocity_now = rep(NA, length(first_digitized_frame[i]:last_digitized_frame[i]))
  ## dt now
  dt_now = as.numeric(WB_array[first_digitized_frame[i],col_time,i])
  ## derive each components
  X_velocity_now = differentiate_central(as.numeric(WB_array[,col_x,i]), dt_now)
  Y_velocity_now = differentiate_central(as.numeric(WB_array[,col_y,i]), dt_now)
  Z_velocity_now = differentiate_central(as.numeric(WB_array[,col_z,i]), dt_now)
  velocity_vector_now = cbind(X_velocity_now, Y_velocity_now, Z_velocity_now)
  velocity_vectors[,,i] = velocity_vector_now
  ## compute instantaneous speed as the norm of the velocity vector (i.e. binded x,y,z)
  for (j in first_digitized_frame[i]:last_digitized_frame[i]){
    WB_array[j,col_body_speed,i] = norma(velocity_vector_now[j,])
    # transfer the body speed values in the normalized array:
    WB_array_norma[1:norma_wb_length[i], col_body_speed, i] = WB_array[(start_end_wb$start[i]) : (start_end_wb$end[i]), col_body_speed, i]
  }
}

## correct for weird jump in body speed at frame 14 in WB_array_norma[,col_body_speed,13] (N_008 seq_003)
WB_array_norma[which(WB_array_norma[,col_body_speed,which(full_IDs=="bodyNwing_kine_seq_003_N_008.mat")]>10),col_body_speed,which(full_IDs=="bodyNwing_kine_seq_003_N_008.mat")] = mean(c(0.255652328243008,0.141832589148993))


## plot variation in body speed over the wingbeat:
p=1
plot(WB_array_norma[,col_wingbeat_cycle, p], WB_array_norma[,col_body_speed, p], typ='l', lwd=2, ylim=c(0,2.5), col='darkblue', xlab=c('wingbeat cycle [%]'), ylab=c('body speed [m.s-1]'), bty='n', main='body speed')
for (p in 2:dim(WB_array_norma)[3]){
  points(WB_array_norma[,col_wingbeat_cycle, p], WB_array_norma[,col_body_speed, p], typ='l', lwd=2, col='darkblue')
}

## to check the wingbeats one by one:
# for (p in dim(WB_array_norma)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_norma[,col_wingbeat_cycle, p], WB_array_norma[,col_body_speed, p], typ='l', lwd=2, ylim=c(0,3), col='darkblue', xlab=c('wingbeat cycle [%]'), ylab=c('body speed [m.s-1]'), bty='n', main=paste0('body speed: ',specimen_now,' - ',sequence_now))
# }



################ ____ body acceleration ############### 

## compute acceleration using the "differentiate_central" function on the previously computed body speed
for (i in 1:dim(WB_array)[3]){
  dt_now = as.numeric(WB_array[first_digitized_frame[i],col_time,i])
  WB_array[,col_body_accel,i] = differentiate_central(as.numeric(WB_array[,col_body_speed,i]), dt_now)
  # transfer the body acceleration values in the normalized array:
  WB_array_norma[1:norma_wb_length[i], col_body_accel, i] = WB_array[(start_end_wb$start[i]) : (start_end_wb$end[i]), col_body_accel, i]
}

## plot variation in body acceleration over the wingbeat:
p=1
plot(WB_array_norma[,col_wingbeat_cycle, p], WB_array_norma[,col_body_accel, p], typ='l', ylim=c(-2000,2000), lwd=2, col='darkred', xlab=c('wingbeat cycle [%]'), ylab=c('body acceleration [m.s-2]'), bty='n', main='body acceleration')
for (p in 2:dim(WB_array_norma)[3]){
  points(WB_array_norma[,col_wingbeat_cycle, p], WB_array_norma[,col_body_accel, p], typ='l', lwd=2, col='darkred')
}

## to check the wingbeats one by one:
# for (p in dim(WB_array_norma)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_norma[,col_wingbeat_cycle, p], WB_array_norma[,col_body_accel, p], typ='l', lwd=2, ylim=c(-10000,10000), col='darkblue', xlab=c('wingbeat cycle [%]'), ylab=c('body acceleration [m.s-2]'), bty='n', main=paste0('body acceleration: ',specimen_now,' - ',sequence_now))
# }

# ACCELERATION IS EXTREMELY NOISY






################ ____ climb angle ############### 

## climb angle = 180/pi*atan2(w, sqrt( u.^2 + v.^2))
## see 'Muijres et al. 2017. Escaping blood-fed malaria mosquitoes minimize tactile detection without compromising on take-off speed'

## computing climb angle
for (i in 1:dim(WB_array)[3]){
  for (j in first_digitized_frame[i]:last_digitized_frame[i]){
    WB_array[j,col_climb_angle,i] = 180/pi*atan2(velocity_vectors[j,,i][3], sqrt( (velocity_vectors[j,,i][1])^2 + (velocity_vectors[j,,i][2])^2))
  }
  # transfer the climb angle values in the normalized array:
  WB_array_norma[1:norma_wb_length[i], col_climb_angle, i] = WB_array[(start_end_wb$start[i]) : (start_end_wb$end[i]), col_climb_angle, i]
}

## plot variation in climb angle over the wingbeat:
p=1
plot(WB_array_norma[,col_wingbeat_cycle, p], WB_array_norma[,col_climb_angle, p], typ='l', ylim=c(-100,100), lwd=2, col='darkorchid4', xlab=c('wingbeat cycle [%]'), ylab=c('climb angle [deg]'), bty='n', main='climb angle')
for (p in 2:dim(WB_array_norma)[3]){
  points(WB_array_norma[,col_wingbeat_cycle, p], WB_array_norma[,col_climb_angle, p], typ='l', lwd=2, col='darkorchid4')
}

## to check the wingbeats one by one:
# for (p in dim(WB_array_norma)[3]:1){
#   specimen_now = substr(full_IDs[p],24,28)
#   sequence_now = substr(full_IDs[p],16,22)
#   plot(WB_array_norma[,col_wingbeat_cycle, p], WB_array_norma[,col_climb_angle, p], typ='l', lwd=2, ylim=c(-100,100), col='darkorchid4', xlab=c('wingbeat cycle [%]'), ylab=c('climb angle [deg]'), bty='n', main=paste0('climb angle: ',specimen_now,' - ',sequence_now))
# }




################ ____ covered distance ############### 

## compute the consecutive distance between all digitized positions
# consec_dist = array(NA, dim= c(dim(WB_array)[1], 1, dim(WB_array)[3]))
# prog_bar <- txtProgressBar(min = 0, max = dim(WB_array)[3], style = 3)
# for (i in 1:dim(WB_array)[3] ){
#   for (j in (1+X):(norma_wb_length[i]+X)){
#     consec_dist[j,1,i] = sqrt( (diff(as.numeric(na.omit(WB_array_norma[,col_x,i])))[j])^2 + 
#                                  (diff(as.numeric(na.omit(WB_array_norma[,col_y,i])))[j])^2 +
#                                  (diff(as.numeric(na.omit(WB_array_norma[,col_z,i])))[j])^2 )
#   }
#   Sys.sleep(0.1)
#   setTxtProgressBar(prog_bar, i)
# }
# close(prog_bar)
# 
# covered_dist = matrix(NA, dim(WB_array)[3], 1)
# for (i in 1:dim(WB_array)[3]){
#   covered_dist[i,] = sum(na.omit(consec_dist[,,i]))
# }








############################# [4] compute wingbeat parameters ############################# 

## we create a new dataframe in which we store the (mean) value over the wingbeat for each parameter
wb_parameters = as.data.frame(matrix(NA, length(full_IDs), 24))
colnames(wb_parameters) = c('full_ID', 'species_ID', 'sequence_ID', 'genus', 'species', 'genus_species',
                            'WB_freq', 'stroke_amplitude', 'rotation_amplitude','deviation_amplitude',
                            'peak_ang_velo_fs', 'peak_ang_velo_bs' ,'AoA_peak_ang_velo_fs', 'AoA_peak_ang_velo_bs',
                            'peak_dev_velo', 'peak_rot_velo', 'peak_ang_accel', 'peak_dev_accel', 'peak_rot_accel',
                            'abs_wing_speed', 'mean_body_speed', 'mean_body_accel', 'mean_climb_angle', 'mean_body_pitch'
)
# *** --> "species_ID" was initially called specimen_ID but the specimen level is actually hidden flight trials variation
wb_parameters$full_ID = as.factor(full_IDs)
wb_parameters$species_ID = as.factor(substr(full_IDs,24,28))
wb_parameters$sequence_ID = as.factor(substr(full_IDs,16,22))
wb_parameters$genus = as.factor(IDs_file_seq_lvl$genus)
wb_parameters$species = as.factor(IDs_file_seq_lvl$species)
for (i in 1:length(full_IDs)){wb_parameters$genus_species[i] = paste0(wb_parameters$genus[i],'_',wb_parameters$species[i])}
wb_parameters$genus_species = as.factor(wb_parameters$genus_species)








################ ____ wingbeat frequency ################ 

wb_parameters$WB_freq = IDs_file_seq_lvl$full_seq_WB_freq
# this freq is more accurate because computed on the full sequence


# for (i in 1:nrow(wb_parameters)){
#   time_step_now = as.numeric(WB_array_norma[2,col_time,i]) - as.numeric(WB_array_norma[1,col_time,i])
#   frame_rate_now = 1/time_step_now
#   # wingbeat frequency is computed as 1/time needed to complete a full wingbeat
#   wb_parameters$WB_freq[i] = 1 / ((norma_wb_length[i]-2*X) * time_step_now)
# }





################ ____ stroke parameters ################ 

######## stroke amplitude
## stroke amplitude is computed as the absolute value of (min_value - max_value) 
for (i in 1:nrow(wb_parameters)){
  wb_parameters$stroke_amplitude[i] = abs( max(as.numeric(na.omit(WB_array_fitted[,'stroke_fitted',i]))) - min(as.numeric(na.omit(WB_array_fitted[,'stroke_fitted',i]))) )
}

mean(wb_parameters$stroke_amplitude)
sd(wb_parameters$stroke_amplitude)

######## peaks stroke rate (dot)
## (1) first method: just take the peak value ( max(stroke_rate) )
# for (i in 1:nrow(wb_parameters)){
#   wb_parameters$peak_ang_velo_fs[i] = max(abs(as.numeric(na.omit(WB_array_fitted[1:51,'stroke_rate_fitted',i]))))
#   wb_parameters$peak_ang_velo_bs[i] = max(abs(as.numeric(na.omit(WB_array_fitted[51:101,'stroke_rate_fitted',i]))))
# }
## (2) second method:
## take the average value from 15% to 25% of the wingbeat (for front-stroke)
## and average value from 75% to 85% of the wingbeat (for front-stroke) I think this method damp the noise
for (i in 1:nrow(wb_parameters)){
  wb_parameters$peak_ang_velo_fs[i] = mean(abs(as.numeric(na.omit(WB_array_fitted[16:26,'stroke_rate_fitted',i]))))
  wb_parameters$peak_ang_velo_bs[i] = max(abs(as.numeric(na.omit(WB_array_fitted[76:86,'stroke_rate_fitted',i]))))
}

######## peak stroke acceleration (dot dot)
# just take the peak value
for (i in 1:nrow(wb_parameters)){
  wb_parameters$peak_ang_accel[i] = max(as.numeric(na.omit(WB_array_fitted[,'stroke_accel_fitted',i])))
}



################ ____ deviation parameters ################ 

######## deviation amplitude
## deviation amplitude is computed as the absolute value of (min_value - max_value) 
for (i in 1:nrow(wb_parameters)){
  wb_parameters$deviation_amplitude[i] = abs( max(as.numeric(na.omit(WB_array_fitted[,'deviation_fitted',i]))) - min(as.numeric(na.omit(WB_array_fitted[,'deviation_fitted',i]))) )
}

######## peak deviation rate (dot)
for (i in 1:nrow(wb_parameters)){
  wb_parameters$peak_dev_velo[i] = max(as.numeric(na.omit(WB_array_fitted[,'deviation_rate_fitted',i])))
}
######## peak deviation acceleration (dot dot)
for (i in 1:nrow(wb_parameters)){
  wb_parameters$peak_dev_accel[i] = max(abs(as.numeric(na.omit(WB_array_fitted[,'deviation_accel_fitted',i]))))
}




################ ____ rotation parameters ################ 

######## rotation amplitude
## rotation amplitude is computed as the absolute value of (min_value - max_value) 
for (i in 1:nrow(wb_parameters)){
  wb_parameters$rotation_amplitude[i] = abs( max(as.numeric(na.omit(WB_array_fitted[,'rotation_fitted',i]))) - min(as.numeric(na.omit(WB_array_fitted[,'rotation_fitted',i]))) )
}

######## peak rotation rate (dot)
for (i in 1:nrow(wb_parameters)){
  wb_parameters$peak_rot_velo[i] = max(as.numeric(na.omit(WB_array_fitted[,'rotation_rate_fitted',i])))
}
######## peak deviation acceleration (dot dot)
for (i in 1:nrow(wb_parameters)){
  wb_parameters$peak_rot_accel[i] = max(abs(as.numeric(na.omit(WB_array_fitted[,'rotation_accel_fitted',i]))))
}



################ ____ absolute wing speed ################ 
## take the average value from 15% to 25% of the wingbeat (for front-stroke)
## and average value from 75% to 85% of the wingbeat (for front-stroke)
for (i in 1:nrow(wb_parameters)){
  wb_parameters$abs_wing_speed[i] = mean(as.numeric(na.omit(WB_array_fitted[16:26,'abs_wing_speed',i])))
  # wb_parameters$abs_wing_speed[i] = max(as.numeric(na.omit(WB_array_fitted[76:86,'abs_wing_speed',i])))
}



################ ____ AoA at peak wing velocity ################

## (1) first method: just take AoA matching the peak value of angular velocity 
## --> only works if first method is used when computing angular velocity
# for (i in 1:nrow(wb_parameters)){
#   wb_parameters$AoA_peak_ang_velo_fs[i] = as.numeric(WB_array_fitted[which(WB_array_fitted[,'stroke_rate_fitted',i] == wb_parameters$peak_ang_velo_fs[i]), 'angle_of_attack_2', i])
#   wb_parameters$AoA_peak_ang_velo_bs[i] = as.numeric(WB_array_fitted[which(WB_array_fitted[,'stroke_rate_fitted',i] == wb_parameters$peak_ang_velo_bs[i]), 'angle_of_attack_2', i])
# }
## (2) second method:
## take the average AoA from 15% to 25% of the wingbeat (for front-stroke)
## and average AoA from 75% to 85% of the wingbeat (for front-stroke)
for (i in 1:nrow(wb_parameters)){
  wb_parameters$AoA_peak_ang_velo_fs[i] = mean(as.numeric(na.omit(WB_array_fitted[16:26,'angle_of_attack_2',i])))
  wb_parameters$AoA_peak_ang_velo_bs[i] = mean(as.numeric(na.omit(WB_array_fitted[76:86,'angle_of_attack_2',i])))
}




################ ____ body parameters ################

########## mean body speed

## first method, take the mean of instantaneous speeds
for (i in 1:nrow(wb_parameters)){
  wb_parameters$mean_body_speed[i] = mean(as.numeric(na.omit(WB_array_norma[(X+1):(norma_wb_length[i]+X),col_body_speed,i])))
  # the '(X+1):(norma_wb_length[i]-X)' makes sure we are considering the normalized wingbeat, i.e. from pronation to supination
}

## alternative method:
## divide the covered distance by time
# for (i in 1:nrow(wb_parameters)){
#   wb_parameters$mean_body_speed[i] = covered_dist[i] / (as.numeric(WB_array_norma[(norma_wb_length[i]+X),col_time,i]) - as.numeric(WB_array_norma[(X+1),col_time,i]))
#   # the '(X+1):(norma_wb_length[i]-X)' makes sure we are considering the normalized wingbeat, i.e. from pronation to supination
# }



########## mean body acceleration
for (i in 1:nrow(wb_parameters)){
  wb_parameters$mean_body_accel[i] = mean(as.numeric(na.omit(WB_array_norma[(X+1):(norma_wb_length[i]-X),col_body_accel,i])))
  # the '(X+1):(norma_wb_length[i]-X)' makes sure we are considering the normalized wingbeat, i.e. from pronation to supination
} # ACCELERATION IS EXTREMELY NOISY


########## mean climb angle
for (i in 1:nrow(wb_parameters)){
  wb_parameters$mean_climb_angle[i] = mean(as.numeric(na.omit(WB_array_norma[(X+1):(norma_wb_length[i]-X),col_climb_angle,i])))
  # the '(X+1):(norma_wb_length[i]-X)' makes sure we are considering the normalized wingbeat, i.e. from pronation to supination
}

##########  mean body pitch
for (i in 1:nrow(wb_parameters)){
  wb_parameters$mean_body_pitch[i] = mean(as.numeric(na.omit(WB_array_fitted[,"body_pitch" ,i])))
  ## we can then express the stroke plane inclinations as body_pitch - 45? (because we set the stroke plane to be at 45? from the body axis)
  wb_parameters$stroke_plane_incl[i] = wb_parameters$mean_body_pitch[i] - 45
}

mean(wb_parameters$mean_body_pitch)
sd(wb_parameters$mean_body_pitch)

##########  mean stroke incl relative to velocity vect
for (i in 1:nrow(wb_parameters)){
  wb_parameters$stroke_plane_incl_to_velo[i] = wb_parameters$mean_body_pitch[i] - wb_parameters$mean_climb_angle[i] - 45
}

mean(wb_parameters$stroke_plane_incl)
sd(wb_parameters$stroke_plane_incl)





################ ____ advance ratio ################

## in Ellington 1984c, The Aerodynamics of Hovering Insect Flight. III. Kinematics,
## the advance ratio is define as J = V / 2*phi*n*R,
## with V the forward speed, phi the stroke amplitude, n the wingbeat frequency, R the wingspan
## "2*phi*n*R" can be replace by direct measurement of wing angular velocity * R, which we have here.
## Ellington proposes the limit of 0.1 for defining hovering flight
## this paper uses the same criteria for instance: Cheng, Xin, and Mao Sun. "Wing-kinematics measurement and aerodynamics in a small insect in hovering flight." Scientific reports 6.1 (2016): 25706.

## ** we don't have individual level measurements of wing length, but let's use the mean wing length per species
## we take the mean wingspan from the full_dataset created in the code "BBhov_analyzing_WBkinematics_vs_morphology.R"
full_dataset_species_level <- read.csv(paste0(local_path,'/full_dataset_species_level.csv'))

for (i in 1:nrow(wb_parameters)){
  wb_parameters$advance_ratio[i] = wb_parameters$mean_body_speed[i] / (wb_parameters$peak_ang_velo_fs[i]*(pi/180)*(full_dataset_species_level$length_cm[which(full_dataset_species_level$genus_species==wb_parameters$genus_species[i])]/100))
} ## here we express the angular velocity in rad/sec, so we multiply it with (pi/180)
format(wb_parameters$advance_ratio, scientific = F)



mean(wb_parameters$advance_ratio)
sd(wb_parameters$advance_ratio)
hist(wb_parameters$advance_ratio, main='distribution of advance ratios', xlab='advance ratio (J=U/ang_velo_rads)')



# write.csv(wb_parameters, paste0(local_path, '/wb_parameters.csv'), row.names = F)

## locate numerical column in "wb_parameters"
wb_parameters_num_cols <- unlist(lapply(wb_parameters, is.numeric)) # identify numeric columns

## check for correlation among parameters
# chart.Correlation(wb_parameters[,wb_parameters_num_cols], method="pearson", histogram=T, pch=16)






############################# [4.5] effect of body dynamic on WB kinematics [multiple regressions] ############################# 
wb_parameters$wing_speed_ftimesA = wb_parameters$WB_freq * wb_parameters$stroke_amplitude

summary(lm(wb_parameters$wing_speed_ftimesA ~ wb_parameters$mean_body_speed + wb_parameters$mean_climb_angle))

summary(lm(wb_parameters$abs_wing_speed ~ wb_parameters$mean_body_speed + wb_parameters$mean_climb_angle))
summary(lm(wb_parameters$AoA_peak_ang_velo_fs ~ wb_parameters$mean_body_speed + wb_parameters$mean_climb_angle))
summary(lm(wb_parameters$WB_freq ~ wb_parameters$mean_body_speed + wb_parameters$mean_climb_angle))
summary(lm(wb_parameters$stroke_amplitude ~ wb_parameters$mean_body_speed + wb_parameters$mean_climb_angle))
summary(lm(wb_parameters$deviation_amplitude ~ wb_parameters$mean_body_speed + wb_parameters$mean_climb_angle))
summary(lm(wb_parameters$rotation_amplitude ~ wb_parameters$mean_body_speed + wb_parameters$mean_climb_angle))
summary(lm(wb_parameters$peak_ang_velo_fs ~ wb_parameters$mean_body_speed + wb_parameters$mean_climb_angle))
summary(lm(wb_parameters$peak_ang_accel ~ wb_parameters$mean_body_speed + wb_parameters$mean_climb_angle))


## mean overall body speed for the studied wingbeat
mean(wb_parameters$mean_body_speed)
sd(wb_parameters$mean_body_speed)
se<- function(vecteur) { sd(vecteur)/sqrt(length(vecteur)) }
se(wb_parameters$mean_body_speed)

## mean overall climb angle for the studied wingbeat
mean(wb_parameters$mean_climb_angle)
sd(wb_parameters$mean_climb_angle)
se(wb_parameters$mean_climb_angle)



############################# [5] between group variation (boxplots) ############################# 


#### change the order of levels(wb_parameters$genus_species) according to a parameter. for example average weight:
## first get mean weight per species:
mean_value_per_species = as.data.frame(matrix(NA, length(levels(wb_parameters$genus_species)), length(colnames(wb_parameters)[wb_parameters_num_cols])))
colnames(mean_value_per_species) = colnames(wb_parameters)[wb_parameters_num_cols]
rownames(mean_value_per_species) = levels(wb_parameters$genus_species)
for (i in 1:length(colnames(wb_parameters)[wb_parameters_num_cols])){
  mean_value_per_species[,i] = tapply(wb_parameters[,which(colnames(wb_parameters)==colnames(mean_value_per_species)[i])], wb_parameters$genus_species, mean)
}
mean_value_per_species$weight = tapply(IDs_file_seq_lvl$mean_weight_mg, IDs_file_seq_lvl$genus_species, mean)
param_now = mean_value_per_species$weight
wb_parameters$genus_species <- factor(wb_parameters$genus_species, levels = c(levels(wb_parameters$genus_species)[order(param_now)]))


### set discrete color for species weight:
discrete_weight_color <- scico(25, palette = 'bamako')[c(1,3,6,8,10,12,13,15,17,19)] #choose gradient among: scico_palette_show()
discrete_weight_color = discrete_weight_color[length(discrete_weight_color):1] # reverse the gradient
plot(rep(1,length(discrete_weight_color)), pch=15, cex=15,col=discrete_weight_color, bty='n', axes=F)

colo_weight <- discrete_weight_color[match(wb_parameters$genus_species,levels(wb_parameters$genus_species))]



################ ____ inter-specific variation ################

## we want to plot, for each wb parameter, a boxplot per species showing data variation.

###### wingbeat frequency
ggplot(wb_parameters, aes(x = genus_species , y = WB_freq, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'wingbeat frequency [Hz]', x='taxa') + ylim(50,400) + labs(fill = '') +
  ggtitle('variation in wb frequency') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))
anova(lm(wb_parameters$WB_freq ~ wb_parameters$species)) %>% 
  
  ###### stroke amplitude
  ggplot(wb_parameters, aes(x = genus_species , y = stroke_amplitude, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'stroke amplitude [deg]', x='taxa') + ylim(0,250) + labs(fill = '') +
  ggtitle('variation in stroke amplitude') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))
anova(lm(wb_parameters$stroke_amplitude ~ wb_parameters$species))

###### rotation amplitude
ggplot(wb_parameters, aes(x = genus_species , y = rotation_amplitude, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'rotation amplitude [deg]', x='taxa') + ylim(0,250) + labs(fill = '') +
  ggtitle('variation in rotation amplitude') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))
anova(lm(wb_parameters$rotation_amplitude ~ wb_parameters$species))

###### deviation amplitude
ggplot(wb_parameters, aes(x = genus_species , y = deviation_amplitude, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'deviation amplitude [deg]', x='taxa') + ylim(0,50) + labs(fill = '') +
  ggtitle('variation in deviation amplitude') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))
anova(lm(wb_parameters$deviation_amplitude ~ wb_parameters$species))

###### peak stroke rate
ggplot(wb_parameters, aes(x = genus_species , y = peak_ang_velo_fs, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'peak stroke rate [deg/sec]', x='taxa') + ylim(0,160000) + labs(fill = '') +
  ggtitle('variation in peak stroke rate ') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))
anova(lm(wb_parameters$peak_ang_velo_fs ~ wb_parameters$species))

###### absolute wing speed
ggplot(wb_parameters, aes(x = genus_species , y = abs_wing_speed, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'absolute wing speed [deg/sec]', x='taxa') + ylim(0,150000) + labs(fill = '') +
  ggtitle('variation in absolute wing speed ') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))
anova(lm(wb_parameters$abs_wing_speed ~ wb_parameters$species))

###### peak deviation rate (deviation velocity)
ggplot(wb_parameters, aes(x = genus_species , y = peak_dev_velo, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'peak deviation rate [deg/sec]', x='taxa') + ylim(0,70000) + labs(fill = '') +
  ggtitle('variation in peak deviation rate') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))
anova(lm(wb_parameters$peak_dev_velo ~ wb_parameters$species))

###### peak rotation rate (rotation velocity)
ggplot(wb_parameters, aes(x = genus_species , y = peak_rot_velo, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'peak rotation rate [deg/sec]', x='taxa') + ylim(0,280000) + labs(fill = '') +
  ggtitle('variation in peak rotation rate') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))
anova(lm(wb_parameters$peak_rot_velo ~ wb_parameters$species))

###### peak stroke acceleration
ggplot(wb_parameters, aes(x = genus_species , y = peak_ang_accel, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'peak stroke acceleration [deg/sec^-2]', x='taxa') + ylim(0,2e+08) + labs(fill = '') +
  ggtitle('variation in peak stroke acceleration ') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))

###### peak deviation acceleration
ggplot(wb_parameters, aes(x = genus_species , y = peak_dev_accel, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'peak deviation acceleration [deg/sec^-2]', x='taxa') + ylim(0,2e+08) + labs(fill = '') +
  ggtitle('variation in peak deviation acceleration ') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))

###### peak rotation acceleration
ggplot(wb_parameters, aes(x = genus_species , y = peak_rot_accel, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'peak rotation acceleration [deg/sec^-2]', x='taxa') + ylim(0,1.5e+09) + labs(fill = '') +
  ggtitle('variation in peak deviation acceleration ') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))

###### AoA at peak wing velocity
ggplot(wb_parameters, aes(x = genus_species , y = AoA_peak_ang_velo_fs, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'AoA at peak wing velocity [deg]', x='taxa') + ylim(0,100) + labs(fill = '') +
  ggtitle('variation in AoA at peak wing velocity') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))
anova(lm(wb_parameters$AoA_peak_ang_velo_fs ~ wb_parameters$species))


###### body speed
ggplot(wb_parameters, aes(x = genus_species , y = mean_body_speed, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'body speed [m.s-1]', x=' ') + ylim(0,1) + labs(fill = '') +
  ggtitle('variation in body speed') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))

###### climb angle
ggplot(wb_parameters, aes(x = genus_species , y = mean_climb_angle, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'climb angle [deg]', x='  ') + ylim(-100,100) + labs(fill = '') +
  ggtitle('variation in climb angle') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))


###### body pitch
ggplot(wb_parameters, aes(x = genus_species , y = mean_body_pitch, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'mean_body_pitch[deg]', x='taxa') + ylim(-70,70) + labs(fill = '') +
  ggtitle('variation in mean_body_pitch') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))

###### stroke plane inclination
ggplot(wb_parameters, aes(x = genus_species , y = stroke_plane_incl, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'stroke_plane_incl [deg]', x='taxa') + ylim(-70,70) + labs(fill = '') +
  ggtitle('variation in stroke_plane_incl') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))


###### stroke plane inclination relative to velocity vector
ggplot(wb_parameters, aes(x = genus_species , y = stroke_plane_incl_to_velo, fill = genus_species)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=discrete_weight_color) +
  # scale_fill_manual(values=c(rep('gray80',nb_species))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'stroke_plane_incl_to_velo [deg]', x='taxa') + ylim(-100,100) + labs(fill = '') +
  ggtitle('variation in stroke_plane_incl_to_velo') +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))



################ _________ create species-level dataset ################


## compute mean parameter value per species
mean_value_per_species = as.data.frame(matrix(NA, length(levels(wb_parameters$genus_species)), length(colnames(wb_parameters)[wb_parameters_num_cols])+1))
colnames(mean_value_per_species) = c('genus_species', colnames(wb_parameters)[wb_parameters_num_cols])
mean_value_per_species$genus_species = as.factor(levels(wb_parameters$genus_species))
for (i in 1:length(colnames(wb_parameters)[wb_parameters_num_cols])){
  mean_value_per_species[,(1+i)] = tapply(wb_parameters[,which(colnames(wb_parameters)==colnames(mean_value_per_species)[1+i])], wb_parameters$genus_species, mean)
}  # 1+i is because we inserted a first column for genus_species

## save 'mean_value_per_species' dataset for flight ~ morphology analyses:
# write.csv(mean_value_per_species, paste0(local_path, '/hoverflies_species_averaged_flight_dataset.csv'), row.names = F)

## compute the interspecific variance for each parameter
interspecific_variance <- as.data.frame(matrix(NA, length(colnames(wb_parameters)[wb_parameters_num_cols]), 2))
colnames(interspecific_variance) = c('parameter', 'interspecific_variation')
interspecific_variance$parameter = colnames(wb_parameters)[wb_parameters_num_cols]
for (i in 1:nrow(interspecific_variance)){
  interspecific_variance[i,2] = var(mean_value_per_species[,which(colnames(mean_value_per_species)==interspecific_variance$parameter[i])])
}

# ggplot(interspecific_variance, aes(x = parameter , y = interspecific_variation, fill = parameter)) +
#   geom_col() +
#   scale_fill_manual(values=c(rep('gray50',nrow(interspecific_variance)))) +
#   labs(y = 'between-species variance', x='parameters') + labs(fill = '') +
#   ggtitle('inter-specific variation among parameters') + coord_cartesian(ylim = c(0, 2e+16)) +
#   theme_classic() + theme(legend.position = "none")









# ################ ____ intra-specific variation ################

## compute the mean variance within species for each parameter
intraspecific_variance <- as.data.frame(matrix(NA, length(levels(wb_parameters$genus_species)) * length(colnames(wb_parameters)[wb_parameters_num_cols]), 3))
colnames(intraspecific_variance) = c('parameter', 'species', 'intraspecific_variation')
intraspecific_variance$parameter = rep(colnames(wb_parameters)[wb_parameters_num_cols], length(levels(wb_parameters$genus_species)))
## fill the species columns:
nb_of_param_now = length(colnames(wb_parameters)[wb_parameters_num_cols])
seq_now = seq(1, nrow(intraspecific_variance), by = nb_of_param_now)
for (i in 1:length(seq_now) ){
  intraspecific_variance$species[seq_now[i]:((seq_now[i]+ nb_of_param_now)-1)] = levels(wb_parameters$genus_species)[i]
}
## fill the variance per parameter for each species:
for (i in 1:nrow(intraspecific_variance)){
  intraspecific_variance$intraspecific_variation[i] = format( var(wb_parameters[which(wb_parameters$genus_species == intraspecific_variance$species[i]),
                                                                                which(colnames(wb_parameters)== intraspecific_variance$parameter[i])]), scientific = F)
}
intraspecific_variance$intraspecific_variation = as.numeric(intraspecific_variance$intraspecific_variation)
intraspecific_variance$species = as.factor(intraspecific_variance$species)
intraspecific_variance$parameter = as.factor(intraspecific_variance$parameter)


ggplot(intraspecific_variance, aes(x = parameter , y = log(intraspecific_variation), fill = parameter)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=c(rep('gray80',length(levels(intraspecific_variance$parameter))))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'intra-specific variance', x='parameters') + ylim(0,50) + labs(fill = '') +
  ggtitle('log(intra-specific variation) among parameters') +
  theme_classic() + theme(legend.position = "none")
print("variance in peak angular velocity trigger the warnings")

## interpret the above plot this way:
## each point is a species. Low variance values means that the species showed consistent (little variation in
## wingbeat kinematics over the 3 replicates (which likely include different individuals)).



ggplot(intraspecific_variance, aes(x = parameter , y = log(intraspecific_variation), fill = parameter)) +
  geom_boxplot(position=position_dodge(0.78), outlier.shape = NA) +
  scale_fill_manual(values=c(rep('gray80',length(levels(intraspecific_variance$parameter))))) +
  geom_jitter(position = position_jitter(0), cex = 2, color ='gray26') +
  labs(y = 'intra-specific variance', x='parameters') + ylim(0,50) + labs(fill = '') +
  ggtitle('log(intra-specific variation) among parameters') +
  theme_classic() + theme(legend.position = "none")








############################# [6] PCA on wingbeat parameters ############################# 

## many variable are highly correlated among each other
data_for_PCA = wb_parameters[,c(
  # which(colnames(wb_parameters) == 'mean_body_speed'),
  # which(colnames(wb_parameters) == 'mean_body_accel'),
  # which(colnames(wb_parameters) == 'mean_climb_angle'),
  which(colnames(wb_parameters) == 'WB_freq'),
  which(colnames(wb_parameters) == 'stroke_amplitude'),
  which(colnames(wb_parameters) == 'rotation_amplitude'),
  which(colnames(wb_parameters) == 'deviation_amplitude'),
  which(colnames(wb_parameters) == 'peak_ang_velo_fs'),
  which(colnames(wb_parameters) == 'peak_ang_velo_bs'),
  which(colnames(wb_parameters) == 'AoA_peak_ang_velo_fs'),
  which(colnames(wb_parameters) == 'AoA_peak_ang_velo_bs'),
  which(colnames(wb_parameters) == 'peak_dev_velo'),
  which(colnames(wb_parameters) == 'peak_rot_velo'),
  which(colnames(wb_parameters) == 'peak_ang_accel'),
  which(colnames(wb_parameters) == 'peak_rot_accel'),
  which(colnames(wb_parameters) == 'peak_dev_accel'),
  which(colnames(wb_parameters) == 'abs_wing_speed')
)]

## check for correlation
chart.Correlation(data_for_PCA, method='pearson', histogram=T, pch=16)

pca = princomp(data_for_PCA, scores=T, cor=T)
eig.val = get_eigenvalue(pca) ; eig.val$variance.percent
# fviz_pca_biplot(pca, axes = c(1, 2), repel = TRUE, label = 'var')
fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T)     # Avoid text overlapping
pca$loadings

plot(pca$scores[,2] ~ pca$scores[,1], pch=16, lwd=3, col=discrete_weight_color, cex=1.6, bty='n', asp=T, axes=F)
axis(1, at=seq(-5 ,5, by=5)) ; axis(2, at=seq(-2 ,4, by=2))
# text(pca$scores[,2] ~ pca$scores[,1], labels=wb_parameters$infraorder , pos=2, cex=0.8)

# plot legend only:
plot(0, col='white',axes=F, xlab='', ylab='');legend('center', pch=16, pt.cex=2.2, c(levels(wb_parameters$genus_species)), col=discrete_weight_color, bty='n', cex=0.9)








############################# compare WB freq and WB freq over the full sequence  ############################# 

## for each sequence from which we extracted a wingbeat, we counted the total number of wingbeats
## and computed the wingbeat frequency. We can thus assess if the frequency measured on 1 wingbeat
## reflects the frequency calculated over a larger number of wingbeats.


plot(IDs_file_seq_lvl$WB_freq, IDs_file_seq_lvl$full_seq_WB_freq, pch=16, lwd=3, col=colo_weight, cex=1.6, bty='n', xlab='WB frequency one wingbeat [Hz]', ylab='WB frequency full flight sequence [Hz]')
cor.test(IDs_file_seq_lvl$WB_freq, IDs_file_seq_lvl$full_seq_WB_freq)

mean(IDs_file_seq_lvl$WB_freq)
sd(IDs_file_seq_lvl$WB_freq)
mean(IDs_file_seq_lvl$full_seq_WB_freq)
sd(IDs_file_seq_lvl$full_seq_WB_freq)

mean(IDs_file_seq_lvl$WB_nb)
se <- function(x) sd(x)/sqrt(length(x))
se(IDs_file_seq_lvl$WB_nb)
sd(IDs_file_seq_lvl$WB_nb)






