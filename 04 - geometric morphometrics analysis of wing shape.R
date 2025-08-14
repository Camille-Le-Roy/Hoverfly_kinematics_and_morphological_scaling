

rm(list=ls())


library(MASS)
library(geomorph)
library(scico)
library(tagcloud)
library(phytools)

local_path = dirname(rstudioapi::getSourceEditorContext()$path)




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



########################### load morphometric data ########################### 

## reading tps files
wing_mat <- readland.tps(paste0(local_path,'/tps_files_digitized_wings/all_files.tps'))
dimnames(wing_mat)[[1]] = seq(1,301)
dimnames(wing_mat)[[2]] = c('X','Y')

## reading slider files 
sliders <- read.table(paste0(local_path,'/sliders_301ldm.txt'), h=T)

## ID files, including morphological measurements
morpho_dataset_ind_level <- read.csv(paste0(substr(local_path, 1, as.numeric(gregexpr(pattern ='/03_hoverfly_wings',local_path)[[1]])), 'hoverflies_weight_and_wing_data_ind_level.csv'))
morpho_dataset_ind_level$genus_species = as.factor(morpho_dataset_ind_level$genus_species)
filmed <- which(morpho_dataset_ind_level$origin=='field')

## make a full ID variable
full_ID <- rep(NA,dim(wing_mat)[[3]] )
for (i in 1:nrow(morpho_dataset_ind_level)){full_ID[i] = paste0(morpho_dataset_ind_level$session_ID[i],'_',morpho_dataset_ind_level$individual_ID[i])}
dimnames(wing_mat)[[3]] = full_ID


## plotting a wing shape
p=13
plot(wing_mat[,,p], pch=16, type = 'l', lwd = 3, col='black', cex=0.7, axes=F, xlab='',ylab='', asp=T, main = '')
# plot(wing_mat[,,p], pch=16, type = 'l', lwd = 3, col='black', cex=0.7, axes=F, xlab='',ylab='', asp=T, main = paste0(morpho_dataset_ind_level$genus_species[i],' ',full_ID[i]))



####### set colors for species

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

colo_species <- colors_species[match(morpho_dataset_ind_level$genus_species,levels(morpho_dataset_ind_level$genus_species))]

# plot legend only:
plot(0, col='white',axes=F, xlab='', ylab='');legend('center', pch=16, pt.cex=1.5, c(levels(morpho_dataset_ind_level$genus_species)), col=colors_species, bty='n', cex=0.6)


### create a color variable where the filmed specimens appear in white just to facilitate there superimposition
colo_species_filmed_invisible = colo_species
colo_species_filmed_invisible[which(colo_species != "grey75")] <- "white"


####### set color according to weight:
### set discrete color for species weight:
weight_color <- scico(25, palette = 'bamako')[1:22] #choose gradient among: scico_palette_show()
weight_color = weight_color[length(weight_color):1] # reverse the gradient
plot(rep(1,length(weight_color)), pch=15, cex=15,col=weight_color, bty='n', axes=F)
### set a gradient
variable_col_gradient = log(morpho_dataset_ind_level$fresh_weight)
gradient_now <- smoothPalette(variable_col_gradient, weight_color)
t <- colorRampPalette(c(weight_color[1], weight_color[length(weight_color)]))
seq_now = seq(min(variable_col_gradient), max(variable_col_gradient))



########################### ___ General Procruste Superimposition [individual-level] ########################### 
gpa <- gpagen(wing_mat, PrinAxes = F, ProcD = F, curves = sliders)
# summary(gpa)
plot(gpa)
plot(gpa$coords[,,1], pch=16, bty='n', cex=0.5, xlim=c(-0.08,0.08), ylim=c(-0.08,0.08), axes=F, asp=T)
for (i in seq(1,nrow(morpho_dataset_ind_level), by=10)){points(gpa$coords[,,i], pch=16, bty='n', cex=0.5, asp=T)}

## perform PCA on the procruste coordinates:
gpa_PCA <- gm.prcomp(gpa$coords)

## plot the result of the PCA
gpa_plot <- plot(gpa_PCA, axis1=1, axis2=2, pch = 16, bty='n', cex=2.1, col=colo_species_filmed_invisible)
points(gpa_plot$PC.points[,1][filmed], gpa_plot$PC.points[,2][filmed], pch=21, cex=2.1, col='black', bg=colo_species[filmed])


## set the body mass as a function of point size (log transformed because too large range of body mass)
gpa_plot <- plot(gpa_PCA, axis1=1, axis2=2, pch = 21, bty='n', cex=log10(morpho_dataset_ind_level$fresh_weight)*1.5, col='black', bg=colo_species)

## set the body mass as a color gradient
gpa_plot <- plot(gpa_PCA, axis1=1, axis2=2, pch = 16, bty='n', cex=2.5, col=gradient_now)
# text(gpa_plot$PC.points[,1], gpa_plot$PC.points[,2], labels = morpho_dataset_ind_level$genus_species, pos=2, cex=0.8)
# text(gpa_plot$PC.points[,1], gpa_plot$PC.points[,2], labels = full_ID, pos=2, cex=0.8)

## interactive function to visualize deformation grid
# picknplot.shape(gpa_plot)
# check https://cran.r-project.org/web/packages/geomorph/vignettes/geomorph.PCA.html


## manually produce deformation grids using 'plotRefToTarget'
## get the mean wing shape among all specimens
shape_consensus = mshape(gpa$coords)

## on PC 1
## Min values vs consensus
plotRefToTarget(shape_consensus, gpa_PCA$shapes$shapes.comp1$min, method = 'TPS')
## Max values vs consensus
plotRefToTarget(shape_consensus, gpa_PCA$shapes$shapes.comp1$max)

## on PC 2
## Min values vs consensus
plotRefToTarget(shape_consensus, gpa_PCA$shapes$shapes.comp2$min)
## Max values vs consensus
plotRefToTarget(shape_consensus, gpa_PCA$shapes$shapes.comp2$max)
# try method = c("TPS", "vector", "points", "surface") in 'plotRefToTarget'



## look at variation on the other PCs
gpa_plot <- plot(gpa_PCA, axis1=3, axis2=4, pch = 16, bty='n')
text(gpa_plot$PC.points[,1], gpa_plot$PC.points[,2], labels = morpho_dataset_ind_level$genus_species, pos=2, cex=0.8)
plotRefToTarget(shape_consensus, gpa_PCA$shapes$shapes.comp3$min)
plotRefToTarget(shape_consensus, gpa_PCA$shapes$shapes.comp3$max)

plotRefToTarget(shape_consensus, gpa_PCA$shapes$shapes.comp4$min)
plotRefToTarget(shape_consensus, gpa_PCA$shapes$shapes.comp4$max)





## correlation between PC 1 and body mass:
plot(gpa_plot$PC.points[,1] ~ log(morpho_dataset_ind_level$fresh_weight), pch = 21, bty='n', cex=2.3, col=colo_species_filmed_invisible, bg=colo_species_filmed_invisible, xlab='log(body mass)', ylab='PC 1: 69.21%')
points(gpa_plot$PC.points[,1][filmed] ~ log(morpho_dataset_ind_level$fresh_weight)[filmed], pch=21, cex=2.3, col='black', bg=colo_species[filmed])

plot(gpa_plot$PC.points[,2] ~ log(morpho_dataset_ind_level$fresh_weight), pch = 21, bty='n', cex=2.3, col=colo_species_filmed_invisible, bg=colo_species_filmed_invisible, xlab='log(body mass)', ylab='PC 2: 14.64%')
points(gpa_plot$PC.points[,2][filmed] ~ log(morpho_dataset_ind_level$fresh_weight)[filmed], pch=21, cex=2.3, col='black', bg=colo_species[filmed])
## cex=log10(morpho_dataset_ind_level$fresh_weight)*2 # to get point size as function of body mass

## correlation between PC 1 and non-dim S2:
plot(gpa_plot$PC.points[,1] ~ log(morpho_dataset_ind_level$non_dimensional_2ndMoment), pch = 21, bty='n', cex=2.3, col=colo_species_filmed_invisible, bg=colo_species_filmed_invisible, xlab='non-dim S2', ylab='PC 1: 69.21%')
points(gpa_plot$PC.points[,1][filmed] ~ log(morpho_dataset_ind_level$non_dimensional_2ndMoment)[filmed], pch=21, cex=2.3, col='black', bg=colo_species[filmed])

# plot(gpa_plot$PC.points[,1] ~ log(morpho_dataset_ind_level$non_dimensional_2ndMoment), pch = 21, bty='n', cex=log10(morpho_dataset_ind_level$fresh_weight)*2, col='black', bg=colo_species, xlab='non-dim S2', ylab='PC 1: 69.21%')
# cor.test(gpa_plot$PC.points[,1], log(morpho_dataset_ind_level$non_dimensional_2ndMoment))




########### ___ create species-level dataset (consensus shape per species) ########### 

### create a wing_ID dataset matching the species level: (order here is the one from levels())

wing_ID_species_lvl = as.data.frame(matrix(NA, length(levels(morpho_dataset_ind_level$genus_species)),1))
colnames(wing_ID_species_lvl) = c('genus_species')
wing_ID_species_lvl$genus_species = levels(morpho_dataset_ind_level$genus_species)

### compute the mean shape among the set of individual from each species
gpa_coords_species_lvl = array(NA, dim = c(301, 2, length(levels(morpho_dataset_ind_level$genus_species))))
dimnames(gpa_coords_species_lvl)[[2]] <-  dimnames(gpa$coords)[[2]]
dimnames(gpa_coords_species_lvl)[[3]] <-  levels(morpho_dataset_ind_level$genus_species)
for (i in 1:dim(gpa_coords_species_lvl)[3]){
  gpa_coords_species_lvl[,,i] = mshape(gpa$coords[,,which(morpho_dataset_ind_level$genus_species == levels(morpho_dataset_ind_level$genus_species)[i])])
}

### save the mean wing shape per species to send to Thomas Engels:
# for (i in 1:length(levels(morpho_dataset_ind_level$genus_species))) {
#   write.csv(gpa_coords_species_lvl[,,i], paste0(local_path,'/mean wing shape per species/',dimnames(gpa_coords_species_lvl)[[3]][i],'.csv'), row.names = F)
# }

## plotting a species mean wing shape
p = which(dimnames(gpa_coords_species_lvl)[[3]] == 'Eupeodes_nielseni')
plot(gpa_coords_species_lvl[,,p], pch=16, type = 'l', lwd = 3, col='black', cex=0.7, axes=F, xlab='',ylab='', asp=T, main = dimnames(gpa_coords_species_lvl)[[3]][p])
# for (p in length(dimnames(gpa_coords_species_lvl)[[3]]) : 1 ) {
#   plot(gpa_coords_species_lvl[,,p], pch=16, type = 'l', lwd = 3, col='black', cex=0.7, axes=F, xlab='',ylab='', asp=T, main = dimnames(gpa_coords_species_lvl)[[3]][p])
# }

### to have a species mean PCA, doing a gpa at the species level does not work well
### --> we compute the average values on PC 1 and 2 from the PCA on individual procruste coordinates

species_mean_PCs <- as.data.frame(matrix(NA, length(levels(morpho_dataset_ind_level$genus_species)), 3))
colnames(species_mean_PCs) = c('genus_species', 'PC1_mean_value', 'PC2_mean_value')
species_mean_PCs$genus_species = levels(morpho_dataset_ind_level$genus_species)
for (i in 1:nrow(species_mean_PCs)) {
  species_mean_PCs$PC1_mean_value[i] = mean(gpa_plot$PC.points[which(morpho_dataset_ind_level$genus_species==species_mean_PCs$genus_species[i]),1])
  species_mean_PCs$PC2_mean_value[i] = mean(gpa_plot$PC.points[which(morpho_dataset_ind_level$genus_species==species_mean_PCs$genus_species[i]),2])
}


######## morphospace (PC1 vs PC2) at the species level (Figure S3 in current ms)
### first check at the individual level
gpa_plot <- plot(gpa_PCA, axis1=1, axis2=2, pch = 16, bty='n', cex=2.1, col=colo_species_filmed_invisible)
points(gpa_plot$PC.points[,1][filmed], gpa_plot$PC.points[,2][filmed], pch=21, cex=2.1, col='black', bg=colo_species[filmed])
## then at the species level
plot(species_mean_PCs$PC1_mean_value, species_mean_PCs$PC2_mean_value, pch = 21, bty='n', cex=3, col='black', bg=colors_species, xlab='PC 1 (69.21% var. expl.)', ylab='PC 2 (14.64% var. expl.)')


############ phylogenetic PCA on wing shape ############ 

### load the hoverfly phylogeny (from Wong et al. 2023. The phylogeny and evolutionary ecology of hoverflies)
### the tree was pruned and formatted in the code "format_hoverfly_phylogeny.R"
### tree matching the large morphological dataset (8 filmed species + 20 Leiden museum species:
load(paste0(substr(local_path, 1, as.numeric(gregexpr(pattern ='/03_hoverfly_wings',local_path)[[1]])), 'phylogeny from Daniel Wong et al. 2023/hoverfly_tree_(morphology_dataset).RData'))
plot(phy, cex=0.6, edge.width = 0.5, use.edge.length = F)

## re-order 'species_mean_PCs' to match phylogeny
species_mean_PCs_tip_order <- species_mean_PCs[match(phy$tip.label, species_mean_PCs$genus_species), ]
rownames(species_mean_PCs_tip_order) <- phy$tip.label

## order the colors to follow tip order as well:
colors_named <- setNames(colors_species, species_mean_PCs$genus_species)
colors_species_tip_order <- colors_named[phy$tip.label]

## perform the phylogenetic PCA
phylo_pca = phyl.pca(phy, species_mean_PCs_tip_order[,2:3], method="BM")

# plot the first two PCs
plot(-phylo_pca$S[,1], -phylo_pca$S[,2], pch = 21, bty='n', cex=3, col='black', bg=colors_species_tip_order, main = "Phylogenetic PCA")
# percentage of variance explained by each PC
eigenvalues <- phylo_pca$Eval ; total_variance <- sum(eigenvalues)
variance_explained <- (eigenvalues / total_variance) * 100

#### the results of the phylogenetic PCA is highly similar to the traditional PCA
#### plots are basically the same





##################### PC vs body mass and S2 ##################### 

## load the species level datasets built in the code "BBhov_analyzing_WBkinematics_vs_morphology.R"
# full_dataset <- read.csv(paste0(substr(local_path, 1, as.numeric(gregexpr(pattern ='/03_hoverfly_wings',local_path)[[1]])),'full_dataset_species_level.csv'))
morpho_dataset_species_level <- read.csv(paste0(substr(local_path, 1, as.numeric(gregexpr(pattern ='/03_hoverfly_wings',local_path)[[1]])),'hoverflies_weight_and_wing_data_species_level.csv'))



##################### ____ species level ##################### 

## PC1 vs body mass
plot(log(morpho_dataset_species_level$fresh_weight), species_mean_PCs$PC1_mean_value, pch = 21, bty='n', cex=3, col='black', bg=colors_species, ylab='PC 1', xlab='log(body mass)')
## PC2 vs body mass
plot(log(morpho_dataset_species_level$fresh_weight), species_mean_PCs$PC2_mean_value, pch = 21, bty='n', cex=3, col='black', bg=colors_species, ylab='PC 2', xlab='log(body mass)')

## test the relationships
cor.test(morpho_dataset_species_level$fresh_weight, species_mean_PCs$PC1_mean_value) 
cor.test(morpho_dataset_species_level$non_dimensional_2ndMoment, species_mean_PCs$PC1_mean_value) 

## correlation with body mass and dimensionless S2 (Figure "wing shape body and S2" in ms)
plot(species_mean_PCs$PC1_mean_value, morpho_dataset_species_level$non_dimensional_2ndMoment, pch = 21, bty='n', cex=3, col='black', bg=colors_species, xlab='PC 1 (69.21% var. expl.)', ylab='non-dim S2')
plot(species_mean_PCs$PC1_mean_value, log(morpho_dataset_species_level$fresh_weight), pch = 21, bty='n', cex=3, col='black', bg=colors_species, xlab='PC 1 (69.21% var. expl.)', ylab='log(body mass)')

plot(log(morpho_dataset_species_level$non_dimensional_2ndMoment), log(morpho_dataset_species_level$fresh_weight), pch = 21, bty='n', cex=3, col='black', bg=colors_species, xlab='log(non-dim S2)', ylab='log(body mass)')

### (results looks very similar as if we do a gpa on species mean shape).


#### plot in log scale but with linear labels
plot(abs(species_mean_PCs$PC1_mean_value), morpho_dataset_species_level$fresh_weight, pch = 21, bty='n', cex=3, col='black', bg=colors_species, xlab='PC 1 (69.21% var. expl.)', ylab='body mass[mg]', log='xy', axes=F)
axis(1, at=log10Tck('x','minor'), tcl= 0.2)
axis(2, at=log10Tck('y','minor'), tcl= 0.2)




##################### ____ individual level ##################### 

### correlation at the individual level
plot(gpa_plot$PC.points[,1], log(morpho_dataset_ind_level$fresh_weight), pch = 16, bty='n', cex=2.5, col=colo_species_filmed_invisible, xlab='PC 1: 69.21%', ylab='log(body mass)')
points(gpa_plot$PC.points[,1][filmed], log(morpho_dataset_ind_level$fresh_weight)[filmed], pch = 21, col='black', bg=colo_species, cex=2.5)

plot(log(morpho_dataset_ind_level$non_dimensional_2ndMoment), log(morpho_dataset_ind_level$fresh_weight), pch = 16, bty='n', cex=2.5, col=colo_species_filmed_invisible, xlab='log(non-dim S2)', ylab='log(body mass)')
points(log(morpho_dataset_ind_level$non_dimensional_2ndMoment)[filmed], log(morpho_dataset_ind_level$fresh_weight)[filmed], pch = 21, col='black', bg=colo_species, cex=2.5)


### for figure, we plot the species level and the individual level plot of body mass vs. PC1:
plot(gpa_plot$PC.points[,1], log(morpho_dataset_ind_level$fresh_weight), pch = 21, bty='n', cex=2, col='black', bg=colo_species, xlab='PC 1: 69.21%', ylab='log(body mass)')
points(species_mean_PCs$PC1_mean_value, log(morpho_dataset_species_level$fresh_weight), pch = 21, cex=3, col='black', bg=colors_species)

plot(log(morpho_dataset_ind_level$non_dimensional_2ndMoment), log(morpho_dataset_ind_level$fresh_weight), pch = 21, bty='n', cex=2, col='black', bg=colo_species, xlab='log(non-dim S2)', ylab='log(body mass)')
points(log(morpho_dataset_species_level$non_dimensional_2ndMoment), log(morpho_dataset_species_level$fresh_weight), pch = 21, cex=3, col='black', bg=colors_species)

plot(gpa_plot$PC.points[,1], morpho_dataset_ind_level$non_dimensional_2ndMoment, pch = 21, bty='n', cex=2, col='black', bg=colo_species, xlab='PC 1: 69.21%', ylab='S2*')
points(species_mean_PCs$PC1_mean_value, morpho_dataset_species_level$non_dimensional_2ndMoment, pch = 21, cex=3, col='black', bg=colors_species)








