#This is an R script that was used to analyse the data for a study Wutkowska et al., 'Dead or alive...'  

R.version.string #[1] "R version 3.4.4 (2018-03-15)"

###################################################################################
####Read OTU table, taxonomy/functional assigment table and edaphic parameters ####
###################################################################################

otus<-read.delim("/Users/magdalenawutkowska/Dropbox/PhD/PROJECT_Snowfence/3_stats/publication_reanalysisfromscratch/otu_table_rarified_sorted.txt",sep="\t",dec=".",header=TRUE,row.names=1)
otus_t <- t(otus)
otus_t <- data.frame(otus_t)
otus_t$rowsums <- rowSums(otus_t)
otus_t$rowsums

tax_fun<-read.delim("/Users/magdalenawutkowska/Dropbox/PhD/PROJECT_Snowfence/3_stats/publication_reanalysisfromscratch/tax&fun_assignments.txt", sep="\t", dec= ".", header=TRUE, row.names=1)
tax_fun<-tax_fun[1:837,] #removing empty row
as.factor(tax_fun$Simplified_Trophic_mode) 
levels(tax_fun$Simplified_Trophic_mode) #checking levels of trophic mode

sample_names <- row.names(otus_t)
sample_names

env<-read.delim("/Users/magdalenawutkowska/Dropbox/PhD/PROJECT_Snowfence/3_stats/publication_reanalysisfromscratch/env.csv", sep=",",dec=".",header=TRUE,row.names=1)
env_dna <- env[1:23,]
env_rna <- env[24:42,]
all(rownames(otus)==rownames(tax_fun)) #check order of OTUs in both files


#######################################################################################################################################
####GNMDS was based on script https://datadryad.org/bitstream/handle/10255/dryad.74203/Script%20for%20NMDS%20analysis.R?sequence=1 ####
#######################################################################################################################################
####FOR ALL OTUS abundance

library(vegan) #This is vegan 2.5-1
library(MASS) #This is MASS 7.3-47
library(stats) #This is stats 3.3.2

#otu tables are already prepared in the working space
attach(otus)
names(otus)
str(otus)

#environmental data
attach(env)
names(env)
str(env)

#making Bray-Curtis dissimilarity matrix:
dist.all.abu<-vegdist(t(otus), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
dist.all.abu 

#define a general, empty object called mds:
mds.all.abu<-NULL
#making 100 "mds"s from initial starting configurations, allocating them into the mds object:
for(i in 1:100)
{mds.all.abu[[i]]<-isoMDS(dist.all.abu,initMDS(dist.all.abu, k=2), k=2, maxit=1000,tol=1e-7)}
#k determines the number of dimensions in the ordination
#the mds object is now a list consisting of 100 "subobjects" being lists
# extracting the stress values as a vector - a measure of the correspondence between dissimilarities
# and distance for each ordination:
mds.stress.all.abu<-unlist(lapply(mds.all.abu,function(v){v[[2]]})) 

dist.all.abu.clust<-hclust(dist.all.abu,"single")      # hierarchical clustering 
plot(dist.all.abu.clust,ylab="Dissimilarity",xlab="",sub="")  # plot dendrogram

# looking at the stress values for 100 mds:
mds.stress.all.abu
#ordering the stress values for the 100 mds:
order(mds.stress.all.abu)
#Saving the order in a vector
ordered.all.abu<-order(mds.stress.all.abu)
ordered.all.abu

#find the stress of the solutions with the lowest and second lowest stress:
mds.stress.all.abu[ordered.all.abu[1]]
mds.stress.all.abu[ordered.all.abu[2]]

#scaling of axes to half change units and varimax rotation
mds.best.all.abu<-postMDS(mds.all.abu[[ordered.all.abu[1]]],dist.all.abu)
mds.best.all.abu
mds.secbest.all.abu<-postMDS(mds.all.abu[[ordered.all.abu[2]]],dist.all.abu)
mds.secbest.all.abu

#Procrustes comparisons
procrustes(mds.best.all.abu,mds.secbest.all.abu,permutations=999)
protest(mds.best.all.abu,mds.secbest.all.abu,permutations=999)
plot(procrustes(mds.best.all.abu,mds.secbest.all.abu,permutations=999))

#making variables from GNMDS axes 1 and 2 for plotting
gnmds1_all.abu<-mds.best.all.abu$points[,1]
gnmds2_all.abu<-mds.best.all.abu$points[,2]

##Biplot ##both using origin=TRUE here
mds.df.all.abu<-data.frame(gnmds1_all.abu,gnmds2_all.abu)
fit.all.abu<-envfit(mds.df.all.abu, env[,c(1, 2, 6, 7, 8, 9, 10, 12, 13, 14)], 999)
fit.all.abu

# ***VECTORS
# 
# gnmds1_all.abu gnmds2_all.abu     r2 Pr(>r)    
# pH                 -0.97843        0.20659 0.7459  0.001 ***
# moisture           -0.28034       -0.95990 0.0813  0.191    
# conductivity       -0.78960        0.61363 0.1832  0.019 *  
# OM                 -0.98228       -0.18743 0.2718  0.002 ** 
# total_Nplus1       -0.81539       -0.57892 0.4361  0.001 ***
# C                  -0.87940       -0.47608 0.3474  0.001 ***
# cn_ratio            0.69538        0.71865 0.2247  0.011 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 999
# 
# ***FACTORS:
#   
#   Centroids:
#   gnmds1_all.abu gnmds2_all.abu
# vegetation_typeheath         -0.1837         0.0054
# vegetation_typemeadow         0.2020        -0.0059
# snow_regimeambient           -0.0060         0.0104
# snow_regimedeep               0.0054        -0.0094
# lumpedDNA_H                  -0.2143        -0.0143
# lumpedDNA_M                   0.2428        -0.0171
# lumpedRNA_H                  -0.1530         0.0250
# lumpedRNA_M                   0.1409         0.0110
# 
# Goodness of fit:
#   r2 Pr(>r)    
# vegetation_type 0.3779  0.001 ***
#   snow_regime     0.0013  0.940    
# lumped          0.3980  0.001 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 999


#Making the ordination diagram using the best solution:
par(mfrow=c(1,1))
plot(main="GNMDS_all.abu", gnmds1_all.abu,gnmds2_all.abu,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)
#points(gnmds1_all.abu[nucleic_acid=="DNA"], gnmds2_all.abu[nucleic_acid=="DNA"], pch=19, col='black')
#points(gnmds1_all.abu[nucleic_acid=="RNA"], gnmds2_all.abu[nucleic_acid=="RNA"], pch=21, col='black')
points(gnmds1_all.abu[treatment=="RNA_H_A"], gnmds2_all.abu[treatment=="RNA_H_A"], pch=19, col='black')
points(gnmds1_all.abu[treatment=="RNA_H_D"], gnmds2_all.abu[treatment=="RNA_H_D"], pch=21, col='black')
points(gnmds1_all.abu[treatment=="RNA_M_A"], gnmds2_all.abu[treatment=="RNA_M_A"], pch=19, col='grey')
points(gnmds1_all.abu[treatment=="RNA_M_D"], gnmds2_all.abu[treatment=="RNA_M_D"], pch=21, col='grey')
points(gnmds1_all.abu[treatment=="DNA_H_A"], gnmds2_all.abu[treatment=="DNA_H_A"], pch=15, col='black')
points(gnmds1_all.abu[treatment=="DNA_H_D"], gnmds2_all.abu[treatment=="DNA_H_D"], pch=22, col='black')
points(gnmds1_all.abu[treatment=="DNA_M_A"], gnmds2_all.abu[treatment=="DNA_M_A"], pch=15, col='grey')
points(gnmds1_all.abu[treatment=="DNA_M_D"], gnmds2_all.abu[treatment=="DNA_M_D"], pch=22, col='grey')
#ordiellipse(mds.df.all.abu, env$treatment, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="black", border = "black")
ordiellipse(mds.df.all.abu, env$vegetation_type, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "black")
#ordiellipse(mds.df.all.abu, env$snow_regime, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="black", border = "black")
ordiellipse(mds.df.all.abu, env$nucleic_acid, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="red", border = "black")
#legend("bottomright", legend = c("DNA", "RNA") , pch=c(19,21), col=c("black", "black"), bty="n")
legend("topleft", legend = c("RNA_H_A", "RNA_H_D", "RNA_M_A", "RNA_M_D", "DNA_H_A", "DNA_H_D", "DNA_M_A", "DNA_M_D") , pch=c(19,21,19,21,15,22,15,22), col=c("black", "black", "grey", "grey", "black", "black", "grey", "grey"), bty="n")
plot(fit.all.abu)

#Lumped treatments:
par(mfrow=c(1,1))
plot(main="GNMDS_all.abu", gnmds1_all.abu,gnmds2_all.abu,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)
points(gnmds1_all.abu[lumped=="RNA_H"], gnmds2_all.abu[lumped=="RNA_H"], pch=19, col='black') 
points(gnmds1_all.abu[lumped=="RNA_M"], gnmds2_all.abu[lumped=="RNA_M"], pch=21, col='black')
points(gnmds1_all.abu[lumped=="DNA_H"], gnmds2_all.abu[lumped=="DNA_H"], pch=19, col='grey') 
points(gnmds1_all.abu[lumped=="DNA_M"], gnmds2_all.abu[lumped=="DNA_M"], pch=21, col='grey')
#ordiellipse(mds.df.all.abu, env$lumped, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="blue", border = "black")
ordiellipse(mds.df.all.abu, env$vegetation_type, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "red")
#ordiellipse(mds.df.all.abu, env$snow_regime, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="black", border = "black")
ordiellipse(mds.df.all.abu, env$nucleic_acid, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="red", border = "blue")
legend("topleft", legend = c("RNA_H", "RNA_M", "DNA_H", "DNA_M") , pch=c(19,21,19,21), col=c("black", "black", "grey", "grey"), bty="n")
#plot(fit.all.abu)

detach(otus)

#######################################################################################################################################
####GNMDS was based on script https://datadryad.org/bitstream/handle/10255/dryad.74203/Script%20for%20NMDS%20analysis.R?sequence=1 ####
#######################################################################################################################################
####FOR ALL OTUS presence - absence

library(vegan) #This is vegan 2.5-1
library(MASS) #This is MASS 7.3-47
library(stats) #This is stats 3.3.2

#converting matrix to presence-absence
otus_pa <- otus
otus_pa[otus_pa>0] <-1

#otu tables are already prepared in the working space
attach(otus_pa)
names(otus_pa)
str(otus_pa)

#making Bray-Curtis dissimilarity matrix:
dist.all.pa<-vegdist(t(otus_pa), method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
dist.all.pa 

#define a general, empty object called mds:
mds.all.pa<-NULL
#making 100 "mds"s from initial starting configurations, allocating them into the mds object:
for(i in 1:100)
{mds.all.pa[[i]]<-isoMDS(dist.all.pa,initMDS(dist.all.pa, k=2), k=2, maxit=1000,tol=1e-7)}
#k determines the number of dimensions in the ordination
#the mds object is now a list consisting of 100 "subobjects" being lists
# extracting the stress values as a vector - a measure of the correspondence between dissimilarities
# and distance for each ordination:
mds.stress.all.pa<-unlist(lapply(mds.all.pa,function(v){v[[2]]})) 

dist.all.pa.clust<-hclust(dist.all.pa,"single")      # hierarchical clustering 
plot(dist.all.pa.clust,ylab="Dissimilarity",xlab="",sub="")  # plot dendrogram

# looking at the stress values for 100 mds:
mds.stress.all.pa
#ordering the stress values for the 100 mds:
order(mds.stress.all.pa)
#Saving the order in a vector
ordered.all.pa<-order(mds.stress.all.pa)
ordered.all.pa

#find the stress of the solutions with the lowest and second lowest stress:
mds.stress.all.pa[ordered.all.pa[1]]
mds.stress.all.pa[ordered.all.pa[2]]

#scaling of axes to half change units and varimax rotation
mds.best.all.pa<-postMDS(mds.all.pa[[ordered.all.pa[1]]],dist.all.pa)
mds.best.all.pa
mds.secbest.all.pa<-postMDS(mds.all.pa[[ordered.all.pa[2]]],dist.all.pa)
mds.secbest.all.pa

#Procrustes comparisons
procrustes(mds.best.all.pa,mds.secbest.all.pa,permutations=999)
protest(mds.best.all.pa,mds.secbest.all.pa,permutations=999)
plot(procrustes(mds.best.all.pa,mds.secbest.all.pa,permutations=999))

#making variables from GNMDS axes 1 and 2 for plotting
gnmds1_all.pa<-mds.best.all.pa$points[,1]
gnmds2_all.pa<-mds.best.all.pa$points[,2]

##Biplot ##both using origin=TRUE here
mds.df.all.pa<-data.frame(gnmds1_all.pa,gnmds2_all.pa)
fit.all.pa<-envfit(mds.df.all.pa, env[,c(1, 2, 6, 7, 8, 9, 10, 12, 13, 14)], 999)
fit.all.pa

# ***VECTORS
# 
# gnmds1_all.pa gnmds2_all.pa     r2 Pr(>r)    
# pH                -1.00000      -0.00269 0.7561  0.001 ***
# moisture           0.53810      -0.84288 0.0388  0.468    
# conductivity      -0.89133      -0.45336 0.1326  0.060 .  
# OM                -0.73242      -0.68085 0.2871  0.002 ** 
# total_Nplus1      -0.86674      -0.49876 0.3287  0.003 ** 
# C                 -0.76163      -0.64802 0.3184  0.003 ** 
# cn_ratio           0.65055      -0.75947 0.1772  0.026 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 999
# 
# ***FACTORS:
#   
#   Centroids:
#   gnmds1_all.pa gnmds2_all.pa
# vegetation_typeheath        -0.2318        0.0105
# vegetation_typemeadow        0.2550       -0.0115
# snow_regimeambient           0.0324       -0.0074
# snow_regimedeep             -0.0295        0.0067
# lumpedDNA_H                 -0.2696       -0.0649
# lumpedDNA_M                  0.2716       -0.1274
# lumpedRNA_H                 -0.1940        0.0858
# lumpedRNA_M                  0.2300        0.1624
# 
# Goodness of fit:
#   r2 Pr(>r)    
# vegetation_type 0.5645  0.001 ***
#   snow_regime     0.0096  0.686    
# lumped          0.6933  0.001 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 999

plot(fit.all.pa)

#Making the ordination diagram using the best solution:
par(mfrow=c(1,1))
plot(main="GNMDS_all.pa", gnmds1_all.pa,gnmds2_all.pa,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)
#points(gnmds1_all.abu[nucleic_acid=="DNA"], gnmds2_all.abu[nucleic_acid=="DNA"], pch=19, col='black')
#points(gnmds1_all.abu[nucleic_acid=="RNA"], gnmds2_all.abu[nucleic_acid=="RNA"], pch=21, col='black')
points(gnmds1_all.pa[treatment=="RNA_H_A"], gnmds2_all.pa[treatment=="RNA_H_A"], pch=19, col='black') 
points(gnmds1_all.pa[treatment=="RNA_H_D"], gnmds2_all.pa[treatment=="RNA_H_D"], pch=21, col='black')
points(gnmds1_all.pa[treatment=="RNA_M_A"], gnmds2_all.pa[treatment=="RNA_M_A"], pch=19, col='grey')
points(gnmds1_all.pa[treatment=="RNA_M_D"], gnmds2_all.pa[treatment=="RNA_M_D"], pch=21, col='grey')
points(gnmds1_all.pa[treatment=="DNA_H_A"], gnmds2_all.pa[treatment=="DNA_H_A"], pch=15, col='black') 
points(gnmds1_all.pa[treatment=="DNA_H_D"], gnmds2_all.pa[treatment=="DNA_H_D"], pch=22, col='black')
points(gnmds1_all.pa[treatment=="DNA_M_A"], gnmds2_all.pa[treatment=="DNA_M_A"], pch=15, col='grey')
points(gnmds1_all.pa[treatment=="DNA_M_D"], gnmds2_all.pa[treatment=="DNA_M_D"], pch=22, col='grey')
#ordiellipse(mds.df.all.abu, env$treatment, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = FALSE, lwd=1, col="black", border = "black")
ordiellipse(mds.df.all.abu, env$vegetation_type, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "black")
#ordiellipse(mds.df.all.abu, env$snow_regime, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="black", border = "black")
ordiellipse(mds.df.all.abu, env$nucleic_acid, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="red", border = "black")
#legend("bottomright", legend = c("DNA", "RNA") , pch=c(19,21), col=c("black", "black"), bty="n")
legend("topleft", legend = c("RNA_H_A", "RNA_H_D", "RNA_M_A", "RNA_M_D", "DNA_H_A", "DNA_H_D", "DNA_M_A", "DNA_M_D") , pch=c(19,21,19,21,15,22,15,22), col=c("black", "black", "grey", "grey", "black", "black", "grey", "grey"), bty="n")
plot(fit.all.pa)

#Lumped treatments:
par(mfrow=c(1,1))
plot(main="GNMDS_all.pa", gnmds1_all.pa,gnmds2_all.pa,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)
points(gnmds1_all.pa[lumped=="RNA_H"], gnmds2_all.pa[lumped=="RNA_H"], pch=19, col='black') 
points(gnmds1_all.pa[lumped=="RNA_M"], gnmds2_all.pa[lumped=="RNA_M"], pch=21, col='black')
points(gnmds1_all.pa[lumped=="DNA_H"], gnmds2_all.pa[lumped=="DNA_H"], pch=19, col='grey') 
points(gnmds1_all.pa[lumped=="DNA_M"], gnmds2_all.pa[lumped=="DNA_M"], pch=21, col='grey')
#ordiellipse(mds.df.all.pa, env$lumped, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="black", border = "black")
ordiellipse(mds.df.all.pa, env$vegetation_type, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "red")
#ordiellipse(mds.df.all.pa, env$snow_regime, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="black", border = "black")
ordiellipse(mds.df.all.pa, env$nucleic_acid, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="red", border = "blue")
legend("bottomright", legend = c("RNA_H", "RNA_M", "DNA_H", "DNA_M") , pch=c(19,21,19,21), col=c("black", "black", "grey", "grey"), bty="n")

#Lumped, but divided into dna/rna + veg
par(mfrow=c(1,1))
plot(main="GNMDS_all.pa", gnmds1_all.pa,gnmds2_all.pa,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)
points(gnmds1_all.pa[lumped=="RNA_H"], gnmds2_all.pa[lumped=="RNA_H"], pch=19, col='black') 
points(gnmds1_all.pa[lumped=="RNA_M"], gnmds2_all.pa[lumped=="RNA_M"], pch=21, col='black')
points(gnmds1_all.pa[lumped=="DNA_H"], gnmds2_all.pa[lumped=="DNA_H"], pch=19, col='grey') 
points(gnmds1_all.pa[lumped=="DNA_M"], gnmds2_all.pa[lumped=="DNA_M"], pch=21, col='grey')
ordiellipse(mds.df.all.pa, env$lumped, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="red", border = "red")
legend("bottomright", legend = c("RNA_H", "RNA_M", "DNA_H", "DNA_M") , pch=c(19,21,19,21), col=c("black", "black", "grey", "grey"), bty="n")


detach(otus_pa)
detach(env)


########################################
####RELATIVE ABUNDANCE OF ALL READS ####
########################################

library(ggplot2) #this is ggplot 2.2.1
library(tidyverse) #this is tidyverse 1.1.1
library(dplyr) #this is dplyr 0.7.4
library(reshape) #this is reshape 0.8.6
library(RColorBrewer) #this is RColorBrewer 1.1-2

otus_t <- t(otus)
all(colnames(otus_t)==rownames(tax_fun)) #checking if the numbers are matching
all(rownames(otus_t)==rownames(env))

env$richness<-specnumber(otus_t) #number of OTUs per sample

colnames(tax_fun)
tax_fun$curated.trophic.mode <- as.character(tax_fun$Simplified_Trophic_mode)
#tax_fun$curated.trophic.mode <- ifelse(is.na(tax_fun$agreed.trophic.mode), 'unassigned', tax_fun$curated.trophic.mode)
tax_fun$curated.trophic.mode

#Summarizing sequences by curated.trophic.mode
total.por.otu<-colSums(otus_t)
total<-data.frame(total.por.otu,tax_fun$curated.trophic.mode)

total.agg1<-aggregate(total[,1],list(total[,2]),sum) #sums all sequences belonging to a specific curated.trophic.level
colnames(total.agg1)<-c("Trophic mode","Totals")
total.agg1$proportion<-total.agg1[,2]*100/sum(total.agg1[,2])
pos<-order(total.agg1$proportion,decreasing=T)
totagg<-total.agg1[pos,]

rownames(otus)==rownames(tax_fun)    
otus$tax<-tax_fun$curated.trophic.mode
head(otus)  
otustax<-otus
dim(otustax)
otustax<-droplevels(otustax)

otustax$tax
test_grouped<-group_by(otustax,tax) #grouping by tax
str(test_grouped)
test<-summarise_all(test_grouped, sum)
str(test)
test<-data.frame(test)
test

sumDNA_heath <- rowSums(test[, c(2:12)])
sumDNA_meadow <- rowSums(test[, c(13:24)])
sumRNA_heath <- rowSums(test[, c(25:35)])
sumRNA_meadow <- rowSums(test[, c(36:43)])
abund <- c(sumDNA_heath, sumDNA_meadow, sumRNA_heath, sumRNA_meadow)
tax <- rep(c("Pathotroph", "Saprotroph", "Symbiotroph", "unassigned"), times=4)
nucleic_acid <- c(rep("DNA", 8), rep("RNA", 8))
vegetation <- c(rep("heath", 4), rep("meadow", 4), rep("heath", 4), rep("meadow", 4))
test_sum <- cbind(tax, abund, nucleic_acid, vegetation)

test_sum <- data.frame(test_sum)
test_sum$abund <- as.numeric(as.character(test_sum$abund))

###plotting proportion of all reads by trophic modes divided into rDNA/rRNA, snow regimes and separate vegetation types
ggplot(test_sum, aes(x=nucleic_acid, y=abund, fill=tax))+
  geom_bar(position="fill",stat="identity")+ 
  facet_wrap(~vegetation)+
  theme(panel.grid.minor.y=element_blank(),panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),axis.ticks=element_line(size=0.2),
        aspect.ratio=1,text=element_text(size=11),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11,angle=90,hjust=1),
        plot.title=element_text(size=10, face='bold',margin=margin(0,0,0,0)))+
  labs(title=NULL,x=NULL, y="Proportion of all reads")+
  theme_bw()+
  theme(legend.position="bottom")
  #labs(x = " ", y="Richness of functionally unassigned OTUs")?

#counting relative abundances
melted_combi <- split(test_sum, list(test_sum[["vegetation"]], test_sum[["nucleic_acid"]]))
melted_funct <- split(test_sum, list(test_sum[["vegetation"]], test_sum[["nucleic_acid"]], test_sum[["tax"]]))

(sum(melted_funct$heath.DNA.Symbiotroph$abund)/sum(melted_combi$heath.DNA$abund))*100 #91.67679
(sum(melted_funct$heath.RNA.Symbiotroph$abund)/sum(melted_combi$heath.RNA$abund))*100 #76.27565

(sum(melted_funct$meadow.DNA.Symbiotroph$abund)/sum(melted_combi$meadow.DNA$abund))*100 #87.30053
(sum(melted_funct$meadow.RNA.Symbiotroph$abund)/sum(melted_combi$meadow.RNA$abund))*100 #80.72132



(sum(melted_funct$heath.DNA.Saprotroph$abund)/sum(melted_combi$heath.DNA$abund))*100 #3.502379
(sum(melted_funct$heath.RNA.Saprotroph$abund)/sum(melted_combi$heath.RNA$abund))*100 #9.492734

(sum(melted_funct$meadow.DNA.Saprotroph$abund)/sum(melted_combi$meadow.DNA$abund))*100 #3.680255
(sum(melted_funct$meadow.RNA.Saprotroph$abund)/sum(melted_combi$meadow.RNA$abund))*100 #7.004625



(sum(melted_funct$heath.DNA.Pathotroph$abund)/sum(melted_combi$heath.DNA$abund))*100 #0.2978381
(sum(melted_funct$heath.RNA.Pathotroph$abund)/sum(melted_combi$heath.RNA$abund))*100 #0.1093357

(sum(melted_funct$meadow.DNA.Pathotroph$abund)/sum(melted_combi$meadow.DNA$abund))*100 #0.156907
(sum(melted_funct$meadow.RNA.Pathotroph$abund)/sum(melted_combi$meadow.RNA$abund))*100 #0.4724863



(sum(melted_funct$heath.DNA.unassigned$abund)/sum(melted_combi$heath.DNA$abund))*100 #4.522988
(sum(melted_funct$heath.RNA.unassigned$abund)/sum(melted_combi$heath.RNA$abund))*100 #14.12228

(sum(melted_funct$meadow.DNA.unassigned$abund)/sum(melted_combi$meadow.DNA$abund))*100 #8.862306
(sum(melted_funct$meadow.RNA.unassigned$abund)/sum(melted_combi$meadow.RNA$abund))*100 #11.80157







###################################################
###   GNMDS for DNA presence absence results    ###
###################################################

library(vegan) #This is vegan 2.5-1
library(MASS) #This is MASS 7.3-47
library(stats) #This is stats 3.3.2

#subsampling DNA samples from presence-absence dataset
otus_DNA <- data.frame(otus_pa[1:23])

#otu tables are already prepared in the working space
attach(otus_DNA)
names(otus_DNA)
str(otus_DNA)
attach(env_dna)

#making Bray-Curtis dissimilarity matrix:
dist.DNA<-vegdist(t(otus_DNA), method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
dist.DNA 

#define a general, empty object called mds:
mds.DNA<-NULL
#making 100 "mds"s from initial starting configurations, allocating them into the mds object:
for(i in 1:100)
{mds.DNA[[i]]<-isoMDS(dist.DNA,initMDS(dist.DNA, k=2), k=2, maxit=1000,tol=1e-7)}
#k determines the number of dimensions in the ordination
#the mds object is now a list consisting of 100 "subobjects" being lists
# extracting the stress values as a vector - a measure of the correspondence between dissimilarities
# and distance for each ordination:
mds.stress.DNA<-unlist(lapply(mds.DNA,function(v){v[[2]]})) 

dist.DNA.clust<-hclust(dist.DNA,"single")      # hierarchical clustering 
plot(dist.DNA.clust,ylab="Dissimilarity",xlab="",sub="")  # plot dendrogram

# looking at the stress values for 100 mds:
mds.stress.DNA
#ordering the stress values for the 100 mds:
order(mds.stress.DNA)
#Saving the order in a vector
ordered.DNA<-order(mds.stress.DNA)
ordered.DNA

#find the stress of the solutions with the lowest and second lowest stress:
mds.stress.DNA[ordered.DNA[1]]
mds.stress.DNA[ordered.DNA[2]]

#scaling of axes to half change units and varimax rotation
mds.best.DNA<-postMDS(mds.DNA[[ordered.DNA[1]]],dist.DNA)
mds.best.DNA
mds.secbest.DNA<-postMDS(mds.DNA[[ordered.DNA[2]]],dist.DNA)
mds.secbest.DNA

#Procrustes comparisons
procrustes(mds.best.DNA,mds.secbest.DNA,permutations=999)
protest(mds.best.DNA,mds.secbest.DNA,permutations=999)
plot(procrustes(mds.best.DNA,mds.secbest.DNA,permutations=999))

#making variables from GNMDS axes 1 and 2 for plotting
gnmds1_DNA<-mds.best.DNA$points[,1]
gnmds2_DNA<-mds.best.DNA$points[,2]

##Biplot ##both using origin=TRUE here
mds.df.DNA<-data.frame(gnmds1_DNA,gnmds2_DNA)
fit.DNA<-envfit(mds.df.DNA, env_dna[,c(1,7, 8, 9, 10, 11, 13, 14)], 999)
fit.DNA

# gnmds1_DNA gnmds2_DNA     r2 Pr(>r)    
# pH             -0.98264    0.18551 0.7694  0.001 ***
# moisture        0.87636    0.48166 0.1204  0.245    
# conductivity   -0.47767    0.87854 0.2193  0.062 .  
# OM             -0.43856    0.89870 0.1215  0.256    
# N              -0.99725   -0.07413 0.1733  0.144    
# C              -0.81961    0.57293 0.1305  0.239    
# cn_ratio        0.50303    0.86427 0.3913  0.011 *   

#Lumped treatments:
par(mfrow=c(1,1))
plot(main="GNMDS_DNA", gnmds1_DNA,gnmds2_DNA,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)
# points(gnmds1_DNA[treatment=="DNA_H_A"], gnmds2_DNA[treatment=="DNA_H_A"], pch=19, col='black') 
# points(gnmds1_DNA[treatment=="DNA_H_D"], gnmds2_DNA[treatment=="DNA_H_D"], pch=21, col='black')
# points(gnmds1_DNA[treatment=="DNA_M_A"], gnmds2_DNA[treatment=="DNA_M_A"], pch=19, col='grey') 
# points(gnmds1_DNA[treatment=="DNA_M_D"], gnmds2_DNA[treatment=="DNA_M_D"], pch=21, col='grey')
points(gnmds1_DNA[lumped=="DNA_H"], gnmds2_DNA[lumped=="DNA_H"], pch=19, col='black') 
points(gnmds1_DNA[lumped=="DNA_M"], gnmds2_DNA[lumped=="DNA_M"], pch=21, col='black')
ordiellipse(mds.df.DNA, env_dna$vegetation_type, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "red")
legend("topleft", legend = c("heath (DNA)", "meadow (DNA)") , pch=c(19,21), col=c("black", "black"), bty="n")

detach(otus_DNA)
detach(env_dna)


###################################################
###   GNMDS for RNA presence absence results    ###
###################################################

#subsampling RNA samples from presence-absence dataset
otus_RNA <- data.frame(otus_pa[24:42])

#otu tables are already prepared in the working space
attach(otus_RNA)
names(otus_DNA)
str(otus_DNA)
attach(env_rna)

#making Bray-Curtis dissimilarity matrix:
dist.RNA<-vegdist(t(otus_RNA), method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
dist.RNA 

#define a general, empty object called mds:
mds.RNA<-NULL
#making 100 "mds"s from initial starting configurations, allocating them into the mds object:
for(i in 1:100)
{mds.RNA[[i]]<-isoMDS(dist.RNA,initMDS(dist.RNA, k=2), k=2, maxit=1000,tol=1e-7)}
#k determines the number of dimensions in the ordination
#the mds object is now a list consisting of 100 "subobjects" being lists
# extracting the stress values as a vector - a measure of the correspondence between dissimilarities
# and distance for each ordination:
mds.stress.RNA<-unlist(lapply(mds.RNA,function(v){v[[2]]})) 

dist.RNA.clust<-hclust(dist.RNA,"single")      # hierarchical clustering 
plot(dist.RNA.clust,ylab="Dissimilarity",xlab="",sub="")  # plot dendrogram

# looking at the stress values for 100 mds:
mds.stress.RNA
#ordering the stress values for the 100 mds:
order(mds.stress.RNA)
#Saving the order in a vector
ordered.RNA<-order(mds.stress.RNA)
ordered.RNA

#find the stress of the solutions with the lowest and second lowest stress:
mds.stress.RNA[ordered.RNA[1]]
mds.stress.RNA[ordered.RNA[2]]

#scaling of axes to half change units and varimax rotation
mds.best.RNA<-postMDS(mds.RNA[[ordered.RNA[1]]],dist.RNA)
mds.best.RNA
mds.secbest.RNA<-postMDS(mds.RNA[[ordered.RNA[2]]],dist.RNA)
mds.secbest.RNA

#Procrustes comparisons
procrustes(mds.best.RNA,mds.secbest.RNA,permutations=999)
protest(mds.best.RNA,mds.secbest.RNA,permutations=999)
plot(procrustes(mds.best.RNA,mds.secbest.RNA,permutations=999))

#making variables from GNMDS axes 1 and 2 for plotting
gnmds1_RNA<-mds.best.RNA$points[,1]
gnmds2_RNA<-mds.best.RNA$points[,2]

##Biplot ##both using origin=TRUE here
mds.df.RNA<-data.frame(gnmds1_RNA,gnmds2_RNA)
fit.RNA<-envfit(mds.df.RNA, env_rna[,c(1,7, 8, 9, 10, 11, 13, 14)], 999)
fit.RNA

# ***VECTORS
# 
# gnmds1_RNA gnmds2_RNA     r2 Pr(>r)    
# pH             -0.99905    0.04349 0.6736  0.001 ***
# moisture       -0.30575    0.95211 0.1640  0.248    
# conductivity   -0.47014    0.88259 0.3560  0.028 *  
# OM             -0.83109    0.55614 0.6919  0.001 ***
# N              -0.75287    0.65817 0.6090  0.001 ***
# C              -0.69783    0.71627 0.6360  0.001 ***
# cn_ratio        0.71228    0.70190 0.0290  0.790      


#Lumped treatments:
par(mfrow=c(1,1))
plot(main="GNMDS_RNA", gnmds1_RNA,gnmds2_RNA,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)
# points(gnmds1_RNA[treatment=="RNA_H_A"], gnmds2_RNA[treatment=="RNA_H_A"], pch=19, col='black') 
# points(gnmds1_RNA[treatment=="RNA_H_D"], gnmds2_RNA[treatment=="RNA_H_D"], pch=21, col='black')
# points(gnmds1_RNA[treatment=="RNA_M_A"], gnmds2_RNA[treatment=="RNA_M_A"], pch=19, col='grey') 
# points(gnmds1_RNA[treatment=="RNA_M_D"], gnmds2_RNA[treatment=="RNA_M_D"], pch=21, col='grey')
points(gnmds1_RNA[lumped=="RNA_H"], gnmds2_RNA[lumped=="RNA_H"], pch=19, col='grey') 
points(gnmds1_RNA[lumped=="RNA_M"], gnmds2_RNA[lumped=="RNA_M"], pch=21, col='grey')
ordiellipse(mds.df.RNA, env_rna$vegetation_type, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "red")
legend("topleft", legend = c("heath (RNA)", "meadow (RNA)") , pch=c(19,21), col=c("grey", "grey"), bty="n")

detach(otus_RNA)
detach(env_rna)



###################################################################
###   GNMDS for symbiotrophs in DNA presence absence results    ###
###################################################################

#subsampling DNA samples from presence-absence dataset
otus_symbio <- data.frame(otus_symbio)
otus_symbio_DNA <- otus_symbio[1:23,]
otus_symbio_DNA[otus_symbio_DNA>0] <-1

#otu tables are already prepared in the working space
attach(otus_symbio_DNA)
names(otus_symbio_DNA)
str(otus_symbio_DNA)
attach(env_dna)

#making Bray-Curtis dissimilarity matrix:
dist.symbio_DNA<-vegdist(t(otus_DNA), method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
dist.symbio_DNA 

#define a general, empty object called mds:
mds.symbio_DNA<-NULL
#making 100 "mds"s from initial starting configurations, allocating them into the mds object:
for(i in 1:100)
{mds.symbio_DNA[[i]]<-isoMDS(dist.symbio_DNA,initMDS(dist.symbio_DNA, k=2), k=2, maxit=1000,tol=1e-7)}
#k determines the number of dimensions in the ordination
#the mds object is now a list consisting of 100 "subobjects" being lists
# extracting the stress values as a vector - a measure of the correspondence between dissimilarities
# and distance for each ordination:
mds.stress.symbio_DNA<-unlist(lapply(mds.symbio_DNA,function(v){v[[2]]})) 

dist.symbio_DNA.clust<-hclust(dist.symbio_DNA,"single")      # hierarchical clustering 
plot(dist.symbio_DNA.clust,ylab="Dissimilarity",xlab="",sub="")  # plot dendrogram

# looking at the stress values for 100 mds:
mds.stress.symbio_DNA
#ordering the stress values for the 100 mds:
order(mds.stress.symbio_DNA)
#Saving the order in a vector
ordered.symbio_DNA<-order(mds.stress.symbio_DNA)
ordered.symbio_DNA

#find the stress of the solutions with the lowest and second lowest stress:
mds.stress.symbio_DNA[ordered.symbio_DNA[1]]
mds.stress.symbio_DNA[ordered.symbio_DNA[2]]

#scaling of axes to half change units and varimax rotation
mds.best.symbio_DNA<-postMDS(mds.symbio_DNA[[ordered.symbio_DNA[1]]],dist.symbio_DNA)
mds.best.symbio_DNA
mds.secbest.symbio_DNA<-postMDS(mds.symbio_DNA[[ordered.symbio_DNA[2]]],dist.symbio_DNA)
mds.secbest.symbio_DNA

#Procrustes comparisons
procrustes(mds.best.symbio_DNA,mds.secbest.symbio_DNA,permutations=999)
protest(mds.best.symbio_DNA,mds.secbest.symbio_DNA,permutations=999)
plot(procrustes(mds.best.symbio_DNA,mds.secbest.symbio_DNA,permutations=999))

#making variables from GNMDS axes 1 and 2 for plotting
gnmds1_symbio_DNA<-mds.best.symbio_DNA$points[,1]
gnmds2_symbio_DNA<-mds.best.symbio_DNA$points[,2]

##Biplot ##both using origin=TRUE here
mds.df.symbio_DNA<-data.frame(gnmds1_symbio_DNA,gnmds2_symbio_DNA)
fit.symbio_DNA<-envfit(mds.df.symbio_DNA, env_dna[,c(1,7, 8, 9, 10, 11, 13, 14)], 999)
fit.symbio_DNA

# gnmds1_symbio_DNA gnmds2_symbio_DNA     r2 Pr(>r)    
# pH                    -0.98264           0.18552 0.7694  0.001 ***
# moisture               0.87631           0.48174 0.1204  0.282    
# conductivity          -0.47767           0.87854 0.2193  0.063 .  
# OM                    -0.43853           0.89872 0.1215  0.263    
# N                     -0.99726          -0.07394 0.1733  0.150    
# C                     -0.81954           0.57302 0.1305  0.248    
# cn_ratio               0.50308           0.86424 0.3913  0.010 **     
 

#Lumped treatments:
par(mfrow=c(1,1))
plot(main="GNMDS_symbio_DNA", gnmds1_symbio_DNA,gnmds2_symbio_DNA,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)
#points(gnmds1_symbio_DNA[treatment=="DNA_H_A"], gnmds2_symbio_DNA[treatment=="DNA_H_A"], pch=19, col='black') 
#points(gnmds1_symbio_DNA[treatment=="DNA_H_D"], gnmds2_symbio_DNA[treatment=="DNA_H_D"], pch=21, col='black')
#points(gnmds1_symbio_DNA[treatment=="DNA_M_A"], gnmds2_symbio_DNA[treatment=="DNA_M_A"], pch=19, col='grey') 
#points(gnmds1_symbio_DNA[treatment=="DNA_M_D"], gnmds2_symbio_DNA[treatment=="DNA_M_D"], pch=21, col='grey')
points(gnmds1_symbio_DNA[lumped=="DNA_H"], gnmds2_symbio_DNA[lumped=="DNA_H"], pch=19, col='black') 
points(gnmds1_symbio_DNA[lumped=="DNA_M"], gnmds2_symbio_DNA[lumped=="DNA_M"], pch=21, col='black')
ordiellipse(mds.df.symbio_DNA, env_dna$vegetation_type, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "red")
legend("topleft", legend = c("heath (DNA)", "meadow (DNA)") , pch=c(19,21), col=c("black", "black"), bty="n")

detach(otus_symbio_DNA)
detach(env_dna)


###################################################################
###   GNMDS for symbiotrophs in RNA presence absence results    ###
###################################################################


#subsampling DNA samples from presence-absence dataset
otus_symbio_RNA <- otus_symbio[24:42,]
otus_symbio_DNA[otus_symbio_DNA>0] <-1

#otu tables are already prepared in the working space
attach(otus_symbio_RNA)
names(otus_symbio_RNA)
str(otus_symbio_RNA)
attach(env_rna)

#making Bray-Curtis dissimilarity matrix:
dist.symbio_RNA<-vegdist(otus_symbio_RNA, method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
dist.symbio_RNA

#define a general, empty object called mds:
mds.symbio_RNA<-NULL
#making 100 "mds"s from initial starting configurations, allocating them into the mds object:
for(i in 1:100)
{mds.symbio_RNA[[i]]<-isoMDS(dist.symbio_RNA,initMDS(dist.symbio_RNA, k=2), k=2, maxit=1000,tol=1e-7)}
#k determines the number of dimensions in the ordination
#the mds object is now a list consisting of 100 "subobjects" being lists
# extracting the stress values as a vector - a measure of the correspondence between dissimilarities
# and distance for each ordination:
mds.stress.symbio_RNA<-unlist(lapply(mds.symbio_RNA,function(v){v[[2]]})) 

dist.symbio_RNA.clust<-hclust(dist.symbio_RNA,"single")      # hierarchical clustering 
plot(dist.symbio_RNA.clust,ylab="Dissimilarity",xlab="",sub="")  # plot dendrogram

# looking at the stress values for 100 mds:
mds.stress.symbio_RNA
#ordering the stress values for the 100 mds:
order(mds.stress.symbio_RNA)
#Saving the order in a vector
ordered.symbio_RNA<-order(mds.stress.symbio_RNA)
ordered.symbio_RNA

#find the stress of the solutions with the lowest and second lowest stress:
mds.stress.symbio_RNA[ordered.symbio_RNA[1]]
mds.stress.symbio_RNA[ordered.symbio_RNA[2]]

#scaling of axes to half change units and varimax rotation
mds.best.symbio_RNA<-postMDS(mds.symbio_RNA[[ordered.symbio_RNA[1]]],dist.symbio_RNA)
mds.best.symbio_RNA
mds.secbest.symbio_RNA<-postMDS(mds.symbio_RNA[[ordered.symbio_RNA[2]]],dist.symbio_RNA)
mds.secbest.symbio_RNA

#Procrustes comparisons
procrustes(mds.best.symbio_RNA,mds.secbest.symbio_RNA,permutations=999)
protest(mds.best.symbio_RNA,mds.secbest.symbio_RNA,permutations=999)
plot(procrustes(mds.best.symbio_RNA,mds.secbest.symbio_RNA,permutations=999))

#making variables from GNMDS axes 1 and 2 for plotting
gnmds1_symbio_RNA<-mds.best.symbio_RNA$points[,1]
gnmds2_symbio_RNA<-mds.best.symbio_RNA$points[,2]

##Biplot ##both using origin=TRUE here
mds.df.symbio_RNA<-data.frame(gnmds1_symbio_RNA,gnmds2_symbio_RNA)
fit.symbio_RNA<-envfit(mds.df.symbio_RNA, env_rna[,c(1,7, 8, 9, 10, 11, 13, 14)], 999)
fit.symbio_RNA

# gnmds1_symbio_RNA gnmds2_symbio_RNA     r2 Pr(>r)    
# pH                    -0.64995          -0.75998 0.6302  0.002 ** 
# moisture               0.17635          -0.98433 0.1264  0.342    
# conductivity          -0.07671          -0.99705 0.3632  0.018 *  
# OM                    -0.30455          -0.95250 0.6475  0.001 ***
# N                     -0.31215          -0.95003 0.5040  0.005 ** 
# C                     -0.25140          -0.96788 0.4925  0.007 ** 
# cn_ratio               0.97711           0.21274 0.0611  0.621 

#Lumped treatments:
par(mfrow=c(1,1))
plot(main="GNMDS_symbio_RNA", gnmds1_symbio_RNA,gnmds2_symbio_RNA,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)
# points(gnmds1_symbio_RNA[treatment=="RNA_H_A"], gnmds2_symbio_RNA[treatment=="RNA_H_A"], pch=19, col='black') 
# points(gnmds1_symbio_RNA[treatment=="RNA_H_D"], gnmds2_symbio_RNA[treatment=="RNA_H_D"], pch=21, col='black')
# points(gnmds1_symbio_RNA[treatment=="RNA_M_A"], gnmds2_symbio_RNA[treatment=="RNA_M_A"], pch=19, col='grey') 
# points(gnmds1_symbio_RNA[treatment=="RNA_M_D"], gnmds2_symbio_RNA[treatment=="RNA_M_D"], pch=21, col='grey')
points(gnmds1_symbio_RNA[lumped=="RNA_H"], gnmds2_symbio_RNA[lumped=="RNA_H"], pch=19, col='grey') 
points(gnmds1_symbio_RNA[lumped=="RNA_M"], gnmds2_symbio_RNA[lumped=="RNA_M"], pch=21, col='grey')
ordiellipse(mds.df.symbio_RNA, env_rna$vegetation_type, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "red")
legend("topleft", legend = c("heath (RNA)", "meadow (RNA)") , pch=c(19,21), col=c("grey", "grey"), bty="n")

detach(otus_symbio_RNA)
detach(env_rna)

###################################################################
###   GNMDS for SAPROTROPHS in DNA presence absence results    ###
###################################################################

#subsampling DNA samples from presence-absence dataset
otus_sapro <- data.frame(otus_sapro)
otus_sapro_DNA <- otus_sapro[1:23,]
otus_sapro_DNA[otus_sapro_DNA>0] <-1

#otu tables are already prepared in the working space
attach(otus_sapro_DNA)
names(otus_sapro_DNA)
str(otus_sapro_DNA)
attach(env_dna)

#making Bray-Curtis dissimilarity matrix:
dist.sapro_DNA<-vegdist(otus_sapro_DNA, method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
dist.sapro_DNA 

#define a general, empty object called mds:
mds.sapro_DNA<-NULL
#making 100 "mds"s from initial starting configurations, allocating them into the mds object:
for(i in 1:100)
{mds.sapro_DNA[[i]]<-isoMDS(dist.sapro_DNA,initMDS(dist.sapro_DNA, k=2), k=2, maxit=1000,tol=1e-7)}
#k determines the number of dimensions in the ordination
#the mds object is now a list consisting of 100 "subobjects" being lists
# extracting the stress values as a vector - a measure of the correspondence between dissimilarities
# and distance for each ordination:
mds.stress.sapro_DNA<-unlist(lapply(mds.sapro_DNA,function(v){v[[2]]})) 

dist.sapro_DNA.clust<-hclust(dist.sapro_DNA,"single")      # hierarchical clustering 
plot(dist.sapro_DNA.clust,ylab="Dissimilarity",xlab="",sub="")  # plot dendrogram

# looking at the stress values for 100 mds:
mds.stress.sapro_DNA
#ordering the stress values for the 100 mds:
order(mds.stress.sapro_DNA)
#Saving the order in a vector
ordered.sapro_DNA<-order(mds.stress.sapro_DNA)
ordered.sapro_DNA

#find the stress of the solutions with the lowest and second lowest stress:
mds.stress.sapro_DNA[ordered.sapro_DNA[1]]
mds.stress.sapro_DNA[ordered.sapro_DNA[2]]

#scaling of axes to half change units and varimax rotation
mds.best.sapro_DNA<-postMDS(mds.sapro_DNA[[ordered.sapro_DNA[1]]],dist.sapro_DNA)
mds.best.sapro_DNA
mds.secbest.sapro_DNA<-postMDS(mds.sapro_DNA[[ordered.sapro_DNA[2]]],dist.sapro_DNA)
mds.secbest.sapro_DNA

#Procrustes comparisons
procrustes(mds.best.sapro_DNA,mds.secbest.sapro_DNA,permutations=999)
protest(mds.best.sapro_DNA,mds.secbest.sapro_DNA,permutations=999)
plot(procrustes(mds.best.sapro_DNA,mds.secbest.sapro_DNA,permutations=999))

#making variables from GNMDS axes 1 and 2 for plotting
gnmds1_sapro_DNA<-mds.best.sapro_DNA$points[,1]
gnmds2_sapro_DNA<-mds.best.sapro_DNA$points[,2]

##Biplot ##both using origin=TRUE here
mds.df.sapro_DNA<-data.frame(gnmds1_sapro_DNA,gnmds2_sapro_DNA)
fit.sapro_DNA<-envfit(mds.df.sapro_DNA, env_dna[,c(1,7, 8, 9, 10, 11, 13, 14)], 999)
fit.sapro_DNA

# gnmds1_sapro_DNA gnmds2_sapro_DNA     r2 Pr(>r)    
# pH                   -0.99224         -0.12435 0.5750  0.001 ***
# moisture              0.53363          0.84572 0.1024  0.340    
# conductivity         -0.79652          0.60462 0.1045  0.347    
# OM                   -0.82112          0.57075 0.1029  0.323    
# N                    -0.99198         -0.12638 0.2494  0.061 .  
# C                    -0.97531          0.22086 0.1892  0.114    
# cn_ratio              0.80126          0.59832 0.2749  0.043 *

#Lumped treatments:
par(mfrow=c(1,1))
plot(main="GNMDS_sapro_DNA", gnmds1_sapro_DNA,gnmds2_sapro_DNA,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)
# points(gnmds1_sapro_DNA[treatment=="DNA_H_A"], gnmds2_sapro_DNA[treatment=="DNA_H_A"], pch=19, col='black') 
# points(gnmds1_sapro_DNA[treatment=="DNA_H_D"], gnmds2_sapro_DNA[treatment=="DNA_H_D"], pch=21, col='black')
# points(gnmds1_sapro_DNA[treatment=="DNA_M_A"], gnmds2_sapro_DNA[treatment=="DNA_M_A"], pch=19, col='grey') 
# points(gnmds1_sapro_DNA[treatment=="DNA_M_D"], gnmds2_sapro_DNA[treatment=="DNA_M_D"], pch=21, col='grey')
points(gnmds1_sapro_DNA[lumped=="DNA_H"], gnmds2_sapro_DNA[lumped=="DNA_H"], pch=19, col='black') 
points(gnmds1_sapro_DNA[lumped=="DNA_M"], gnmds2_sapro_DNA[lumped=="DNA_M"], pch=21, col='black')
ordiellipse(mds.df.sapro_DNA, env_dna$vegetation_type, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "red")
legend("topleft", legend = c("heath (DNA)", "meadow (DNA)") , pch=c(19,21), col=c("black", "black"), bty="n")

detach(otus_sapro_DNA)
detach(env_dna)


###################################################################
###   GNMDS for SAPROTROPHS in RNA presence absence results    ###
###################################################################
#subsampling RNA samples from presence-absence dataset
otus_sapro_RNA <- otus_sapro[24:42,]
otus_sapro_RNA[otus_sapro_RNA>0] <-1

#otu tables are already prepared in the working space
attach(otus_sapro_RNA)
names(otus_sapro_RNA)
str(otus_sapro_RNA)
attach(env_rna)

#making Bray-Curtis dissimilarity matrix:
dist.sapro_RNA<-vegdist(otus_sapro_RNA, method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
dist.sapro_RNA

#define a general, empty object called mds:
mds.sapro_RNA<-NULL
#making 100 "mds"s from initial starting configurations, allocating them into the mds object:
for(i in 1:100)
{mds.sapro_RNA[[i]]<-isoMDS(dist.sapro_RNA,initMDS(dist.sapro_RNA, k=2), k=2, maxit=1000,tol=1e-7)}
#k determines the number of dimensions in the ordination
#the mds object is now a list consisting of 100 "subobjects" being lists
# extracting the stress values as a vector - a measure of the correspondence between dissimilarities
# and distance for each ordination:
mds.stress.sapro_RNA<-unlist(lapply(mds.sapro_RNA,function(v){v[[2]]})) 

dist.sapro_RNA.clust<-hclust(dist.sapro_RNA,"single")      # hierarchical clustering 
plot(dist.sapro_RNA.clust,ylab="Dissimilarity",xlab="",sub="")  # plot dendrogram

# looking at the stress values for 100 mds:
mds.stress.sapro_RNA
#ordering the stress values for the 100 mds:
order(mds.stress.sapro_RNA)
#Saving the order in a vector
ordered.sapro_RNA<-order(mds.stress.sapro_RNA)
ordered.sapro_RNA

#find the stress of the solutions with the lowest and second lowest stress:
mds.stress.sapro_RNA[ordered.sapro_RNA[1]]
mds.stress.sapro_RNA[ordered.sapro_RNA[2]]

#scaling of axes to half change units and varimax rotation
mds.best.sapro_RNA<-postMDS(mds.sapro_RNA[[ordered.sapro_RNA[1]]],dist.sapro_RNA)
mds.best.sapro_RNA
mds.secbest.sapro_RNA<-postMDS(mds.sapro_RNA[[ordered.sapro_RNA[2]]],dist.sapro_RNA)
mds.secbest.sapro_RNA

#Procrustes comparisons
procrustes(mds.best.sapro_RNA,mds.secbest.sapro_RNA,permutations=999)
protest(mds.best.sapro_RNA,mds.secbest.sapro_RNA,permutations=999)
plot(procrustes(mds.best.sapro_RNA,mds.secbest.sapro_RNA,permutations=999))

#making variables from GNMDS axes 1 and 2 for plotting
gnmds1_sapro_RNA<-mds.best.sapro_RNA$points[,1]
gnmds2_sapro_RNA<-mds.best.sapro_RNA$points[,2]

##Biplot ##both using origin=TRUE here
mds.df.sapro_RNA<-data.frame(gnmds1_sapro_RNA,gnmds2_sapro_RNA)
fit.sapro_RNA<-envfit(mds.df.sapro_RNA, env_rna[,c(1, 7, 8, 9, 10, 11, 13, 14)], 999)
fit.sapro_RNA

# gnmds1_sapro_RNA gnmds2_sapro_RNA     r2 Pr(>r)    
# pH                    0.88994          0.45607 0.6018  0.001 ***
# moisture              0.99585          0.09102 0.1674  0.237    
# conductivity          0.89240          0.45124 0.1148  0.384    
# OM                    0.94193          0.33581 0.5765  0.002 ** 
# N                     0.95540          0.29531 0.4328  0.015 *  
# C                     0.97385          0.22721 0.5453  0.003 ** 
# cn_ratio              0.48022         -0.87715 0.0589  0.598

#Lumped treatments:
par(mfrow=c(1,1))
plot(main="GNMDS_sapro_RNA", gnmds1_sapro_RNA,gnmds2_sapro_RNA,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)
# points(gnmds1_sapro_RNA[treatment=="RNA_H_A"], gnmds2_sapro_RNA[treatment=="RNA_H_A"], pch=19, col='black') 
# points(gnmds1_sapro_RNA[treatment=="RNA_H_D"], gnmds2_sapro_RNA[treatment=="RNA_H_D"], pch=21, col='black')
# points(gnmds1_sapro_RNA[treatment=="RNA_M_A"], gnmds2_sapro_RNA[treatment=="RNA_M_A"], pch=19, col='grey') 
# points(gnmds1_sapro_RNA[treatment=="RNA_M_D"], gnmds2_sapro_RNA[treatment=="RNA_M_D"], pch=21, col='grey')
points(gnmds1_sapro_RNA[lumped=="RNA_H"], gnmds2_sapro_RNA[lumped=="RNA_H"], pch=19, col='grey') 
points(gnmds1_sapro_RNA[lumped=="RNA_M"], gnmds2_sapro_RNA[lumped=="RNA_M"], pch=21, col='grey')
ordiellipse(mds.df.sapro_RNA, env_rna$vegetation_type, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "red")
legend("topleft", legend = c("heath (RNA)", "meadow (RNA)") , pch=c(19,21), col=c("grey", "grey"), bty="n")

detach(otus_sapro_RNA)
detach(env_rna)

###################################################################
###   GNMDS for PATHOTROPHS  in DNA presence absence results    ###
###################################################################

#subsampling DNA samples from presence-absence dataset
otus_patho <- data.frame(otus_patho)
otus_patho_DNA <- otus_patho[,1:23]
otus_patho_DNA[otus_patho_DNA>0] <-1

#otu tables are already prepared in the working space
attach(otus_patho_DNA)
names(otus_patho_DNA)
str(otus_patho_DNA)
attach(env_dna)

#making Bray-Curtis dissimilarity matrix:
dist.patho_DNA<-vegdist(t(otus_patho_DNA), method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
dist.patho_DNA 

#define a general, empty object called mds:
mds.patho_DNA<-NULL
#making 100 "mds"s from initial starting configurations, allocating them into the mds object:
for(i in 1:100)
{mds.patho_DNA[[i]]<-isoMDS(dist.patho_DNA,initMDS(dist.patho_DNA, k=2), k=2, maxit=1000,tol=1e-7)}
#k determines the number of dimensions in the ordination
#the mds object is now a list consisting of 100 "subobjects" being lists
# extracting the stress values as a vector - a measure of the correspondence between dissimilarities
# and distance for each ordination:
mds.stress.patho_DNA<-unlist(lapply(mds.patho_DNA,function(v){v[[2]]})) 

dist.patho_DNA.clust<-hclust(dist.patho_DNA,"single")      # hierarchical clustering 
plot(dist.patho_DNA.clust,ylab="Dissimilarity",xlab="",sub="")  # plot dendrogram

# looking at the stress values for 100 mds:
mds.stress.patho_DNA
#ordering the stress values for the 100 mds:
order(mds.stress.patho_DNA)
#Saving the order in a vector
ordered.patho_DNA<-order(mds.stress.patho_DNA)
ordered.patho_DNA

#find the stress of the solutions with the lowest and second lowest stress:
mds.stress.patho_DNA[ordered.patho_DNA[1]]
mds.stress.patho_DNA[ordered.patho_DNA[2]]

#scaling of axes to half change units and varimax rotation
mds.best.patho_DNA<-postMDS(mds.patho_DNA[[ordered.patho_DNA[1]]],dist.patho_DNA)
mds.best.patho_DNA
mds.secbest.patho_DNA<-postMDS(mds.patho_DNA[[ordered.patho_DNA[2]]],dist.patho_DNA)
mds.secbest.patho_DNA

#Procrustes comparisons
procrustes(mds.best.patho_DNA,mds.secbest.patho_DNA,permutations=999)
protest(mds.best.patho_DNA,mds.secbest.patho_DNA,permutations=999)
plot(procrustes(mds.best.patho_DNA,mds.secbest.patho_DNA,permutations=999))

#making variables from GNMDS axes 1 and 2 for plotting
gnmds1_patho_DNA<-mds.best.patho_DNA$points[,1]
gnmds2_patho_DNA<-mds.best.patho_DNA$points[,2]

##Biplot ##both using origin=TRUE here
mds.df.patho_DNA<-data.frame(gnmds1_patho_DNA,gnmds2_patho_DNA)
fit.patho_DNA<-envfit(mds.df.patho_DNA, env_dna[,c(1, 7, 8, 9, 10, 11, 13, 14)], 999)
fit.patho_DNA

# gnmds1_patho_DNA gnmds2_patho_DNA     r2 Pr(>r)
# pH                    0.95321         -0.30231 0.1424  0.218
# moisture             -0.94692          0.32147 0.0732  0.452
# conductivity          0.97153          0.23693 0.0088  0.910
# OM                    0.95799         -0.28681 0.0063  0.918
# N                     0.94503         -0.32697 0.0037  0.963
# C                     0.22732         -0.97382 0.0093  0.899
# cn_ratio             -0.37961         -0.92514 0.0309  0.733

#Lumped treatments:
par(mfrow=c(1,1))
plot(main="GNMDS_patho_DNA", gnmds1_patho_DNA,gnmds2_patho_DNA,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)
# points(gnmds1_patho_DNA[treatment=="DNA_H_A"], gnmds2_patho_DNA[treatment=="DNA_H_A"], pch=19, col='black') 
# points(gnmds1_patho_DNA[treatment=="DNA_H_D"], gnmds2_patho_DNA[treatment=="DNA_H_D"], pch=21, col='black')
# points(gnmds1_patho_DNA[treatment=="DNA_M_A"], gnmds2_patho_DNA[treatment=="DNA_M_A"], pch=19, col='grey') 
# points(gnmds1_patho_DNA[treatment=="DNA_M_D"], gnmds2_patho_DNA[treatment=="DNA_M_D"], pch=21, col='grey')
points(gnmds1_patho_DNA[lumped=="DNA_H"], gnmds2_patho_DNA[lumped=="DNA_H"], pch=19, col='black') 
points(gnmds1_patho_DNA[lumped=="DNA_M"], gnmds2_patho_DNA[lumped=="DNA_M"], pch=21, col='black')
ordiellipse(mds.df.patho_DNA, env_dna$vegetation_type, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "red")
legend("topleft", legend = c("heath (DNA)", "meadow (DNA)") , pch=c(19,21), col=c("black", "black"), bty="n")

detach(otus_patho_DNA)
detach(env_dna)


###################################################################
###   GNMDS for PATHOTROPHS in RNA presence absence results    ###
###################################################################
#subsampling RNA samples from presence-absence dataset
otus_patho_RNA <- otus_patho[,24:42]
otus_patho_RNA[otus_patho_RNA>0] <-1

#otu tables are already prepared in the working space
attach(otus_patho_RNA)
names(otus_patho_RNA)
str(otus_patho_RNA)
attach(env_rna)

#making Bray-Curtis dissimilarity matrix:
dist.patho_RNA<-vegdist(t(otus_patho_RNA), method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
dist.patho_RNA

#define a general, empty object called mds:
mds.patho_RNA<-NULL
#making 100 "mds"s from initial starting configurations, allocating them into the mds object:
for(i in 1:100)
{mds.patho_RNA[[i]]<-isoMDS(dist.patho_RNA,initMDS(dist.patho_RNA, k=2), k=2, maxit=1000,tol=1e-7)}
#k determines the number of dimensions in the ordination
#the mds object is now a list consisting of 100 "subobjects" being lists
# extracting the stress values as a vector - a measure of the correspondence between dissimilarities
# and distance for each ordination:
mds.stress.patho_RNA<-unlist(lapply(mds.patho_RNA,function(v){v[[2]]})) 

dist.patho_RNA.clust<-hclust(dist.patho_RNA,"single")      # hierarchical clustering 
plot(dist.patho_RNA.clust,ylab="Dissimilarity",xlab="",sub="")  # plot dendrogram

# looking at the stress values for 100 mds:
mds.stress.patho_RNA
#ordering the stress values for the 100 mds:
order(mds.stress.patho_RNA)
#Saving the order in a vector
ordered.patho_RNA<-order(mds.stress.patho_RNA)
ordered.patho_RNA

#find the stress of the solutions with the lowest and second lowest stress:
mds.stress.patho_RNA[ordered.patho_RNA[1]]
mds.stress.patho_RNA[ordered.patho_RNA[2]]

#scaling of axes to half change units and varimax rotation
mds.best.patho_RNA<-postMDS(mds.patho_RNA[[ordered.patho_RNA[1]]],dist.patho_RNA)
mds.best.patho_RNA
mds.secbest.patho_RNA<-postMDS(mds.patho_RNA[[ordered.patho_RNA[2]]],dist.patho_RNA)
mds.secbest.patho_RNA

#Procrustes comparisons
procrustes(mds.best.patho_RNA,mds.secbest.patho_RNA,permutations=999)
protest(mds.best.patho_RNA,mds.secbest.patho_RNA,permutations=999)
plot(procrustes(mds.best.patho_RNA,mds.secbest.patho_RNA,permutations=999))

#making variables from GNMDS axes 1 and 2 for plotting
gnmds1_patho_RNA<-mds.best.patho_RNA$points[,1]
gnmds2_patho_RNA<-mds.best.patho_RNA$points[,2]

##Biplot ##both using origin=TRUE here
mds.df.patho_RNA<-data.frame(gnmds1_patho_RNA,gnmds2_patho_RNA)
fit.patho_RNA<-envfit(mds.df.patho_RNA, env_rna[,c(1, 7, 8, 9, 10, 11, 13, 14)], 999)
fit.patho_RNA

# gnmds1_patho_RNA gnmds2_patho_RNA     r2 Pr(>r)
# pH                    0.12636         -0.99198 0.0309  0.779
# moisture             -0.90322          0.42918 0.0087  0.946
# conductivity          0.38394          0.92336 0.0834  0.486
# OM                   -0.02816         -0.99960 0.0404  0.736
# N                     0.07268         -0.99736 0.0370  0.771
# C                     0.19256         -0.98128 0.0299  0.806
# cn_ratio              0.48733         -0.87322 0.0122  0.914

#Lumped treatments:
par(mfrow=c(1,1))
plot(main="GNMDS_patho_RNA", gnmds1_patho_RNA,gnmds2_patho_RNA,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)
# points(gnmds1_patho_RNA[treatment=="RNA_H_A"], gnmds2_patho_RNA[treatment=="RNA_H_A"], pch=19, col='black') 
# points(gnmds1_patho_RNA[treatment=="RNA_H_D"], gnmds2_patho_RNA[treatment=="RNA_H_D"], pch=21, col='black')
# points(gnmds1_patho_RNA[treatment=="RNA_M_A"], gnmds2_patho_RNA[treatment=="RNA_M_A"], pch=19, col='grey') 
# points(gnmds1_patho_RNA[treatment=="RNA_M_D"], gnmds2_patho_RNA[treatment=="RNA_M_D"], pch=21, col='grey')
points(gnmds1_patho_RNA[lumped=="RNA_H"], gnmds2_patho_RNA[lumped=="RNA_H"], pch=19, col='grey') 
points(gnmds1_patho_RNA[lumped=="RNA_M"], gnmds2_patho_RNA[lumped=="RNA_M"], pch=21, col='grey')
ordiellipse(mds.df.patho_RNA, env_rna$vegetation_type, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "red")
legend("topleft", legend = c("heath (RNA)", "meadow (RNA)") , pch=c(19,21), col=c("grey", "grey"), bty="n")

detach(otus_patho_RNA)
detach(env_rna)


###################################################
####COMPARISON OF ORDINATIONS FOR DNA AND RNA  ####
###################################################

comparison_DNA <- otus_DNA[, c(1, 2, 3, 5, 17, 19, 20, 22, 23)]#chosing only those samples that have both DNA and RNA data
comparison_RNA <- otus_RNA[, c(1, 5, 6, 7, 12, 13, 15, 18, 19)]
##double-checking colnames(comparison_DNA)

#calculating and ploting distances:
dist.comparison_DNA <-vegdist(t(comparison_DNA), method = "jaccard", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)
dist.comparison_RNA <-vegdist(t(comparison_RNA), method = "jaccard", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)
plot(dist.comparison_DNA, dist.comparison_RNA, xlab = "Dissimilarity matrix based on rDNA", ylab = "Disssimilarity matrix based on rRNA", main = "Dissimilarity matrices of coextracted rDNA and rRNA samples", xlim = c(0.5, 1), ylim = c(0.5, 1))

#Mantel's test comparing both distance matrices
library(ade4)
mantel_ade4 <- mantel.rtest(dist.comparison_DNA, dist.comparison_RNA, nrepet = 9999)
mantel_ade4

# Monte-Carlo test
# Call: mantel.rtest(m1 = dist.comparison_DNA, m2 = dist.comparison_RNA, 
#                    nrepet = 9999)
# 
# Observation: 0.7269997 
# 
# Based on 9999 replicates
# Simulated p-value: 1e-04 
# Alternative hypothesis: greater 
# 
# Std.Obs   Expectation      Variance 
# 4.2758248614 -0.0003383711  0.0289356209 

plot(dist.comparison_DNA, dist.comparison_RNA, main = "Dissimilarity matrices of coextracted rDNA and rRNA samples")
plot(mantel_ade4, main = "Mantel's test")



#############################
####PERMANOVA            ####
#############################

library(vegan)
##DNA_all############################################################################
adonis(t(otus_DNA) ~ snow_regime, data=env_dna, method = "jaccard", pemutations=999)
adonis(t(otus_DNA) ~ vegetation_type, data=env_dna, method = "jaccard", pemutations=999)
#adonis(t(otus_DNA) ~ vegetation_type * snow_regime, data=env_dna, pemutations=999, method = "jaccard")
adonis(t(otus_DNA) ~ vegetation_type * snow_regime, data=env_dna, pemutations=999, method = "jaccard", strata = env_dna$block)

##RNA_all#############################################################################
adonis(t(otus_RNA) ~ snow_regime, data=env_rna, method = "jaccard", pemutations=999)
adonis(t(otus_RNA) ~ vegetation_type, data=env_rna, method = "jaccard", pemutations=999)
#adonis(t(otus_RNA) ~ vegetation_type * snow_regime, data=env_rna, pemutations=999, method = "jaccard")
adonis(t(otus_RNA) ~ vegetation_type * snow_regime, data=env_rna, pemutations=999, method = "jaccard", strata = env_rna$block) 

##symbio_DNA###########################################################################################
adonis(otus_symbio_DNA ~ snow_regime, data=env_dna, method = "jaccard", pemutations=999)
adonis(otus_symbio_DNA ~ vegetation_type, data=env_dna, method = "jaccard", pemutations=999)
#adonis(otus_symbio_DNA ~ vegetation_type*snow_regime, data=env_dna, pemutations=999, method = "jaccard")
adonis(otus_symbio_DNA ~ vegetation_type*snow_regime, data=env_dna, pemutations=999, method = "jaccard", strata = env_dna$block)

##symbio_RNA###########################################################################################
adonis(otus_symbio_RNA ~ snow_regime, data=env_rna, method = "jaccard", pemutations=999)
adonis(otus_symbio_RNA ~ vegetation_type, data=env_rna, method = "jaccard", pemutations=999)
#adonis(otus_symbio_RNA ~ vegetation_type * snow_regime, data=env_rna, method = "jaccard", pemutations=999)
adonis(otus_symbio_RNA ~ vegetation_type * snow_regime, data=env_rna, method = "jaccard", pemutations=999, strata=env_rna$block)

##sapro_DNA###########################################################################################
adonis(otus_sapro_DNA ~ snow_regime, data=env_dna, method='jaccard',pemutations=999)
adonis(otus_sapro_DNA ~ vegetation_type, data=env_dna, method='jaccard',pemutations=999)
#adonis(otus_sapro_DNA ~ vegetation_type * snow_regime,data=env_dna,method='jaccard',pemutations=999)
adonis(otus_sapro_DNA ~ vegetation_type * snow_regime,data=env_dna,method='jaccard',pemutations=999,strata = env_dna$block)

##sapro_RNA###########################################################################################
adonis(otus_sapro_RNA ~ snow_regime, data=env_rna, method='jaccard',pemutations=999)
adonis(otus_sapro_RNA ~ vegetation_type, data=env_rna, method='jaccard',pemutations=999)
#adonis(otus_sapro_RNA ~ vegetation_type * snow_regime, data=env_rna, method='jaccard',pemutations=999)
adonis(otus_sapro_RNA ~ vegetation_type * snow_regime, data=env_rna, method='jaccard',pemutations=999, strata = env_rna$block)

##patho_DNA###########################################################################################
adonis(t(otus_patho_DNA) ~ snow_regime, data=env_dna, method='jaccard',pemutations=999)
adonis(t(otus_patho_DNA) ~ vegetation_type, data=env_dna, method='jaccard',pemutations=999)
#adonis(t(otus_patho_DNA) ~ vegetation_type * snow_regime,data=env_dna,method='jaccard',pemutations=999)
adonis(t(otus_patho_DNA) ~ vegetation_type * snow_regime,data=env_dna,method='jaccard',pemutations=999,strata = env_dna$block)

##patho_RNA###########################################################################################
adonis(t(otus_patho_RNA) ~ snow_regime, data=env_rna, method='jaccard',pemutations=999)
adonis(t(otus_patho_RNA) ~ vegetation_type, data=env_rna, method='jaccard',pemutations=999)
#adonis(t(otus_patho_RNA) ~ vegetation_type * snow_regime, data=env_rna, method='jaccard',pemutations=999)
adonis(t(otus_patho_RNA) ~ vegetation_type * snow_regime, data=env_rna, method='jaccard',pemutations=999, strata = env_rna$block)


###################
#Removing outliers#
###################

otus_DNA_nooutliers <- otus_DNA[,c(1:17, 19:23)]
env_dna_nooutliers <- env_dna[c(1:17, 19:23),]
otus_RNA_nooutliers <- otus_RNA[,c(1:13, 15:19)]
env_rna_nooutliers <- env_rna[c(1:13, 15:19),]
env_nooutliers <- env[c(1:17, 19:36, 38:42),]

###############################
####OTU richness and ANOVA ####
###############################

library(car) #This is car 2.1-4
library(lme4) #This is lme4 1.1-13
library(plyr)

#checking for overall differences in OTU richness between rDNA and rRNA
env$richness<-specnumber(otus_t)
meanRichness <- ddply(env, .(nucleic_acid), summarize, mean=mean(richness))
# DNA 102.1739
# RNA 115.3684
sdRichness <- ddply(env, .(nucleic_acid), summarize, sd=sd(richness))
# DNA 24.36401
# RNA 20.47549

#checking for overall differences in OTU richness between rDNA and rRNA with no outliers
meanRichness <- ddply(env_nooutliers, .(nucleic_acid), summarize, mean=mean(richness))
#DNA  98.77273
#RNA 112.94444
sdRichness <- ddply(env_nooutliers, .(nucleic_acid), summarize, sd=sd(richness))
#DNA 18.52394
#RNA 18.04723

#checking for differences in symbiotrophic OTU richness between rDNA and rRNA 
otus_tax <- cbind(otus, tax_fun)
otus_symbio <- subset(otus_tax, otus_tax$Simplified_Trophic_mode=="Symbiotroph") #subsetting symbiotrophs in OTU table
otus_symbio <- otus_symbio[,1:42]
otus_symbio <- t(otus_symbio)
env$richnessSymbio<-specnumber(otus_symbio)
meanRichness <- ddply(env, .(nucleic_acid), summarize, mean=mean(richnessSymbio))
# DNA 35.217392
# RNA 43.47368
sdRichness <- ddply(env, .(nucleic_acid), summarize, sd=sd(richnessSymbio))
#DNA 9.0702142
#RNA 10.007599

#checking for differences in symbiotrophic OTU richness between rDNA and rRNA without outliers
env_nooutliers <- env[c(1:17, 19:36, 38:42),]
meanRichness <- ddply(env_nooutliers, .(nucleic_acid), summarize, mean=mean(richnessSymbio))
#DNA 35.09091
#RNA 42.66667
sdRichness <- ddply(env_nooutliers, .(nucleic_acid), summarize, sd=sd(richnessSymbio))
#DNA 9.262876
#RNA 9.640600


#checking for overall differences in saprotrophic OTU richness between rDNA and rRNA 
otus_sapro <- subset(otus_tax, otus_tax$Simplified_Trophic_mode=="Saprotroph") #subsetting saprotrophs in OTU table
otus_sapro <- otus_sapro[,1:42]
otus_sapro <- t(otus_sapro)
env$richnessSapro<-specnumber(otus_sapro)
meanRichness <- ddply(env, .(nucleic_acid), summarize, mean=mean(richnessSapro))
# DNA 13.956522
# RNA 15.00000
sdRichness <- ddply(env, .(nucleic_acid), summarize, sd=sd(richnessSapro))
# DNA 5.3298332
# RNA 4.409586

#checking for overall differences in saprotrophic OTU richness between rDNA and rRNA without outliers
env_nooutliers <- env[c(1:17, 19:36, 38:42),]
meanRichness <- ddply(env_nooutliers, .(nucleic_acid), summarize, mean=mean(richnessSapro))
# DNA 13.13636
# RNA 14.55556
sdRichness <- ddply(env_nooutliers, .(nucleic_acid), summarize, sd=sd(richnessSapro))
# DNA 3.681297
# RNA 4.076074


#checking for overall differences in pathotrophic OTU richness between rDNA and rRNA 
otus_patho <- subset(otus_tax, otus_tax$Simplified_Trophic_mode=="Pathotroph") #subsetting pathotrophs in OTU table
otus_patho <- otus_patho[,1:42]
otus_patho <- t(otus_patho)
env$richnessPatho<-specnumber(otus_patho)
meanRichness <- ddply(env, .(nucleic_acid), summarize, mean=mean(richnessPatho))
# DNA 5.5652172
# RNA 4.947368
sdRichness <- ddply(env, .(nucleic_acid), summarize, sd=sd(richnessPatho))
# DNA 2.9049502
# RNA 2.296705

#checking for overall differences in pathotrophic OTU richness between rDNA and rRNA without outliers
env_nooutliers <- env[c(1:17, 19:36, 38:42),]
meanRichness <- ddply(env_nooutliers, .(nucleic_acid), summarize, mean=mean(richnessPatho))
# DNA 5.136364
# RNA 4.888889
sdRichness <- ddply(env_nooutliers, .(nucleic_acid), summarize, sd=sd(richnessPatho))
# DNA 2.099784
# RNA 2.348689

#checking for overall differences in unassigned OTU richness between rDNA and rRNA 
otus_unass <- subset(otus_tax, otus_tax$Simplified_Trophic_mode=="unassigned") #subsetting pathotrophs in OTU table
otus_unass <- otus_unass[,1:42]
otus_unass <- t(otus_unass)
env$richnessUnass<-specnumber(otus_unass)
meanRichness <- ddply(env, .(nucleic_acid), summarize, mean=mean(richnessUnass))
# DNA 46.434782
# RNA 50.94737
sdRichness <- ddply(env, .(nucleic_acid), summarize, sd=sd(richnessUnass))
# DNA 14.761212
# RNA 10.15696

#checking for overall differences in unassigned OTU richness between rDNA and rRNA without outliers
env_nooutliers <- env[c(1:17, 19:36, 38:42),]
meanRichness <- ddply(env_nooutliers, .(nucleic_acid), summarize, mean=mean(richnessUnass))
# DNA 44.40909
# RNA 49.83333
sdRichness <- ddply(env_nooutliers, .(nucleic_acid), summarize, sd=sd(richnessUnass))
# DNA 11.375241
# RNA  9.179581

#environmental variables vs richness
plot(richness ~ pH, data = env)
plot(richness ~ moisture, data = env)
plot(richness ~ conductivity, data = env)
plot(richness ~ OM, data = env)
plot(richness ~ N, data = env)
plot(richness ~ C, data = env)
plot(richness ~ cn_ratio, data = env)
plot(richness ~ nucleic_acid, data = env)
plot(richness ~ vegetation_type, data = env)
plot(richness ~ snow_regime, data = env)

#overall differences between mean richness acording to single factor
library(ggplot2)
library(dplyr)
library(ggpubr)

overal1 <- ggplot(env, aes(x=nucleic_acid, y=richness))+
  geom_boxplot()+
  stat_boxplot(geom ='errorbar')+ 
  geom_jitter(width = 0.2)+
  stat_summary(fun.y=mean,geom="line",lwd=2,aes(group=1))+
  theme_bw()+
  labs(x = " ", y="Richness of all OTUs")

overal2 <- ggplot(env, aes(x=vegetation_type, y=richness))+
  geom_boxplot()+
  stat_boxplot(geom ='errorbar')+
  geom_jitter(width = 0.2)+
  stat_summary(fun.y=mean,geom="line",lwd=2,aes(group=1))+
  theme_bw()+
  labs(x = " ", y=" ")

overal3 <- ggplot(env, aes(x=snow_regime, y=richness))+
  geom_boxplot()+
  stat_boxplot(geom ='errorbar')+
  geom_jitter(width = 0.2)+
  stat_summary(fun.y=mean,geom="line",lwd=2,aes(group=1))+
  theme_bw()+
  labs(x = " ", y=" ")

richness_factor <- ggarrange(overal1, overal2, overal3 + font("x.text", size = 10),
                    ncol = 3, nrow = 1)
annotate_figure(richness_factor, bottom = text_grob("Experimental factor", color = "black", size = 10)) #fig.lab = "Figure 1", fig.lab.face = "bold")
richness_factor


meanRichness <- ddply(env, .(vegetation_type), summarize, mean=mean(richness))
#heath 106.8182
#meadow 107.5000

sdRichness <- ddply(env, .(vegetation_type), summarize, sd=sd(richness))
# heath 17.41274
# meadow 29.05983

meanRichness <- ddply(env_nooutliers, .(vegetation_type), summarize, mean=mean(richness)) #nooutliers
#heath 107.8182
#meadow 101.8889

sdRichness <- ddply(env_nooutliers, .(vegetation_type), summarize, sd=sd(richness)) #nooutliers
#heath 17.41274
#meadow 21.71443

meanRichness <- ddply(env, .(snow_regime), summarize, mean=mean(richness))
# ambient 113.35
# deep 101.50

sdRichness <- ddply(env, .(snow_regime), summarize, sd=sd(richness))
# ambient 24.52770
# deep 21.28659

meanRichness <- ddply(env_nooutliers, .(snow_regime), summarize, mean=mean(richness)) #nooutliers
# ambient 108.3889
# deep 102.5000

sdRichness <- ddply(env_nooutliers, .(snow_regime), summarize, sd=sd(richness)) #nooutliers
# ambient 16.92940
# deep 21.28659


#testing if experimental factors show any effect of overall OTU richness
mod_rich_lmer <- lmer(richness ~ nucleic_acid * vegetation_type + (1|block/fence), data=env) #almost 0 random effects for block/fence
summary(mod_rich_lmer)
Anova(mod_rich_lmer, test.statistic = "F", type = "II")

mod_rich_lm_int <- lm(richness ~ nucleic_acid * vegetation_type, data=env)
summary(mod_rich_lm_int)

#testing if experimental factors show any effect of overall OTU richness (NO OUTLIERS)
mod_rich_lmer <- lmer(richness ~ nucleic_acid * vegetation_type + (1|block/fence), data=env_nooutliers) #almost 0 random effects for block/fence
summary(mod_rich_lmer)
Anova(mod_rich_lmer, test.statistic = "F", type = "II")


#testing if experimental factors show any effect of symbiotrophic OTU richness
mod_richSymbio <- lmer(richnessSymbio ~ nucleic_acid * vegetation_type + (1|block/fence), data=env)
summary(mod_richSymbio)
Anova(mod_richSymbio, test.statistic = "F", type = "II")

#testing if experimental factors show any effect of symbiotrophic OTU richness (NO OUTLIERS)
mod_richSymbio <- lmer(richnessSymbio ~ nucleic_acid * vegetation_type + (1|block/fence), data=env_nooutliers)
summary(mod_richSymbio)
Anova(mod_richSymbio, test.statistic = "F", type = "II")


#testing if experimental factors show any effect of saprotrophic OTU richness, and if there is an interaction
mod_richSapro_int <- lmer(richnessSapro ~ nucleic_acid * vegetation_type +  (1|block/fence), data=env)
summary(mod_richSapro_int) #random variability was 0 for both block and fence, so:

mod_richSapro_int <- lm(richnessSapro ~ nucleic_acid * vegetation_type, data=env)
summary(mod_richSapro_int)
mod_richSapro_add <- lm(richnessSapro ~ nucleic_acid + vegetation_type, data=env)
summary(mod_richSapro_add)

anova(mod_richSapro_int, mod_richSapro_add) #no statistical difference between additive linear model and lm with interaction, showing the one with interaction

#testing if experimental factors show any effect of saprotrophic OTU richness, and if there is an interaction (NO OUTLIERS)
mod_richSapro_int <- lmer(richnessSapro ~ nucleic_acid * vegetation_type +  (1|block/fence), data=env_nooutliers)
summary(mod_richSapro_int) #random variability was 0 for both block and fence, so:

mod_richSapro_int <- lm(richnessSapro ~ nucleic_acid * vegetation_type, data=env_nooutliers)
summary(mod_richSapro_int)
mod_richSapro_add <- lm(richnessSapro ~ nucleic_acid + vegetation_type, data=env_nooutliers)
summary(mod_richSapro_add)

anova(mod_richSapro_int, mod_richSapro_add) #no statistical difference between additive linear model and lm with interaction, showing the one with interaction

#testing if experimental factors show any effect of pathotrophic OTU richness, and if there is an interaction
mod_richPatho <- lmer(richnessPatho ~ nucleic_acid * vegetation_type +  (1|block/fence), data=env)
summary(mod_richPatho)
Anova(mod_richPatho, test.statistic = "F", type = "II")

#testing if experimental factors show any effect of pathotrophic OTU richness, and if there is an interaction
mod_richPatho <- lmer(richnessPatho ~ nucleic_acid * vegetation_type +  (1|block/fence), data=env_nooutliers)
summary(mod_richPatho)
Anova(mod_richPatho, test.statistic = "F", type = "II")

mod_richPatho_int <- lm(richnessPatho ~ nucleic_acid * vegetation_type, data=env_nooutliers)
summary(mod_richPatho_int)
mod_richPatho_add <- lm(richnessPatho ~ nucleic_acid + vegetation_type, data=env_nooutliers)
summary(mod_richPatho_add)

anova(mod_richSapro_int, mod_richSapro_add) #no statistical differences between int and add models


#testing if experimental factors show any effect of unassigned OTU richness, and if there is an interaction
mod_richUnassigned <- lmer(richnessUnass ~ nucleic_acid * vegetation_type +  (1|block/fence), data=env)
summary(mod_richUnassigned) #random effects ~0, thus linear model:

mod_richUnass_int <- lm(richnessUnass ~ nucleic_acid * vegetation_type, data=env)
summary(mod_richUnass_int)
mod_richUnass_add <- lm(richnessUnass ~ nucleic_acid + vegetation_type, data=env)
summary(mod_richUnass_add)

anova(mod_richUnass_int, mod_richUnass_add) #no statistical difference between additive linear model and lm with interaction, showing the one with interaction

#testing if experimental factors show any effect of unassigned OTU richness, and if there is an interaction (NO OUTLIERS)
mod_richUnassigned <- lmer(richnessUnass ~ nucleic_acid * vegetation_type +  (1|block/fence), data=env_nooutliers)
summary(mod_richUnassigned) #random effects ~0, thus linear model:
Anova(mod_richUnassigned, test.statistic = "F", type = "II")


#testing if 
mod_rich_lmer <- lmer(richness ~ nucleic_acid + snow_regime + (1|block/fence), data=env) #almost 0 random effects for block/fence
summary(mod_rich_lmer)
Anova(mod_rich_lmer)
mod_rich_lm_int <- lm(richness ~ nucleic_acid * vegetation_type, data=env)
summary(mod_rich_lm_int)
mod_rich_lm_add <- lm(richness ~ nucleic_acid + vegetation_type, data=env)
summary(mod_rich_lm_add)


mod_rich_DNARNA <- lm(richness ~ nucleic_acid, data = env)
summary(mod_rich_DNARNA)
mod_rich_snowdepth <- lm(richness ~ snow_regime, data=env)
summary(mod_rich_snowdepth)
mod_rich_veg <- lm(richness ~ vegetation_type, data=env)
summary(mod_rich_veg)

#mean richness
meanRichness <- ddply(env, .(lumped), summarize, mean=mean(richness))
# DNA_H 100.0909
# DNA_M 104.0833
# RNA_H 115.5455
# RNA_M 115.1250
sdRichness <- ddply(env, .(lumped), summarize, sd=sd(richness))
# DNA_H 17.05552
# DNA_M 30.23982
# RNA_H 14.64489
# RNA_M 27.77685

meanRichness <- ddply(env_nooutliers, .(lumped), summarize, mean=mean(richness)) #NO OUTLIERS
# DNA_H 100.09091
# DNA_M  97.45455
# RNA_H 115.54545
# RNA_M 108.85714
sdRichness <- ddply(env_nooutliers, .(lumped), summarize, sd=sd(richness)) #NO OUTLIERS
# DNA_H 17.05552
# DNA_M 20.63668
# RNA_H 14.64489
# RNA_M 23.09710

meanRichness <- ddply(env, .(lumped), summarize, mean=mean(richnessSymbio))
# DNA_H 39.36364
# DNA_M 31.41667
# RNA_H 45.09091
# RNA_M 41.25000
sdRichness <- ddply(env, .(lumped), summarize, sd=sd(richnessSymbio))
# DNA_H 9.821128
# DNA_M 6.625822
# RNA_H 10.222079
# RNA_M 9.924717

meanRichness <- ddply(env_nooutliers, .(lumped), summarize, mean=mean(richnessSymbio)) #NO OUTLIERS
# DNA_H 39.36364
# DNA_M 30.81818
# RNA_H 45.09091
# RNA_M 38.85714
sdRichness <- ddply(env_nooutliers, .(lumped), summarize, sd=sd(richnessSymbio))
#  DNA_H  9.821128
#  DNA_M  6.600275
#  RNA_H 10.222079
#  RNA_M  7.840675

meanRichness <- ddply(env, .(lumped), summarize, mean=mean(richnessSapro))
# DNA_H 12.63636
# DNA_M 15.16667
# RNA_H 13.63636
# RNA_M 16.87500
sdRichness <- ddply(env, .(lumped), summarize, sd=sd(richnessSapro))
# DNA_H 2.500909
# DNA_M 6.912878
# RNA_H 3.009077
# RNA_M 5.488625

meanRichness <- ddply(env_nooutliers, .(lumped), summarize, mean=mean(richnessSapro)) #NO OUTLIERS
# DNA_H 12.63636
# DNA_M 13.63636
# RNA_H 13.63636
# RNA_M 16.00000
sdRichness <- ddply(env_nooutliers, .(lumped), summarize, sd=sd(richnessSapro)) #NO OUTLIERS
# DNA_H 2.500909
# DNA_M 4.653444
# RNA_H 3.009077
# RNA_M 5.291503

meanRichness <- ddply(env, .(lumped), summarize, mean=mean(richnessPatho))
# DNA_H 4.818182
# DNA_M 6.250000
# RNA_H 5.454545
# RNA_M 4.250000
sdRichness <- ddply(env, .(lumped), summarize, sd=sd(richnessPatho))
# DNA_H 1.250454
# DNA_M 3.792936
# RNA_H 2.381749
# RNA_M 2.121320

meanRichness <- ddply(env_nooutliers, .(lumped), summarize, mean=mean(richnessPatho)) #NO OUTLIERS
# DNA_H 4.818182
# DNA_M 5.454545
# RNA_H 5.454545
# RNA_M 4.000000
sdRichness <- ddply(env_nooutliers, .(lumped), summarize, sd=sd(richnessPatho)) #NO OUTLIERS
# DNA_H 1.250454
# DNA_M 2.733629
# RNA_H 2.381749
# RNA_M 2.160247

meanRichness <- ddply(env, .(lumped), summarize, mean=mean(richnessUnass))
# DNA_H 42.27273
# DNA_M 50.25000
# RNA_H 50.36364
# RNA_M 51.75000
sdRichness <- ddply(env, .(lumped), summarize, sd=sd(richnessUnass))
# DNA_H 9.371136
# DNA_M 17.965244
# RNA_H 7.486958
# RNA_M 13.562027

meanRichness <- ddply(env_nooutliers, .(lumped), summarize, mean=mean(richnessUnass)) #NO OUTLIERS
# DNA_H 42.27273
# DNA_M 46.54545
# RNA_H 50.36364
# RNA_M 49.00000
sdRichness <- ddply(env_nooutliers, .(lumped), summarize, sd=sd(richnessUnass)) #NO OUTLIERS
# 1  DNA_H  9.371136
# 2  DNA_M 13.186081
# 3  RNA_H  7.486958
# 4  RNA_M 12.000000

table(env_nooutliers$lumped)


library(ggplot2)
library(dplyr)
library(ggpubr)

#plotting overall richness divided into diff vegetation types / nucleic acid
p1 <- ggplot(env, aes(x=nucleic_acid, y=richness))+
  ylim(0,180) +
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  facet_wrap(~vegetation_type)+
  theme_bw()+
  labs(x = " ", y="Richness of all OTUs")

#plotting symbiotrophic richness divided into diff vegetation types / nucleic acid
p2 <- ggplot(env, aes(x=nucleic_acid, y=richnessSymbio))+
  ylim(0,180) +
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  facet_wrap(~vegetation_type)+
  theme_bw()+
  labs(x = " ", y="Richness of symbiotrophs")

#plotting saprotrophic richness divided into diff vegetation types / nucleic acid
p3 <- ggplot(env, aes(x=nucleic_acid, y=richnessSapro))+
  ylim(0,180) +
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  facet_wrap(~vegetation_type)+
  theme_bw()+
  labs(x = " ", y="Richness of saprotrophs")

#plotting pathotrophic richness divided into diff vegetation types / nucleic acid
p4 <- ggplot(env, aes(x=nucleic_acid, y=richnessPatho))+
  ylim(0,180) +
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  facet_wrap(~vegetation_type)+
  theme_bw()+
  labs(x = " ", y="Richness of pathotrophs")

#plotting pathotrophic richness divided into diff vegetation types / nucleic acid
p5 <- ggplot(env, aes(x=nucleic_acid, y=richnessUnass))+
  ylim(0,180) +
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  facet_wrap(~vegetation_type)+
  theme_bw()+
  labs(x = " ", y="Richness of unassigned reads")

p_all <- ggarrange(p1, p2, p3, p4, p5 + font("x.text", size = 10), ncol = 5, nrow = 1)
annotate_figure(p_all, bottom = text_grob("Experimental factor", color = "black", size = 10)) #fig.lab = "Figure 1", fig.lab.face = "bold")

#plotting the same figure without outliers
p1_no <- ggplot(env_nooutliers, aes(x=nucleic_acid, y=richness))+
  ylim(0,150) +
  geom_boxplot()+
  stat_summary(fun.y=mean,geom="line",lwd=1,aes(group=1, colour="red"))+
  geom_jitter(width = 0.2)+
  facet_wrap(~vegetation_type)+
  stat_summary(fun.y=mean,geom="line",lwd=1,aes(group=1, colour="red"))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = " ", y="Richness of all OTUs")

#plotting symbiotrophic richness divided into diff vegetation types / nucleic acid
p2_no <- ggplot(env_nooutliers, aes(x=nucleic_acid, y=richnessSymbio))+
  ylim(0,150) +
  geom_boxplot()+
  stat_summary(fun.y=mean,geom="line",lwd=1,aes(group=1, colour="red"))+
  geom_jitter(width = 0.2)+
  facet_wrap(~vegetation_type)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = " ", y="Richness of symbiotrophs")

#plotting saprotrophic richness divided into diff vegetation types / nucleic acid
p3_no <- ggplot(env_nooutliers, aes(x=nucleic_acid, y=richnessSapro))+
  ylim(0,150) +
  geom_boxplot()+
  stat_summary(fun.y=mean,geom="line",lwd=1,aes(group=1, colour="red"))+
  geom_jitter(width = 0.2)+
  facet_wrap(~vegetation_type)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = " ", y="Richness of saprotrophs")

#plotting pathotrophic richness divided into diff vegetation types / nucleic acid
p4_no <- ggplot(env_nooutliers, aes(x=nucleic_acid, y=richnessPatho))+
  ylim(0,150) +
  geom_boxplot()+
  stat_summary(fun.y=mean,geom="line",lwd=1,aes(group=1, colour="red"))+
  geom_jitter(width = 0.2)+
  facet_wrap(~vegetation_type)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = " ", y="Richness of pathotrophs")

#plotting pathotrophic richness divided into diff vegetation types / nucleic acid
p5_no <- ggplot(env_nooutliers, aes(x=nucleic_acid, y=richnessUnass))+
  ylim(0,150) +
  geom_boxplot()+
  stat_summary(fun.y=mean,geom="line",lwd=1,aes(group=1, colour="red"))+
  geom_jitter(width = 0.2)+
  facet_wrap(~vegetation_type)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = " ", y="Richness of unassigned reads")

p_no_all <- ggarrange(p1_no, p2_no, p3_no, p4_no, p5_no + font("x.text", size = 10), ncol = 5, nrow = 1)
annotate_figure(p_no_all, bottom = text_grob("Experimental factor", color = "black", size = 10)) #fig.lab = "Figure 1", fig.lab.face = "bold")


#test for differences in means
env_heath <- subset(env, vegetation_type=="heath")
lm_heath <- lm(richness ~ nucleic_acid, data=env_heath)
summary(lm_heath)
# Call:
#   lm(formula = richness ~ nucleic_acid, data = env_heath)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -37.091  -8.295  -1.318   7.909  25.455 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       99.091      4.793   20.68 5.73e-15 ***
#   nucleic_acidRNA   15.455      6.778    2.28   0.0337 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 15.9 on 20 degrees of freedom
# Multiple R-squared:  0.2063,	Adjusted R-squared:  0.1666 
# F-statistic: 5.199 on 1 and 20 DF,  p-value: 0.03371

env_meadow <- subset(env, vegetation_type=="meadow")
lm_meadow <- lm(richness ~ nucleic_acid, data=env_meadow)
summary(lm_meadow)
# Call:
#   lm(formula = richness ~ nucleic_acid, data = env_meadow)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -42.125 -19.594  -6.583  16.917  72.917 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       103.08       8.46  12.185 3.95e-10 ***
#   nucleic_acidRNA    11.04      13.38   0.825     0.42    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 29.31 on 18 degrees of freedom
# Multiple R-squared:  0.03647,	Adjusted R-squared:  -0.01706 
# F-statistic: 0.6814 on 1 and 18 DF,  p-value: 0.4199


env_heath <- subset(env, vegetation_type=="heath")
lm_heath <- lm(richnessSymbio ~ nucleic_acid, data=env_heath)
summary(lm_heath)

env_meadow <- subset(env, vegetation_type=="meadow")
lm_meadow <- lm(richnessSymbio ~ nucleic_acid, data=env_meadow)
summary(lm_meadow)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -11.4167  -5.6667   0.1667   6.5833  16.7500 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       31.417      2.330  13.485 7.55e-11 ***
#   nucleic_acidRNA    9.833      3.684   2.669   0.0156 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 8.071 on 18 degrees of freedom
# Multiple R-squared:  0.2836,	Adjusted R-squared:  0.2438 
# F-statistic: 7.126 on 1 and 18 DF,  p-value: 0.01563


env_heath <- subset(env, vegetation_type=="heath")
lm_heath <- lm(richnessSapro ~ nucleic_acid, data=env_heath)
summary(lm_heath)

env_meadow <- subset(env, vegetation_type=="meadow")
lm_meadow <- lm(richnessSapro ~ nucleic_acid, data=env_meadow)
summary(lm_meadow)

env_heath <- subset(env, vegetation_type=="heath")
lm_heath <- lm(richnessPatho ~ nucleic_acid, data=env_heath)
summary(lm_heath)

env_meadow <- subset(env, vegetation_type=="meadow")
lm_meadow <- lm(richnessPatho ~ nucleic_acid, data=env_meadow)
summary(lm_meadow)

env_heath <- subset(env, vegetation_type=="heath")
lm_heath <- lm(richnessUnass ~ nucleic_acid, data=env_heath)
summary(lm_heath)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -19.273  -7.841   3.182   5.409  12.727 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       42.273      2.557  16.530 3.96e-13 ***
#   nucleic_acidRNA    8.091      3.617   2.237   0.0368 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 8.482 on 20 degrees of freedom
# Multiple R-squared:  0.2002,	Adjusted R-squared:  0.1602 
# F-statistic: 5.005 on 1 and 20 DF,  p-value: 0.03682


env_meadow <- subset(env, vegetation_type=="meadow")
lm_meadow <- lm(richnessUnass ~ nucleic_acid, data=env_meadow)
summary(lm_meadow)


#THE SAME STUFF BUT FOR SNOW REGIME:
  #plotting overall richness divided into diff vegetation types / nucleic acid
  s1 <- ggplot(env, aes(x=nucleic_acid, y=richness))+
  ylim(0,180) +
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  facet_wrap(~snow_regime)+
  theme_bw()+
  labs(x = " ", y="Richness of all OTUs")

#plotting symbiotrophic richness divided into diff vegetation types / nucleic acid
s2 <- ggplot(env, aes(x=nucleic_acid, y=richnessSymbio))+
  ylim(0,180) +
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  facet_wrap(~snow_regime)+
  theme_bw()+
  labs(x = " ", y="Richness of symbiotrophs")

#plotting saprotrophic richness divided into diff vegetation types / nucleic acid
s3 <- ggplot(env, aes(x=nucleic_acid, y=richnessSapro))+
  ylim(0,180) +
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  facet_wrap(~snow_regime)+
  theme_bw()+
  labs(x = " ", y="Richness of saprotrophs")

#plotting pathotrophic richness divided into diff vegetation types / nucleic acid
s4 <- ggplot(env, aes(x=nucleic_acid, y=richnessPatho))+
  ylim(0,180) +
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  facet_wrap(~snow_regime)+
  theme_bw()+
  labs(x = " ", y="Richness of pathotrophs")

#plotting pathotrophic richness divided into diff vegetation types / nucleic acid
s5 <- ggplot(env, aes(x=nucleic_acid, y=richnessUnass))+
  ylim(0,180) +
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  facet_wrap(~snow_regime)+
  theme_bw()+
  labs(x = " ", y="Richness of unassigned reads")

s_all <- ggarrange(s1, s2, s3, s4, s5 + font("x.text", size = 10), ncol = 5, nrow = 1)
annotate_figure(s_all, bottom = text_grob("Experimental factor", color = "black", size = 10)) #fig.lab = "Figure 1", fig.lab.face = "bold")

#test for differences in means 
env_ambient <- subset(env, snow_regime=="ambient")
lm_ambient <- lm(richness ~ nucleic_acid, data=env_ambient)
summary(lm_ambient)

env_deep <- subset(env, snow_regime=="deep")
lm_deep <- lm(richness ~ nucleic_acid, data=env_deep)
summary(lm_deep)

lm_ambient <- lm(richnessSymbio ~ nucleic_acid, data=env_ambient)
summary(lm_ambient)
# Call:
#   lm(formula = richnessSymbio ~ nucleic_acid, data = env_ambient)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -16.200  -4.675  -0.100   3.625  20.800 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       36.200      3.261   11.10 1.75e-09 ***
#   nucleic_acidRNA   11.900      4.612    2.58   0.0189 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 10.31 on 18 degrees of freedom
# Multiple R-squared:   0.27,	Adjusted R-squared:  0.2294 
# F-statistic: 6.657 on 1 and 18 DF,  p-value: 0.01887

lm_deep <- lm(richnessSymbio ~ nucleic_acid, data=env_deep)
summary(lm_deep)

lm_ambient <- lm(richnessSapro ~ nucleic_acid, data=env_heath)
summary(lm_heath)

lm_deep <- lm(richnessSapro ~ nucleic_acid, data=env_deep)
summary(lm_deep)

lm_ambient <- lm(richnessPatho ~ nucleic_acid, data=env_ambient)
summary(lm_ambient)

lm_deep <- lm(richnessPatho ~ nucleic_acid, data=env_deep)
summary(lm_deep)

lm_ambient <- lm(richnessUnass ~ nucleic_acid, data=env_ambient)
summary(lm_ambient)

lm_deep <- lm(richnessUnass ~ nucleic_acid, data=env_deep)
summary(lm_deep)

table(env$treatment) #coding for combination of all experimental factors (nucleic acid, vegetation type and treatment)
# DNA_H_A DNA_H_D DNA_M_A DNA_M_D RNA_H_A RNA_H_D RNA_M_A RNA_M_D 
# 3       8       7       5       5       6       5       3
meanRichness <- ddply(env, .(treatment), summarize, mean=mean(richness))
# 1   DNA_H_A 105.3333
# 2   DNA_H_D  98.1250
# 3   DNA_M_A 110.4286
# 4   DNA_M_D  95.2000
# 5   RNA_H_A 117.2000
# 6   RNA_H_D 114.1667
# 7   RNA_M_A 122.4000
# 8   RNA_M_D 103.0000


sdRichness <- ddply(env, .(treatment), summarize, sd=sd(richness))
# 1   DNA_H_A 17.039171
# 2   DNA_H_D 17.787937
# 3   DNA_M_A 31.884464
# 4   DNA_M_D 28.647862
# 5   RNA_H_A 21.452273
# 6   RNA_H_D  7.467708
# 7   RNA_M_A 23.415807
# 8   RNA_M_D 35.369478

#plotting overall richness divided into diff vegetation types / nucleic acid
x1 <- ggplot(env, aes(x=nucleic_acid, y=richness))+
  geom_point(shape=21)+
  facet_grid(vegetation_type~snow_regime)+
  stat_summary(fun.y=mean,geom="line",lwd=1,aes(group=1, colour="red"))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = " ", y="Richness of all OTUs")

#plotting symbiotrophic richness divided into diff vegetation types / nucleic acid
x2 <- ggplot(env, aes(x=nucleic_acid, y=richnessSymbio)) +
  geom_point(shape=21)+
  facet_grid(vegetation_type~snow_regime)+
  stat_summary(fun.y=mean,geom="line",lwd=1,aes(group=1, colour="red"))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = " ", y="Richness of symbiotrophs")

#plotting saprotrophic richness divided into diff vegetation types / nucleic acid
x3 <- ggplot(env, aes(x=nucleic_acid, y=richnessSapro))+
  geom_point(shape=21)+
  facet_grid(vegetation_type~snow_regime)+
  stat_summary(fun.y=mean,geom="line",lwd=1,aes(group=1, colour="red"))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = " ", y="Richness of saprotrophs")

#plotting pathotrophic richness divided into diff vegetation types / nucleic acid
x4 <- ggplot(env, aes(x=nucleic_acid, y=richnessPatho))+
  geom_point(shape=21)+
  facet_grid(vegetation_type ~ snow_regime)+
  stat_summary(fun.y=mean,geom="line",lwd=1,aes(group=1, colour="red"))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = " ", y="Richness of pathotrophs")

#plotting pathotrophic richness divided into diff vegetation types / nucleic acid
x5 <- ggplot(env, aes(x=nucleic_acid, y=richnessUnass))+
  geom_point(shape=21)+
  facet_grid(vegetation_type~snow_regime)+
  stat_summary(fun.y=mean,geom="line",lwd=1,aes(group=1, colour="red"))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = " ", y="Richness of unassigned reads")

x_all <- ggarrange(x1, x2, x3, x4, x5 + font('x.text', size = 10), ncol = 5, nrow = 1)
x_all
annotate_figure(s_all, bottom = text_grob("Experimental factor", color = "black", size = 10)) #fig.lab = "Figure 1", fig.lab.face = "bold")


#comparison of DNA and RNA samples 
richness_DNA <- as.numeric(specnumber(t(comparison_DNA)))
richness_RNA <- as.numeric(specnumber(t(comparison_RNA)))
richness_comparison <- rbind(richness_DNA, richness_RNA)
richness_comparison <- data.frame(richness_comparison)
richness_comparison$nucleic_acid <- c("DNA", "RNA")
library(reshape2)
richness_comparison_long <- melt(richness_comparison, v.names = "OTU_number", idvar = c("nucleic_acid"))
lm(value ~ nucleic_acid, data=richness_comparison_long)
summary(lm(value ~ nucleic_acid, data=richness_comparison_long))
plot(richness_DNA ~ richness_RNA, data = richness_comparison)


##################################
###   TAXONOMY - SAPROTROPHS   ###
##################################

saprofortax <- read.delim("saprofortax.txt", sep="\t", header=TRUE)
class(saprofortax)
class(saprofortax$Order)

library(ggplot2)
ggplot(saprofortax, aes(fill=saprofortax$Phylum, x=saprofortax$combined, y=saprofortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(saprofortax, aes(fill=saprofortax$Class, x=saprofortax$combined, y=saprofortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(saprofortax, aes(fill=saprofortax$Class, x=saprofortax$nucleic_acid, y=saprofortax$abundance)) +
  geom_bar(stat = "identity", position = "fill")

ggplot(saprofortax, aes(fill=saprofortax$nucleic_acid, x=saprofortax$Order, y=saprofortax$abundance)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

########################
ggplot(saprofortax, aes(fill=saprofortax$nucleic_acid, x=saprofortax$Phylum, y=saprofortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #saprofortax_phylum_dnarna.pdf

ggplot(saprofortax, aes(fill=saprofortax$nucleic_acid, x=saprofortax$Class, y=saprofortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #saprofortax_class_dnarna.pdf

ggplot(saprofortax, aes(fill=saprofortax$nucleic_acid, x=saprofortax$Order, y=saprofortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #saprofortax_order_dnarna.pdf

ggplot(saprofortax, aes(fill=saprofortax$vegetation, x=saprofortax$Order, y=saprofortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #saprofortax_order_vegetation.pdf

ggplot(saprofortax, aes(fill=saprofortax$nucleic_acid, x=saprofortax$Order, y=saprofortax$abundance)) +
  geom_bar(stat = "identity") +
  facet_wrap(~saprofortax$combined) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #saprofortax_order_combined_abundance.pdf
  
ggplot(saprofortax, aes(fill=saprofortax$nucleic_acid, x=saprofortax$Order, y=saprofortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #saprofortax_order_combined_abundance.pdf

ggplot(saprofortax, aes(fill=saprofortax$nucleic_acid, x=saprofortax$Order, y=saprofortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~saprofortax$combined) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #saprofortax_order_combined_presence.pdf

ggplot(saprofortax, aes(fill=saprofortax$vegetation, x=saprofortax$Order, y=saprofortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~saprofortax$nucleic_acid) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #saprofortax_order_dnarna_relativeabundance_meadow.pdf

ggplot(saprofortax, aes(fill=saprofortax$vegetation, x=saprofortax$Base_for_functional_assignment, y=saprofortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~saprofortax$nucleic_acid, nrow=2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #saprofortax_baseforfunct_dnarna_vegetation.pdf


###############################################
###   TAXONOMY - FUNCTIONALLY UNANNOTATED   ###
###############################################

unassfortax <- read.delim("unassfortax.txt", sep="\t", header=TRUE)
class(unassfortax)
class(unassfortax$Order)

ggplot(unassfortax, aes(fill=unassfortax$Phylum, x=unassfortax$combined, y=unassfortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #unassfortax_combined_phylum_abundance.pdf

ggplot(unassfortax, aes(fill=unassfortax$Class, x=unassfortax$combined, y=unassfortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #unassfortax_combined_class_abundance.pdf

ggplot(unassfortax, aes(fill=unassfortax$Order, x=unassfortax$combined, y=unassfortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #unassfortax_combined_order_abundance.pdf

ggplot(unassfortax, aes(fill=unassfortax$Class, x=unassfortax$nucleic_acid, y=unassfortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") #unass_dnarna_class_abundance.pdf

ggplot(unassfortax, aes(fill=unassfortax$nucleic_acid, x=unassfortax$Order, y=unassfortax$abundance)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

########################
ggplot(unassfortax, aes(fill=unassfortax$nucleic_acid, x=unassfortax$Phylum, y=unassfortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #unassfortax_phylum_dnarna.pdf

ggplot(unassfortax, aes(fill=unassfortax$nucleic_acid, x=unassfortax$Class, y=unassfortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #unassfortax_class_dnarna.pdf

ggplot(unassfortax, aes(fill=unassfortax$nucleic_acid, x=unassfortax$Order, y=unassfortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #unassfortax_order_dnarna.pdf

ggplot(unassfortax, aes(fill=unassfortax$vegetation, x=unassfortax$Order, y=unassfortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #unassfortax_order_vegetation.pdf

ggplot(unassfortax, aes(fill=unassfortax$nucleic_acid, x=unassfortax$Order, y=unassfortax$abundance)) +
  geom_bar(stat = "identity") +
  facet_wrap(~unassfortax$combined) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #unassfortax_order_combined_abundance.pdf

ggplot(unassfortax, aes(fill=unassfortax$nucleic_acid, x=unassfortax$Order, y=unassfortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~unassfortax$combined) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #unassfortax_order_combined_presence.pdf

ggplot(unassfortax, aes(fill=unassfortax$vegetation, x=unassfortax$Order, y=unassfortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~unassfortax$nucleic_acid) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #unassfortax_order_dnarna_relativeabundance_vegetation.pdf

ggplot(unassfortax, aes(fill=vegetation, x=unassfortax$Base_for_functional_assignment, y=unassfortax$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~unassfortax$nucleic_acid, nrow=2) +
  theme(text = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1)) #unassfortax_baseforfunct_dnarna_vegetation.pdf


##################################
###     TAXONOMY - 9SAMPLES    ###
##################################

library(ggplot2)

allOTUs_9samples <- read.delim("9samples_allOTUS.txt", sep="\t", header=TRUE)
class(allOTUs_9samples)
class(allOTUs_9samples$Order)

allOTUs_9samples_symbio <- allOTUs_9samples[which(allOTUs_9samples$Simplified_Trophic_mode=="Symbiotroph"), ]
allOTUs_9samples_sapro <- allOTUs_9samples[which(allOTUs_9samples$Simplified_Trophic_mode=="Saprotroph"),]
allOTUs_9samples_patho <- allOTUs_9samples[which(allOTUs_9samples$Simplified_Trophic_mode=="Pathotroph"),]
allOTUs_9samples_unassigned <- allOTUs_9samples[which(allOTUs_9samples$Simplified_Trophic_mode=="unassigned"),]



ggplot(allOTUs_9samples, aes(fill=allOTUs_9samples$Phylum, x=allOTUs_9samples$vegetation, y=allOTUs_9samples$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~allOTUs_9samples$nucleic_acid) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #9samples_all_phylum.pdf

ggplot(allOTUs_9samples, aes(fill=allOTUs_9samples$Class, x=allOTUs_9samples$vegetation, y=allOTUs_9samples$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~allOTUs_9samples$nucleic_acid) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #9samples_all_class.pdf

ggplot(allOTUs_9samples, aes(fill=allOTUs_9samples$Order, x=allOTUs_9samples$vegetation, y=allOTUs_9samples$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~allOTUs_9samples$nucleic_acid) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #9samples_all_order.pdf

table_tax <-with(allOTUs_9samples, tapply(abundance, list(Order, combined), sum))


library(reshape)
table_tax_long <- melt(table_tax, id.vars=colnames)
ggplot(table_tax_long, aes(x=table_tax_long$X2, y=table_tax_long$value, group=table_tax_long$X1, col=table_tax_long$X1)) +
  geom_point() +
  geom_line()

ggplot(table_tax_long, aes(log(value)+0.01, X1)) + 
  geom_point() + 
  facet_grid(~ X2)

#########################################################################################################################
################sapro in 9 coextracted samples
#########################################################################################################################
#allOTUs_9samples_sapro <- allOTUs_9samples[which(allOTUs_9samples$Simplified_Trophic_mode=="Saprotroph"),]

library(ggplot2)
library(ggrepel)
library(ggpubr)

ggplot(allOTUs_9samples_sapro, aes(fill=allOTUs_9samples_sapro$vegetation, x=allOTUs_9samples_sapro$Order, y=allOTUs_9samples_sapro$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~allOTUs_9samples_sapro$nucleic_acid, ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()

ggplot(allOTUs_9samples_sapro, aes(fill=allOTUs_9samples_sapro$vegetation, x=allOTUs_9samples_sapro$Base_for_functional_assignment, y=allOTUs_9samples_sapro$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~allOTUs_9samples_sapro$nucleic_acid, ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()

ggplot(allOTUs_9samples_sapro, aes(fill=allOTUs_9samples_sapro$vegetation, x=allOTUs_9samples_sapro$Order, y=log(allOTUs_9samples_sapro$abundance)+0.01)) +
  geom_bar(stat = "identity") +
  facet_wrap(~allOTUs_9samples_sapro$combined, nrow = 4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

table_tax_sapro <-with(allOTUs_9samples_sapro, tapply(abundance, list(Order, combined), sum))
table_tax_sapro
class(table_tax_sapro) #matrix
tax_sapro <- data.frame(table_tax_sapro) #now it is df

tax_sapro_heath <- tax_sapro[is.finite(tax_sapro$dna_heath) & is.finite(tax_sapro$rna_heath), ]
tax_sapro_heath <- subset(tax_sapro_heath, select = c(dna_heath, rna_heath))
tax_sapro_heath.no0 = tax_sapro_heath[ rowSums(tax_sapro_heath)!=0, ]
tax_sapro_meadow <- tax_sapro[is.finite(tax_sapro$dna_meadow) & is.finite(tax_sapro$rna_meadow), ]
tax_sapro_meadow <- subset(tax_sapro_meadow, select = c(dna_meadow, rna_meadow))
tax_sapro_meadow.no0 = tax_sapro_meadow[ rowSums(tax_sapro_meadow)!=0, ]

names_heath <- rownames(tax_sapro_heath.no0)
length(names_heath)
a_heath <- ggplot(tax_sapro_heath.no0, aes(x=log(dna_heath), y=log(rna_heath))) +
  geom_point(color = "blue", size = 3) +
  geom_abline(intercept = 0) +
  xlim(0,10) +
  ylim(0,10) +
  theme_bw() +
  geom_label_repel(aes(label = names_heath),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') #+
  #ggsave("plotheath.png", width = 5, height = 5)

names_meadow <- rownames(tax_sapro_meadow.no0)
length(names_meadow)  
a_meadow <- ggplot(tax_sapro_meadow.no0, aes(x=log(dna_meadow), y=log(rna_meadow))) +
      geom_point(color = "blue", size = 3) +
      geom_abline(intercept = 0) +
      xlim(0,10) +
      ylim(0,10) +
      theme_bw() +
      geom_label_repel(aes(label = names_meadow),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') #+
  #ggsave("plotmeadow.png", width = 5, height = 5)

ggarrange(a_heath, a_meadow, labels = c("heath", "meadow"), common.legend = FALSE)

cor.test(tax_sapro_heath.no0$dna_heath, tax_sapro_heath.no0$rna_heath) #0.5998511
cor.test(tax_sapro_meadow.no0$dna_meadow, tax_sapro_meadow.no0$rna_meadow) #0.9719105 


table_tax_sapro <-with(allOTUs_9samples_sapro, tapply(abundance, list(Order, combined), sum))
table_tax_sapro
class(table_tax_sapro)
table_tax_sapro.noO = table_tax_sapro[ rowSums(table_tax_sapro)!=0, ]
table_tax_sapro.compl <- table_tax_sapro.noO[complete.cases(table_tax_sapro.noO), ]

allOTUs_9samples_sapro_heath <- subset(allOTUs_9samples_sapro, allOTUs_9samples_sapro$vegetation=="heath") #subsetting
a <- ggplot(allOTUs_9samples_sapro_heath, aes(fill=allOTUs_9samples_sapro_heath$nucleic_acid, x=allOTUs_9samples_sapro_heath$Order, y=allOTUs_9samples_sapro_heath$abundance)) +
      geom_bar(stat = "identity", position = "fill") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

allOTUs_9samples_sapro_meadow <- subset(allOTUs_9samples_sapro, allOTUs_9samples_sapro$vegetation=="meadow") #subsetting
b <- ggplot(allOTUs_9samples_sapro_meadow, aes(fill=allOTUs_9samples_sapro_meadow$nucleic_acid, x=allOTUs_9samples_sapro_meadow$Order, y=allOTUs_9samples_sapro_meadow$abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggarrange(a, b, labels = c("heath", "meadow"), common.legend = TRUE, legend = "bottom")




#########################################################################################################################
################symbio in 9 coextracted samples
#########################################################################################################################

table_tax_symbio <-with(allOTUs_9samples_symbio, tapply(abundance, list(Order, combined), sum))
table_tax_symbio
class(table_tax_symbio)
tax_symbio <- data.frame(table_tax_symbio)

tax_symbio_heath <- tax_symbio[is.finite(tax_symbio$dna_heath) & is.finite(tax_symbio$rna_heath), ]
tax_symbio_heath <- subset(tax_symbio_heath, select = c(dna_heath, rna_heath))
tax_symbio_heath.no0 = tax_symbio_heath[ rowSums(tax_symbio_heath)!=0, ]
tax_symbio_meadow <- tax_symbio[is.finite(tax_symbio$dna_meadow) & is.finite(tax_symbio$rna_meadow), ]
tax_symbio_meadow <- subset(tax_symbio_meadow, select = c(dna_meadow, rna_meadow))
tax_symbio_meadow.no0 = tax_symbio_meadow[ rowSums(tax_symbio_meadow)!=0, ]

names_heath <- rownames(tax_symbio_heath.no0)
length(names_heath)
b_heath <- ggplot(tax_symbio_heath.no0, aes(x=log(dna_heath+0.01), y=log(rna_heath+0.01))) +
  geom_point(color = "blue", size = 3) +
  geom_abline(intercept = 0) +
  xlim(0,12) +
  ylim(0,12) +
  theme_bw() +
  geom_label_repel(aes(label = names_heath),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') #+
#ggsave("plotheath.png", width = 5, height = 5)

names_meadow <- rownames(tax_symbio_meadow.no0)
length(names_meadow)  
b_meadow <- ggplot(tax_symbio_meadow.no0, aes(x=log(dna_meadow), y=log(rna_meadow))) +
  geom_point(color = "blue", size = 3) +
  geom_abline(intercept = 0) +
  xlim(0,12) +
  ylim(0,12) +
  theme_bw() +
  geom_label_repel(aes(label = names_meadow),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') #+
#ggsave("plotmeadow.png", width = 5, height = 5)

ggarrange(b_heath, b_meadow, labels = c("heath", "meadow"), common.legend = FALSE)

cor.test(tax_symbio_heath.no0$dna_heath, tax_symbio_heath.no0$rna_heath) #0.851505
cor.test(tax_symbio_meadow.no0$dna_meadow, tax_symbio_meadow.no0$rna_meadow) #0.9393354


#########################################################################################################################
################unassigned in 9 coextracted samples
#########################################################################################################################

table_tax_unassigned <-with(allOTUs_9samples_unassigned, tapply(abundance, list(Order, combined), sum))
table_tax_unassigned
class(table_tax_unassigned)
tax_unassigned <- data.frame(table_tax_unassigned)

tax_unassigned_heath <- tax_unassigned[is.finite(tax_unassigned$dna_heath) & is.finite(tax_unassigned$rna_heath), ]
tax_unassigned_heath <- subset(tax_unassigned_heath, select = c(dna_heath, rna_heath))
tax_unassigned_heath.no0 = tax_unassigned_heath[ rowSums(tax_unassigned_heath)!=0, ]
tax_unassigned_meadow <- tax_unassigned[is.finite(tax_unassigned$dna_meadow) & is.finite(tax_unassigned$rna_meadow), ]
tax_unassigned_meadow <- subset(tax_unassigned_meadow, select = c(dna_meadow, rna_meadow))
tax_unassigned_meadow.no0 = tax_unassigned_meadow[ rowSums(tax_unassigned_meadow)!=0, ]

names_heath <- rownames(tax_unassigned_heath.no0)
length(names_heath)
c_heath <- ggplot(tax_unassigned_heath.no0, aes(x=log(dna_heath), y=log(rna_heath))) +
  geom_point(color = "blue", size = 3) +
  geom_abline(intercept = 0) +
  xlim(0,10) +
  ylim(0,10) +
  theme_bw() +
  geom_label_repel(aes(label = names_heath),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') #+
#ggsave("plotheath.png", width = 5, height = 5)

names_meadow <- rownames(tax_unassigned_meadow.no0)
length(names_meadow)  
c_meadow <- ggplot(tax_unassigned_meadow.no0, aes(x=log(dna_meadow), y=log(rna_meadow))) +
  geom_point(color = "blue", size = 3) +
  geom_abline(intercept = 0) +
  xlim(0,10) +
  ylim(0,10) +
  theme_bw() +
  geom_label_repel(aes(label = names_meadow),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') #+
#ggsave("plotmeadow.png", width = 5, height = 5)

ggarrange(c_heath, c_meadow, labels = c("heath", "meadow"), common.legend = FALSE)

cor.test(tax_unassigned_heath.no0$dna_heath, tax_unassigned_heath.no0$rna_heath) #0.5090624
cor.test(tax_unassigned_meadow.no0$dna_meadow, tax_unassigned_meadow.no0$rna_meadow) #0.5614494

#########################################################################################################################
################single OTU level
#########################################################################################################################

table_tax_singleOTUlevel <-with(allOTUs_9samples, tapply(abundance, list(Base_for_functional_assignment, combined), sum))
table_tax_singleOTUlevel
class(table_tax_singleOTUlevel)
tax_singleOTUlevel <- data.frame(table_tax_singleOTUlevel)
write.table(tax_singleOTUlevel, "tax_singleOTUleve.txt", sep="\t")

tax_singleOTUlevel_heath <- tax_singleOTUlevel[is.finite(tax_singleOTUlevel$dna_heath) & is.finite(tax_singleOTUlevel$rna_heath), ]
tax_singleOTUlevel_heath <- subset(tax_singleOTUlevel_heath, select = c(dna_heath, rna_heath))
tax_singleOTUlevel_heath.no0 = tax_singleOTUlevel_heath[ rowSums(tax_singleOTUlevel_heath)!=0, ]

tax_singleOTUlevel_meadow <- tax_singleOTUlevel[is.finite(tax_singleOTUlevel$dna_meadow) & is.finite(tax_singleOTUlevel$rna_meadow), ]
tax_singleOTUlevel_meadow <- subset(tax_singleOTUlevel_meadow, select = c(dna_meadow, rna_meadow))
tax_singleOTUlevel_meadow.no0 = tax_singleOTUlevel_meadow[ rowSums(tax_singleOTUlevel_meadow)!=0, ]

names_heath <- rownames(tax_singleOTUlevel_heath.no0)
length(names_heath)
  ggplot(tax_singleOTUlevel_heath.no0, aes(x=log(dna_heath+0.01), y=log(rna_heath+0.01))) +
  geom_point(color = "blue", size = 3) +
  geom_abline(intercept = 0) +
  xlim(0,11) +
  ylim(0,11) +
  theme_bw() +
  geom_label_repel(aes(label = names_heath),segment.color = 'grey50') #+
  ggsave("OTUlevel_plotheath.png", width = 20, height = 20)

names_meadow <- rownames(tax_singleOTUlevel_meadow.no0)
length(names_meadow)  
ggplot(tax_singleOTUlevel_meadow.no0, aes(x=log(dna_meadow), y=log(rna_meadow))) +
  geom_point(color = "blue", size = 3) +
  geom_abline(intercept = 0) +
  xlim(0,10) +
  ylim(0,10) +
  theme_bw() +
  geom_label_repel(aes(label = names_meadow), segment.color = 'grey50')+
  ggsave("OTUlevel_plotmeadow.png", width = 20, height = 20)

cor.test(tax_singleOTUlevel_heath.no0$dna_heath, tax_singleOTUlevel_heath.no0$rna_heath) #0.7075866
cor.test(tax_singleOTUlevel_meadow.no0$dna_meadow, tax_singleOTUlevel_meadow.no0$rna_meadow) #0.5614494

######################################################################################
##########################MEAN NO OF READS FROM EACH TROPHIC MODE PER SAMPLE #########
######################################################################################
library(plyr)
#SYMBIOTROPHS:
env$noOfReads_symbio<-rowSums(otus_symbio)
mean_noOfReads <- ddply(env, .(nucleic_acid), summarize, mean=mean(noOfReads_symbio))
#DNA 37981.52
#RNA 33203.32
sdRichness <- ddply(env, .(nucleic_acid), summarize, sd=sd(noOfReads_symbio))
#DNA 5236.142
#RNA 9560.665

env$noOfReads_sapro<-rowSums(otus_sapro)
mean_noOfReads <- ddply(env, .(nucleic_acid), summarize, mean=mean(noOfReads_sapro))
#DNA 1527.522
#RNA 3588.158
sdRichness <- ddply(env, .(nucleic_acid), summarize, sd=sd(noOfReads_sapro))
#DNA 3355.696
#RNA 5334.124

env$noOfReads_patho<-rowSums(otus_patho)
mean_noOfReads <- ddply(env, .(nucleic_acid), summarize, mean=mean(noOfReads_patho))
#DNA  95.30435
#RNA 111.42105
sdRichness <- ddply(env, .(nucleic_acid), summarize, sd=sd(noOfReads_patho))
#DNA 120.7324
#RNA 296.2533

env$noOfReads_unassigned<-rowSums(otus_unass)
mean_noOfReads <- ddply(env, .(nucleic_acid), summarize, mean=mean(noOfReads_unassigned))
#DNA 2883.652
#RNA 5585.105
sdRichness <- ddply(env, .(nucleic_acid), summarize, sd=sd(noOfReads_unassigned))
#DNA 4302.308
#RNA 5297.507

mean(env$noOfReads_patho) #102.5952
sd(env$noOfReads_patho) #215.4502

mean(env$noOfReads_sapro) #2459.714
sd(env$noOfReads_sapro) #4428.478

mean(env$noOfReads_symbio) #35819.95
sd(env$noOfReads_symbio) #7786.864

mean(env$noOfReads_unassigned) #4105.738
sd(env$noOfReads_unassigned)

#richness in 9 samples
library(vegan)
richness_comparison_DNA <- specnumber(t(comparison_DNA))
richness_comparison_RNA <- specnumber(t(comparison_RNA))
cbind(richness_comparison_DNA, richness_comparison_RNA)