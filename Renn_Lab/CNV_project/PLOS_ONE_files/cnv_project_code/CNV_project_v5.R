############################################################################
####             R Code associated with the Manuscript:                 ####
#### Correspondence of aCGH and long-read genome assembly for detection ####
#### of copy number variation: a proof-of-concept with cichlid genomes  ####
####                                                                    ####
####    Gabriel A. Preising, Joshua J. Faber-Hammond, Suzy C.P. Renn    ####
############################################################################

## Last updated 7/15/2021
## Load Packages
library(Hmisc)
library(hexbin)
library(hyperSpec)
library(VennDiagram)
library(ggbiplot)
#install.packages("factoextra")
library(factoextra)
#install.packages("emmeans")
library(emmeans)
library(gridExtra)
library(ggplot2)
#install.packages("compositions")
library(compositions)
#install.packages("RVAideMemoire")
library(RVAideMemoire)
#install.packages("MANOVA.RM")
library(MANOVA.RM)
library(reshape2)
#install.packages("robCompositions")
library(robCompositions)
library(vegan)
#install.packages("scatterplot3d") # Install
library("scatterplot3d") # load
#st.err function for aggregate
f         <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
#install.packages("vioplot")
library(vioplot)
#source("https://bioconductor.org/biocLite.R")
#install.packages("gridExtra")
require(gridExtra)
library (lattice)
library(Hmisc)
library(dplyr)
library(plyr)
library(dplyr)
library(reshape2)
#install.packages("DepthProc")
library(DepthProc)
#install.packages("ddalpha")
library(ddalpha)
library(MASS)
library(eulerr)

## import ndf file, hybe ratio file, and hits per species_seqtech files
ndf <- read.table(file = '120718_Cichlid_EB_CGH_filt.ndf', sep = '\t', header = TRUE)
annot <- read.table(file = 'priority_probe_annotations.txt', sep = '\t', header = TRUE)
hybe <- read.table(file = 'combined_hybe_ratios2.txt', sep = '\t', header = TRUE)
mzi <- read.table(file = 'hits_Mzeb_illumina3.txt', sep = '\t', header = TRUE)
mzp <- read.table(file = 'hits_Mzeb_pacbio3.txt', sep = '\t', header = TRUE)
oni <- read.table(file = 'hits_Oreo_illumina3.txt', sep = '\t', header = TRUE)
onp <- read.table(file = 'hits_Oreo_pacbio3.txt', sep = '\t', header = TRUE)

## merge hits per species_seqtech files
mzeb_hits <- merge(mzi, mzp, all=TRUE, by="PROBE_DESIGN_ID")
onil_hits <- merge(oni, onp, all=TRUE, by="PROBE_DESIGN_ID")
## make final hits per species_seqtech files
ngs_hits <- merge (mzeb_hits, onil_hits, all=TRUE, by="PROBE_DESIGN_ID")

## merge ndf and hybe ratio file
ndf_hybe <- merge(ndf, hybe, all=TRUE, by=c("CHROM", "POSITION"))

## merge previous ndf/hybe ratio file and hit
hybe_v_ngs <- merge (ndf_hybe, ngs_hits, all=TRUE, by="PROBE_DESIGN_ID")

## remove rownames for potential downstream export back to command line
rownames(hybe_v_ngs) <- c()

## exported table as a test
#exptest <- write.table(hybe_v_ngs, file="exptest.tsv", quote=FALSE, sep='\t')

## set NA to 0 in hits per species_seqtech columns
#hybe_v_ngs$HITS_MZEB_ILLUMINA[is.na(hybe_v_ngs$HITS_MZEB_ILLUMINA)] <- 0
#hybe_v_ngs$HITS_MZEB_PACBIO[is.na(hybe_v_ngs$HITS_MZEB_PACBIO)] <- 0
#hybe_v_ngs$HITS_OREO_ILLUMINA[is.na(hybe_v_ngs$HITS_OREO_ILLUMINA)] <- 0
#hybe_v_ngs$HITS_OREO_PACBIO[is.na(hybe_v_ngs$HITS_OREO_PACBIO)] <- 0

## filtering masterfile such that each probe will have at least one hit in each species_technology column
filt_hybe_v_ngs <- hybe_v_ngs[is.na(hybe_v_ngs$HITS_MZEB_ILLUMINA) == FALSE 
                              & is.na(hybe_v_ngs$HITS_MZEB_PACBIO) == FALSE 
                              & is.na(hybe_v_ngs$HITS_OREO_ILLUMINA) == FALSE 
                              & is.na(hybe_v_ngs$HITS_OREO_PACBIO) == FALSE, ]

## create file that contains all non-redundant probe sequences
filt_hybe_v_ngs$dup_check <- duplicated(filt_hybe_v_ngs$PROBE_DESIGN_ID)
filt_hybe_v_ngs_v2 <- filt_hybe_v_ngs[!(filt_hybe_v_ngs$dup_check=="TRUE"),]
## this file contains all probes with more than 1 hit in any assembly
filt_multihit_hybe_v_ngs <- filter(filt_hybe_v_ngs_v2, (HITS_MZEB_PACBIO * HITS_MZEB_ILLUMINA * HITS_OREO_PACBIO * HITS_OREO_ILLUMINA) > 1)
## this file contains all probes with more than 1 hit in either PacBio assembly
filt_multihitPB_hybe_v_ngs <- filter(filt_hybe_v_ngs_v2, (HITS_MZEB_PACBIO * HITS_OREO_PACBIO) > 1)
## this file contains all probes with more than 1 hit either Illumina assembly
filt_multihitIL_hybe_v_ngs <- filter(filt_hybe_v_ngs_v2, (HITS_MZEB_ILLUMINA * HITS_OREO_ILLUMINA) > 1)

## renaming columns
names(hybe_v_ngs)[names(hybe_v_ngs) == "X545383A05_CichlCGH02.WS1_635_qspline"] <- "GCloess_X545383A05_CichlCGH02.WS1_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X545385A05_CichlCGH03.wa05_635_qspline"] <- "GCloess_X545385A05_CichlCGH03.wa05_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X545409A05_CichlCGH05.wa04_635_qspline"] <- "GCloess_X545409A05_CichlCGH05.wa04_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X550971A05_CichlCGH09.wa01_635_qspline"] <- "GCloess_X550971A05_CichlCGH09.wa01_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X545383A05_CichlCGH02.WS1_635_qspline.1"] <- "MA2C_X545383A05_CichlCGH02.WS1_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X545385A05_CichlCGH03.wa05_635_qspline.1"] <- "MA2C_X545385A05_CichlCGH03.wa05_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X545409A05_CichlCGH05.wa04_635_qspline.1"] <- "MA2C_X545409A05_CichlCGH05.wa04_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X550971A05_CichlCGH09.wa01_635_qspline.1"] <- "MA2C_X550971A05_CichlCGH09.wa01_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "chrom"] <- "chrom(dup)"
names(hybe_v_ngs)[names(hybe_v_ngs) == "pos"] <- "pos(dup)"

## Calculate weighted median hybe ratios across samples, using replicate samples as half the weight as singlicate
hybe_v_ngs$GCloess_rep_mean <- apply(hybe_v_ngs[,c("GCloess_X545383A05_CichlCGH02.WS1_635_qspline", "GCloess_X550971A05_CichlCGH09.wa01_635_qspline")],1,mean, na.rm = TRUE)
hybe_v_ngs$MA2C_rep_mean <- apply(hybe_v_ngs[,c("MA2C_X545383A05_CichlCGH02.WS1_635_qspline", "MA2C_X550971A05_CichlCGH09.wa01_635_qspline")],1,mean, na.rm = TRUE)
hybe_v_ngs$GCloess_weighted_median <- apply(hybe_v_ngs[,c("GCloess_rep_mean", "GCloess_X545385A05_CichlCGH03.wa05_635_qspline", "GCloess_X545409A05_CichlCGH05.wa04_635_qspline")],1,median, na.rm = TRUE)
hybe_v_ngs$MA2C_weighted_median <- apply(hybe_v_ngs[,c("MA2C_rep_mean", "MA2C_X545385A05_CichlCGH03.wa05_635_qspline", "MA2C_X545409A05_CichlCGH05.wa04_635_qspline")],1,median, na.rm = TRUE)
### Meta-analysis lead us to use GCloess for primary normalization method in this study ###

## adding columns with Oreo/Mzeb log2 ratios for comparisons between sequencing techs 
## this means log2(HITS_OREO_ILLUMINA/HITS_MZEB_ILLUMINA) and the same for pacbio
hybe_v_ngs$log2_OREOMZEB_ILLUMINA <- log2(hybe_v_ngs$HITS_OREO_ILLUMINA/hybe_v_ngs$HITS_MZEB_ILLUMINA)
hybe_v_ngs$log2_OREOMZEB_PACBIO <- log2(hybe_v_ngs$HITS_OREO_PACBIO/hybe_v_ngs$HITS_MZEB_PACBIO)
## honestly not sure why I added these columns but maybe they'll be useful
hybe_v_ngs$log2_OREO_PACBIOILLUMINA <- log2(hybe_v_ngs$HITS_OREO_PACBIO/hybe_v_ngs$HITS_OREO_ILLUMINA)
hybe_v_ngs$log2_MZEB_PACBIOILLUMINA <- log2(hybe_v_ngs$HITS_MZEB_PACBIO/hybe_v_ngs$HITS_MZEB_ILLUMINA)


## filtering masterfile such that each probe will have at least one hit in each species_technology column
filt_hybe_v_ngs <- hybe_v_ngs[is.na(hybe_v_ngs$HITS_MZEB_ILLUMINA) == FALSE 
                              & is.na(hybe_v_ngs$HITS_MZEB_PACBIO) == FALSE 
                              & is.na(hybe_v_ngs$HITS_OREO_ILLUMINA) == FALSE 
                              & is.na(hybe_v_ngs$HITS_OREO_PACBIO) == FALSE, ]

## Select Probes that have multiple copies in at least one assembly
filt_hybe_v_ngs$dup_check <- duplicated(filt_hybe_v_ngs$PROBE_DESIGN_ID)
filt_hybe_v_ngs_v2 <- filt_hybe_v_ngs[!(filt_hybe_v_ngs$dup_check=="TRUE"),]
filt_multihit_hybe_v_ngs <- filter(filt_hybe_v_ngs_v2, (HITS_MZEB_PACBIO * HITS_MZEB_ILLUMINA * HITS_OREO_PACBIO * HITS_OREO_ILLUMINA) > 1)
## Select Probes that have multiple copies in at least one PacBio assembly
filt_multihitPB_hybe_v_ngs <- filter(filt_hybe_v_ngs_v2, (HITS_MZEB_PACBIO * HITS_OREO_PACBIO) > 1)
## Select Probes that have multiple copies in at least one Illumina assembly
filt_multihitIL_hybe_v_ngs <- filter(filt_hybe_v_ngs_v2, (HITS_MZEB_ILLUMINA * HITS_OREO_ILLUMINA) > 1)

## Select Probes that have multiple copies in all assemblies
#filt_multihit_hybe_v_ngs <- filter(filt_hybe_v_ngs_v2, HITS_MZEB_PACBIO != 1)
#filt_multihit_hybe_v_ngs <- filter(filt_multihit_hybe_v_ngs, HITS_MZEB_ILLUMINA != 1)
#filt_multihit_hybe_v_ngs <- filter(filt_multihit_hybe_v_ngs, HITS_OREO_PACBIO != 1)
#filt_multihit_hybe_v_ngs <- filter(filt_multihit_hybe_v_ngs, HITS_OREO_ILLUMINA != 1)


## Look at correlations of species-bias ratios between technologies
## aCGH vs. PacBio - no cutoff filtering - GCloess
model <- lm(filt_hybe_v_ngs$GCloess_weighted_median~filt_hybe_v_ngs$log2_OREOMZEB_PACBIO)
ggplot(data=filt_hybe_v_ngs, aes(x=GCloess_weighted_median, y=log2_OREOMZEB_PACBIO))+
  geom_point()+
  ggtitle("GCloess_aCGH vs. PacBio - No cutoff filtering")+
  xlab("aCGH log2 ratios")+
  ylab("Pacbio log2 ratios")+
  xlim(-5, 5)+
  ylim(-5, 5)+
  theme_bw()+
  geom_smooth(method=lm)
summary(model)

## repeat with filtered dataset that contains only multi-hit probes
modelB <- lm(filt_multihit_hybe_v_ngs$GCloess_weighted_median~filt_multihit_hybe_v_ngs$log2_OREOMZEB_PACBIO)
ggplot(data=filt_multihit_hybe_v_ngs, aes(x=GCloess_weighted_median, y=log2_OREOMZEB_PACBIO))+
  geom_point()+
  ggtitle("GCloess_aCGH vs. PacBio - No cutoff filtering")+
  xlab("aCGH log2 ratios")+
  ylab("Pacbio log2 ratios")+
  xlim(-5, 5)+
  ylim(-5, 5)+
  theme_bw()+
  geom_smooth(method=lm)
summary(modelB)


## aCGH vs. PacBio: +/- 0.3 cutoff filtering - GCloess
## this removes near-neutral probes
filt_hybe_v_ngs0.3pb_GCloess <- filt_hybe_v_ngs[(filt_hybe_v_ngs$GCloess_weighted_median > 0.3 |
                                         filt_hybe_v_ngs$GCloess_weighted_median < -0.3) &
                                         (filt_hybe_v_ngs$log2_OREOMZEB_PACBIO > 0.3 |
                                         filt_hybe_v_ngs$log2_OREOMZEB_PACBIO < -0.3), ]
model2 <- lm(filt_hybe_v_ngs0.3pb_GCloess$GCloess_weighted_median~filt_hybe_v_ngs0.3pb_GCloess$log2_OREOMZEB_PACBIO)
ggplot(data=filt_hybe_v_ngs0.3pb_GCloess, aes(x=GCloess_weighted_median, y=log2_OREOMZEB_PACBIO))+
  geom_point()+
  ggtitle("GCloess_aCGH vs. PacBio - +/- 0.3")+
  xlab("aCGH log2 ratios")+
  ylab("Pacbio log2 ratios")+
  xlim(-5, 5)+
  ylim(-5, 5)+
  theme_bw()+
  geom_smooth(method=lm)

summary(model2)

## aCGH vs. Illumina - no cutoff filtering - GCloess 
model3 <- lm(filt_hybe_v_ngs$GCloess_weighted_median~filt_hybe_v_ngs$log2_OREOMZEB_ILLUMINA)
ggplot(data=filt_hybe_v_ngs, aes(x=GCloess_weighted_median, y=log2_OREOMZEB_ILLUMINA))+
  geom_point()+
  ggtitle("GCloess_aCGH vs. Illumina - No cutoff filtering")+
  xlab("aCGH log2 ratios")+
  ylab("Illumina log2 ratios")+
  xlim(-5, 5)+
  ylim(-5, 5)+
  theme_bw()+
  geom_smooth(method=lm)
summary(model3)

## repeat with filtered dataset that contains only multi-hit probes
model3B <- lm(filt_multihit_hybe_v_ngs$GCloess_weighted_median~filt_multihit_hybe_v_ngs$log2_OREOMZEB_ILLUMINA)
ggplot(data=filt_multihit_hybe_v_ngs, aes(x=GCloess_weighted_median, y=log2_OREOMZEB_ILLUMINA))+
  geom_point()+
  ggtitle("GCloess_aCGH vs. Illumina - No cutoff filtering")+
  xlab("aCGH log2 ratios")+
  ylab("Illumina log2 ratios")+
  xlim(-5, 5)+
  ylim(-5, 5)+
  theme_bw()+
  geom_smooth(method=lm)
summary(model3B)

## aCGH vs. Illumina: +/- 0.3 cutoff filtering - GCloess
## this removes near-neutral probes
filt_hybe_v_ngs0.3il_GCloess <- filt_hybe_v_ngs[(filt_hybe_v_ngs$GCloess_weighted_median > 0.3 |
                                           filt_hybe_v_ngs$GCloess_weighted_median < -0.3) &
                                          (filt_hybe_v_ngs$log2_OREOMZEB_ILLUMINA > 0.3 |
                                             filt_hybe_v_ngs$log2_OREOMZEB_ILLUMINA < -0.3), ]
model4 <- lm(filt_hybe_v_ngs0.3il_GCloess$GCloess_weighted_median~filt_hybe_v_ngs0.3il_GCloess$log2_OREOMZEB_ILLUMINA)
ggplot(data=filt_hybe_v_ngs0.3il_GCloess, aes(x=GCloess_weighted_median, y=log2_OREOMZEB_ILLUMINA))+
  geom_point()+
  ggtitle("GCloess_aCGH vs. Illumina: +/- 0.3")+
  xlab("aCGH log2 ratios")+
  ylab("Illumina log2 ratios")+
  xlim(-5, 5)+
  ylim(-5, 5)+
  theme_bw()+
  geom_smooth(method=lm)
summary(model4)


## PacBio vs. Illumina - no cutoff filtering 
modelalt <- lm(filt_hybe_v_ngs$log2_OREOMZEB_PACBIO~filt_hybe_v_ngs$log2_OREOMZEB_ILLUMINA)
ggplot(data=filt_hybe_v_ngs, aes(x=log2_OREOMZEB_PACBIO, y=log2_OREOMZEB_ILLUMINA))+
  geom_point()+
  ggtitle("PacBio vs. Illumina - No cutoff filtering")+
  xlab("PacBio log2 ratios")+
  ylab("Illumina log2 ratios")+
  xlim(-5, 5)+
  ylim(-5, 5)+
  theme_bw()+
  geom_smooth(method=lm)
summary(modelalt)

## repeat with filtered dataset that contains only multi-hit probes
modelaltB <- lm(filt_multihit_hybe_v_ngs$log2_OREOMZEB_PACBIO~filt_multihit_hybe_v_ngs$log2_OREOMZEB_ILLUMINA)
ggplot(data=filt_multihit_hybe_v_ngs, aes(x=log2_OREOMZEB_PACBIO, y=log2_OREOMZEB_ILLUMINA))+
  geom_point()+
  ggtitle("PacBio vs. Illumina - No cutoff filtering")+
  xlab("PacBio log2 ratios")+
  ylab("Illumina log2 ratios")+
  xlim(-5, 5)+
  ylim(-5, 5)+
  theme_bw()+
  geom_smooth(method=lm)
summary(modelaltB)


## Onil vs Mzebra - plot between seqtech ratios unfiltered
model5 <- lm(filt_hybe_v_ngs$log2_OREOMZEB_ILLUMINA~filt_hybe_v_ngs$log2_OREOMZEB_PACBIO)
ggplot(data=filt_hybe_v_ngs, aes(x=log2_OREO_PACBIOILLUMINA, y=log2_MZEB_PACBIOILLUMINA))+
  geom_point()+
  ggtitle("ONIL vs. MZEB")+
  xlab("ONIL PacBio vs Illumia")+
  ylab("MZEB PacBio vs Illumia")+
  xlim(-5, 5)+
  ylim(-5, 5)+
  theme_bw()+
  geom_smooth(method=lm)
summary(model5)

## repeat with filtered dataset that contains only multi-hit probes
model5b <- lm(filt_multihit_hybe_v_ngs$log2_OREOMZEB_ILLUMINA~filt_multihit_hybe_v_ngs$log2_OREOMZEB_PACBIO)
ggplot(data=filt_multihit_hybe_v_ngs, aes(x=log2_OREO_PACBIOILLUMINA, y=log2_MZEB_PACBIOILLUMINA))+
  geom_point()+
  ggtitle("ONIL vs. MZEB")+
  xlab("ONIL PacBio vs Illumia")+
  ylab("MZEB PacBio vs Illumia")+
  xlim(-5, 5)+
  ylim(-5, 5)+
  theme_bw()+
  geom_smooth(method=lm)
summary(model5b)


## Apply filters for running platform and species comparisons
## Use +/-0.3 and +/-0.8 filters for hybe ratios
filt_hybe_v_ngs_filt1 <- filt_hybe_v_ngs_v2
filt_hybe_v_ngs_filt2 <- filt_hybe_v_ngs_v2
filt_hybe_v_ngs_filt1[(filt_hybe_v_ngs_filt1 < 0.8) & (filt_hybe_v_ngs_filt1 > -0.8)] <- NA
filt_hybe_v_ngs_filt2[(filt_hybe_v_ngs_filt2 < 0.3) & (filt_hybe_v_ngs_filt2 > -0.3)] <- NA

filt_hybe_v_ngs_ONPB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREO_PACBIOILLUMINA > 0),] #PacBio biased in Onil
filt_hybe_v_ngs_ONIL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREO_PACBIOILLUMINA < 0),] #Illumina biased in Onil
filt_hybe_v_ngs_MZPB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_MZEB_PACBIOILLUMINA > 0),] #PacBio biased in Mzeb
filt_hybe_v_ngs_MZIL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_MZEB_PACBIOILLUMINA < 0),] #Illumina biased in Mzeb

filt_hybe_v_ngs_ON_bias_PB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0),] #Onil biased in PacBio
filt_hybe_v_ngs_ON_bias_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0),] #Onil biased in Illumina
filt_hybe_v_ngs_MZ_bias_PB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < 0),] #Mzeb biased in PacBio
filt_hybe_v_ngs_MZ_bias_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < 0),] #Mzeb biased in Illumina
filt_hybe_v_ngs_ON_bias_aCGH <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median > 0),] #Onil biased in aCGH
filt_hybe_v_ngs_MZ_bias_aCGH <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median < 0),] #Mzeb biased in aCGH

filt_hybe_v_ngs_ON_bias_aCGH_PB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median > 0) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0),] #Onil biased in aCGH and PacBio
filt_hybe_v_ngs_ON_bias_aCGH_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median > 0) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0),] #Onil biased in aCGH and Illumina
filt_hybe_v_ngs_MZ_bias_aCGH_PB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median < 0) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < 0),] #Mzeb biased in aCGH and PacBio
filt_hybe_v_ngs_MZ_bias_aCGH_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median < 0) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < 0),] #Mzeb biased in aCGH and Illumina

filt_hybe_v_ngs_ON_gain_PB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0.8),] #Onil gain in PacBio (past 0.8 threshold)
filt_hybe_v_ngs_ON_gain_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0.8),] #Onil gain in Illumina (past 0.8 threshold)
filt_hybe_v_ngs_MZ_gain_PB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < -0.8),] #Mzeb gain in PacBio (past 0.8 threshold)
filt_hybe_v_ngs_MZ_gain_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < -0.8),] #Mzeb gain in Illumina (past 0.8 threshold)
filt_hybe_v_ngs_ON_gain_aCGH <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median > 0.8),] #Onil gain in aCGH (past 0.8 threshold)
filt_hybe_v_ngs_MZ_gain_aCGH <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median < -0.8),] #Mzeb gain in aCGH (past 0.8 threshold)

filt_hybe_v_ngs_ON_gain_aCGH_PB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0.8),] #Onil gain in aCGH and PacBio (past 0.8 threshold)
filt_hybe_v_ngs_ON_gain_aCGH_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0.8),] #Onil gain in aCGH and Illumina (past 0.8 threshold)
filt_hybe_v_ngs_MZ_gain_aCGH_PB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < -0.8),] #Mzeb gain in aCGH and PacBio (past 0.8 threshold)
filt_hybe_v_ngs_MZ_gain_aCGH_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < -0.8),] #Mzeb gain in aCGH and Illumina (past 0.8 threshold)


## Plot correlations of raw hybes vs PacBio with filters applied
plot(filt_hybe_v_ngs_v2$GCloess_weighted_median, filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO,col = ifelse(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > -0.3 & 0.3 > filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO | filt_hybe_v_ngs_v2$GCloess_weighted_median > -0.3 & 0.3 > filt_hybe_v_ngs_v2$GCloess_weighted_median ,'gray','black'), pch = 16)
par(new=TRUE)
plot(filt_hybe_v_ngs_v2$GCloess_weighted_median, filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO, col = "gray50", pch = 16)
par(new=TRUE)
plot(filt_hybe_v_ngs_v2$GCloess_weighted_median, filt_hybe_v_ngs_filt2$log2_OREOMZEB_PACBIO, col = "gray50", pch = 16)
par(new=TRUE)
plot(filt_hybe_v_ngs_filt2$GCloess_weighted_median, filt_hybe_v_ngs_filt2$log2_OREOMZEB_PACBIO, col = "black", pch = 16)
par(new=TRUE)
plot(filt_hybe_v_ngs_filt1$GCloess_weighted_median, filt_hybe_v_ngs_filt1$log2_OREOMZEB_PACBIO, col = "forestgreen", pch = 16)
abline(h=0, v=0)
abline(lm(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO ~  filt_hybe_v_ngs_v2$GCloess_weighted_median), col="red")
abline(lm(filt_hybe_v_ngs_filt2$log2_OREOMZEB_PACBIO ~  filt_hybe_v_ngs_v2$GCloess_weighted_median), col="blue")
abline(lm(filt_hybe_v_ngs_filt1$log2_OREOMZEB_PACBIO ~  filt_hybe_v_ngs_filt1$GCloess_weighted_median), col="darkgreen")

## Plot correlations of raw hybes vs Illumina with filters applied
plot(filt_hybe_v_ngs_v2$GCloess_weighted_median, filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA,col = ifelse(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > -0.3 & 0.3 > filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA | filt_hybe_v_ngs_v2$GCloess_weighted_median > -0.3 & 0.3 > filt_hybe_v_ngs_v2$GCloess_weighted_median ,'gray','black'), pch = 16)
par(new=TRUE)
plot(filt_hybe_v_ngs_v2$GCloess_weighted_median, filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA, col = "gray50", pch = 16)
par(new=TRUE)
plot(filt_hybe_v_ngs_v2$GCloess_weighted_median, filt_hybe_v_ngs_filt2$log2_OREOMZEB_ILLUMINA, col = "gray50", pch = 16)
par(new=TRUE)
plot(filt_hybe_v_ngs_filt2$GCloess_weighted_median, filt_hybe_v_ngs_filt2$log2_OREOMZEB_ILLUMINA, col = "black", pch = 16)
par(new=TRUE)
plot(filt_hybe_v_ngs_filt1$GCloess_weighted_median, filt_hybe_v_ngs_filt1$log2_OREOMZEB_ILLUMINA, col = "forestgreen", pch = 16)
abline(h=0, v=0)
abline(lm(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA ~  filt_hybe_v_ngs_v2$GCloess_weighted_median), col="red")
abline(lm(filt_hybe_v_ngs_filt2$log2_OREOMZEB_ILLUMINA ~  filt_hybe_v_ngs_v2$GCloess_weighted_median), col="blue")
abline(lm(filt_hybe_v_ngs_filt1$log2_OREOMZEB_ILLUMINA ~  filt_hybe_v_ngs_filt1$GCloess_weighted_median), col="darkgreen")


par(mfrow=c(2,2))
## Plot correlations of raw hybes vs Illumina with PB and Illumina biased (based on ON assembly BLAST hits) probes superimposed
plot(jitter(filt_hybe_v_ngs_v2$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA, factor = 13), col = "gray", pch = 16, ylim=c(-6,6), xlim=c(-6,6))
par(new=TRUE)
plot(jitter(filt_hybe_v_ngs_ONPB$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_ONPB$log2_OREOMZEB_ILLUMINA, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=40, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
par(new=TRUE)
plot(jitter(filt_hybe_v_ngs_ONIL$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_ONIL$log2_OREOMZEB_ILLUMINA, factor = 13), col = rgb(red=100, green=25, blue=25, alpha=35, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
par(new=TRUE)
plot(jitter(filt_hybe_v_ngs_ONPB$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_ONPB$log2_OREOMZEB_ILLUMINA, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=6, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
abline(h=0, v=0)


## Plot correlations of raw hybes vs PacBio with PB and Illumina biased (based on ON assembly BLAST hits) probes superimposed
plot(jitter(filt_hybe_v_ngs_v2$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO, factor = 13), col = "gray", pch = 16, ylim=c(-6,6), xlim=c(-6,6))
par(new=TRUE)
plot(jitter(filt_hybe_v_ngs_ONPB$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_ONPB$log2_OREOMZEB_PACBIO, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=40, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
par(new=TRUE)
plot(jitter(filt_hybe_v_ngs_ONIL$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_ONIL$log2_OREOMZEB_PACBIO, factor = 13), col = rgb(red=100, green=25, blue=25, alpha=35, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
par(new=TRUE)
plot(jitter(filt_hybe_v_ngs_ONPB$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_ONPB$log2_OREOMZEB_PACBIO, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=6, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
abline(h=0, v=0)


## Plot correlations of raw hybes vs Illumina with PB and Illumina biased (based on MZ assembly BLAST hits) probes superimposed
plot(jitter(filt_hybe_v_ngs_v2$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA, factor = 13), col = "gray", pch = 16, ylim=c(-6,6), xlim=c(-6,6))
par(new=TRUE)
plot(jitter(filt_hybe_v_ngs_MZPB$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_MZPB$log2_OREOMZEB_ILLUMINA, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=40, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
par(new=TRUE)
plot(jitter(filt_hybe_v_ngs_MZIL$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_MZIL$log2_OREOMZEB_ILLUMINA, factor = 13), col = rgb(red=100, green=25, blue=25, alpha=35, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
par(new=TRUE)
plot(jitter(filt_hybe_v_ngs_MZPB$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_MZPB$log2_OREOMZEB_ILLUMINA, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=6, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
abline(h=0, v=0)


## Plot correlations of raw hybes vs PacBio with PB and Illumina biased (based on MZ assembly BLAST hits) probes superimposed
plot(jitter(filt_hybe_v_ngs_v2$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO, factor = 13), col = "gray", pch = 16, ylim=c(-6,6), xlim=c(-6,6))
par(new=TRUE)
plot(jitter(filt_hybe_v_ngs_MZPB$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_MZPB$log2_OREOMZEB_PACBIO, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=40, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
par(new=TRUE)
plot(jitter(filt_hybe_v_ngs_MZIL$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_MZIL$log2_OREOMZEB_PACBIO, factor = 13), col = rgb(red=100, green=25, blue=25, alpha=35, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
par(new=TRUE)
plot(jitter(filt_hybe_v_ngs_MZPB$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_MZPB$log2_OREOMZEB_PACBIO, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=6, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
abline(h=0, v=0)
par(mfrow=c(1,1))


## ONPB biased probes in raw hybe vs Illumina comparison (blue on first chart)
dim(filt_hybe_v_ngs_ONPB[(filt_hybe_v_ngs_ONPB$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_ONPB$log2_OREOMZEB_ILLUMINA > 0.8),]) #condordant positive
dim(filt_hybe_v_ngs_ONPB[(filt_hybe_v_ngs_ONPB$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_ONPB$log2_OREOMZEB_ILLUMINA < -0.8),]) #concordant negative
dim(filt_hybe_v_ngs_ONPB[(filt_hybe_v_ngs_ONPB$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_ONPB$log2_OREOMZEB_ILLUMINA < -0.8),]) #discordant aCGH gain / NGS loss
dim(filt_hybe_v_ngs_ONPB[(filt_hybe_v_ngs_ONPB$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_ONPB$log2_OREOMZEB_ILLUMINA > 0.8),]) #discordant aCGH loss / NGS gain

## ONIL biased probes in raw hybe vs Illumina comparison (red on first chart)
dim(filt_hybe_v_ngs_ONIL[(filt_hybe_v_ngs_ONIL$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_ONIL$log2_OREOMZEB_ILLUMINA > 0.8),]) #condordant positive
dim(filt_hybe_v_ngs_ONIL[(filt_hybe_v_ngs_ONIL$GCloess_weighted_median < 0.8) & (filt_hybe_v_ngs_ONIL$log2_OREOMZEB_ILLUMINA < 0.8),]) #concordant negative
dim(filt_hybe_v_ngs_ONIL[(filt_hybe_v_ngs_ONIL$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_ONIL$log2_OREOMZEB_ILLUMINA < 0.8),]) #discordant aCGH gain / NGS loss
dim(filt_hybe_v_ngs_ONIL[(filt_hybe_v_ngs_ONIL$GCloess_weighted_median < 0.8) & (filt_hybe_v_ngs_ONIL$log2_OREOMZEB_ILLUMINA > 0.8),]) #discordant aCGH loss / NGS gain

## ONPB biased probes in raw hybe vs PacBio comparison (blue on second chart)
dim(filt_hybe_v_ngs_ONPB[(filt_hybe_v_ngs_ONPB$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_ONPB$log2_OREOMZEB_PACBIO > 0.8),]) #condordant positive
dim(filt_hybe_v_ngs_ONPB[(filt_hybe_v_ngs_ONPB$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_ONPB$log2_OREOMZEB_PACBIO < -0.8),]) #concordant negative
dim(filt_hybe_v_ngs_ONPB[(filt_hybe_v_ngs_ONPB$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_ONPB$log2_OREOMZEB_PACBIO < -0.8),]) #discordant aCGH gain / NGS loss
dim(filt_hybe_v_ngs_ONPB[(filt_hybe_v_ngs_ONPB$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_ONPB$log2_OREOMZEB_PACBIO > 0.8),]) #discordant aCGH loss / NGS gain

## ONIL biased probes in raw hybe vs Pacbio comparison (red on second chart)
dim(filt_hybe_v_ngs_ONIL[(filt_hybe_v_ngs_ONIL$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_ONIL$log2_OREOMZEB_PACBIO > 0.8),]) #condordant positive
dim(filt_hybe_v_ngs_ONIL[(filt_hybe_v_ngs_ONIL$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_ONIL$log2_OREOMZEB_PACBIO < -0.8),]) #concordant negative
dim(filt_hybe_v_ngs_ONIL[(filt_hybe_v_ngs_ONIL$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_ONIL$log2_OREOMZEB_PACBIO < -0.8),]) #discordant aCGH gain / NGS loss
dim(filt_hybe_v_ngs_ONIL[(filt_hybe_v_ngs_ONIL$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_ONIL$log2_OREOMZEB_PACBIO > 0.8),]) #discordant aCGH loss / NGS gain

## MZPB biased probes in raw hybe vs Illumina comparison (blue on third chart)
dim(filt_hybe_v_ngs_MZPB[(filt_hybe_v_ngs_MZPB$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_MZPB$log2_OREOMZEB_ILLUMINA > 0.8),]) #condordant positive
dim(filt_hybe_v_ngs_MZPB[(filt_hybe_v_ngs_MZPB$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_MZPB$log2_OREOMZEB_ILLUMINA < -0.8),]) #concordant negative
dim(filt_hybe_v_ngs_MZPB[(filt_hybe_v_ngs_MZPB$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_MZPB$log2_OREOMZEB_ILLUMINA < -0.8),]) #discordant aCGH gain / NGS loss
dim(filt_hybe_v_ngs_MZPB[(filt_hybe_v_ngs_MZPB$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_MZPB$log2_OREOMZEB_ILLUMINA > 0.8),]) #discordant aCGH loss / NGS gain

## MZIL biased probes in raw hybe vs Illumina comparison (red on third chart)
dim(filt_hybe_v_ngs_MZIL[(filt_hybe_v_ngs_MZIL$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_MZIL$log2_OREOMZEB_ILLUMINA > 0.8),]) #condordant positive
dim(filt_hybe_v_ngs_MZIL[(filt_hybe_v_ngs_MZIL$GCloess_weighted_median < 0.8) & (filt_hybe_v_ngs_MZIL$log2_OREOMZEB_ILLUMINA < 0.8),]) #concordant negative
dim(filt_hybe_v_ngs_MZIL[(filt_hybe_v_ngs_MZIL$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_MZIL$log2_OREOMZEB_ILLUMINA < 0.8),]) #discordant aCGH gain / NGS loss
dim(filt_hybe_v_ngs_MZIL[(filt_hybe_v_ngs_MZIL$GCloess_weighted_median < 0.8) & (filt_hybe_v_ngs_MZIL$log2_OREOMZEB_ILLUMINA > 0.8),]) #discordant aCGH loss / NGS gain

## MZPB biased probes in raw hybe vs PacBio comparison (blue on fourth chart)
dim(filt_hybe_v_ngs_MZPB[(filt_hybe_v_ngs_MZPB$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_MZPB$log2_OREOMZEB_PACBIO > 0.8),]) #condordant positive
dim(filt_hybe_v_ngs_MZPB[(filt_hybe_v_ngs_MZPB$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_MZPB$log2_OREOMZEB_PACBIO < -0.8),]) #concordant negative
dim(filt_hybe_v_ngs_MZPB[(filt_hybe_v_ngs_MZPB$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_MZPB$log2_OREOMZEB_PACBIO < -0.8),]) #discordant aCGH gain / NGS loss
dim(filt_hybe_v_ngs_MZPB[(filt_hybe_v_ngs_MZPB$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_MZPB$log2_OREOMZEB_PACBIO > 0.8),]) #discordant aCGH loss / NGS gain

## MZIL biased probes in raw hybe vs Pacbio comparison (red on fourth chart)
dim(filt_hybe_v_ngs_MZIL[(filt_hybe_v_ngs_MZIL$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_MZIL$log2_OREOMZEB_PACBIO > 0.8),]) #condordant positive
dim(filt_hybe_v_ngs_MZIL[(filt_hybe_v_ngs_MZIL$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_MZIL$log2_OREOMZEB_PACBIO < -0.8),]) #concordant negative
dim(filt_hybe_v_ngs_MZIL[(filt_hybe_v_ngs_MZIL$GCloess_weighted_median > 0.8) & (filt_hybe_v_ngs_MZIL$log2_OREOMZEB_PACBIO < -0.8),]) #discordant aCGH gain / NGS loss
dim(filt_hybe_v_ngs_MZIL[(filt_hybe_v_ngs_MZIL$GCloess_weighted_median < -0.8) & (filt_hybe_v_ngs_MZIL$log2_OREOMZEB_PACBIO > 0.8),]) #discordant aCGH loss / NGS gain







## Make figure and run Pearson correlations for MS
par(mfrow=c(3,1))
## Plot correlations of aCGH vs Illumina with ON-biased and MZ-biased probes (based on PacBio blast results) superimposed
plot(jitter(filt_multihit_hybe_v_ngs$GCloess_weighted_median, factor = 7), jitter(filt_multihit_hybe_v_ngs$log2_OREOMZEB_ILLUMINA, factor = 13), col = "gray", pch = 16, ylim=c(-6,6), xlim=c(-6,6))
#par(new=TRUE)
#plot(jitter(filt_hybe_v_ngs_ON_bias_PB$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_ON_bias_PB$log2_OREOMZEB_ILLUMINA, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=40, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
#par(new=TRUE)
#plot(jitter(filt_hybe_v_ngs_MZ_bias_PB$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_MZ_bias_PB$log2_OREOMZEB_ILLUMINA, factor = 13), col = rgb(red=100, green=25, blue=25, alpha=35, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
#par(new=TRUE)
#plot(jitter(filt_hybe_v_ngs_ON_bias_PB$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_ON_bias_PB$log2_OREOMZEB_ILLUMINA, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=4, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
abline(h=0, v=0)
abline(lm(filt_multihit_hybe_v_ngs$GCloess_weighted_median ~ filt_multihit_hybe_v_ngs$log2_OREOMZEB_ILLUMINA), col = "red") # plot trendline for all probes with more than one in at least one assembly
#filt_hybe_v_ngs_gains_aCGH <- rbind(filt_hybe_v_ngs_ON_gain_aCGH,filt_hybe_v_ngs_MZ_gain_aCGH)
#filt_hybe_v_ngs_gains_aCGH <- filt_hybe_v_ngs_gains_aCGH[!(filt_hybe_v_ngs_gains_aCGH$log2_OREOMZEB_ILLUMINA==0),]
filt_hybe_v_ngs_filt_gain <- filt_multihit_hybe_v_ngs
filt_hybe_v_ngs_filt_gain <- filt_hybe_v_ngs_filt_gain[(filt_hybe_v_ngs_filt_gain$log2_OREOMZEB_ILLUMINA != 0),]
filt_hybe_v_ngs_filt_gain[(filt_hybe_v_ngs_filt_gain$GCloess_weighted_median < 0.8) & (filt_hybe_v_ngs_filt_gain$GCloess_weighted_median > -0.8),] <- NA
abline(lm(filt_hybe_v_ngs_filt_gain$GCloess_weighted_median ~ filt_hybe_v_ngs_filt_gain$log2_OREOMZEB_ILLUMINA), col = "green4") # plot trendline for all probes identified as gains for two plotted platforms
summary(lm(filt_multihit_hybe_v_ngs$GCloess_weighted_median ~ filt_multihit_hybe_v_ngs$log2_OREOMZEB_ILLUMINA))
summary(lm(filt_hybe_v_ngs_filt_gain$GCloess_weighted_median ~ filt_hybe_v_ngs_filt_gain$log2_OREOMZEB_ILLUMINA))
cor.test( ~ GCloess_weighted_median + log2_OREOMZEB_ILLUMINA,
          data=filt_multihit_hybe_v_ngs,
          method = "pearson")
cor.test( ~ GCloess_weighted_median + log2_OREOMZEB_ILLUMINA,
          data=filt_hybe_v_ngs_filt_gain,
          method = "pearson")

## Plot correlations of aCGH vs PacBio with ON-biased and MZ-biased probes (based on Illumina blast results) superimposed
plot(jitter(filt_multihit_hybe_v_ngs$GCloess_weighted_median, factor = 7), jitter(filt_multihit_hybe_v_ngs$log2_OREOMZEB_PACBIO, factor = 13), col = "gray", pch = 16, ylim=c(-6,6), xlim=c(-6,6))
#par(new=TRUE)
#plot(jitter(filt_hybe_v_ngs_ON_bias_IL$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_ON_bias_IL$log2_OREOMZEB_PACBIO, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=40, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
#par(new=TRUE)
#plot(jitter(filt_hybe_v_ngs_MZ_bias_IL$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_MZ_bias_IL$log2_OREOMZEB_PACBIO, factor = 13), col = rgb(red=100, green=25, blue=25, alpha=35, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
#par(new=TRUE)
#plot(jitter(filt_hybe_v_ngs_ON_bias_IL$GCloess_weighted_median, factor = 7), jitter(filt_hybe_v_ngs_ON_bias_IL$log2_OREOMZEB_PACBIO, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=4, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
abline(h=0, v=0)
abline(lm(filt_multihit_hybe_v_ngs$GCloess_weighted_median ~ filt_multihit_hybe_v_ngs$log2_OREOMZEB_PACBIO), col = "red") # plot trendline for all probes with more than one in at least one assembly
#filt_hybe_v_ngs_gains_aCGH <- rbind(filt_hybe_v_ngs_ON_gain_aCGH,filt_hybe_v_ngs_MZ_gain_aCGH)
#filt_hybe_v_ngs_gains_aCGH <- filt_hybe_v_ngs_gains_aCGH[!(filt_hybe_v_ngs_gains_aCGH$log2_OREOMZEB_PACBIO==0),]
filt_hybe_v_ngs_filt_gain2 <- filt_multihit_hybe_v_ngs
filt_hybe_v_ngs_filt_gain2 <- filt_hybe_v_ngs_filt_gain2[(filt_hybe_v_ngs_filt_gain2$log2_OREOMZEB_PACBIO != 0),]
filt_hybe_v_ngs_filt_gain2[(filt_hybe_v_ngs_filt_gain2$GCloess_weighted_median < 0.8) & (filt_hybe_v_ngs_filt_gain2$GCloess_weighted_median > -0.8),] <- NA
abline(lm(filt_hybe_v_ngs_filt_gain2$GCloess_weighted_median ~ filt_hybe_v_ngs_filt_gain2$log2_OREOMZEB_PACBIO), col = "green4") # plot trendline for all probes identified as gains for two plotted platforms
summary(lm(filt_multihit_hybe_v_ngs$GCloess_weighted_median ~ filt_multihit_hybe_v_ngs$log2_OREOMZEB_PACBIO))
summary(lm(filt_hybe_v_ngs_filt_gain2$GCloess_weighted_median ~ filt_hybe_v_ngs_filt_gain2$log2_OREOMZEB_PACBIO))
cor.test( ~ GCloess_weighted_median + log2_OREOMZEB_PACBIO,
          data=filt_multihit_hybe_v_ngs,
          method = "pearson")
cor.test( ~ GCloess_weighted_median + log2_OREOMZEB_PACBIO,
          data=filt_hybe_v_ngs_filt_gain2,
          method = "pearson")

## Plot correlations of PacBio vs Illumina with ON-biased and MZ-biased probes (based on aCGH results) superimposed
plot(jitter(filt_multihit_hybe_v_ngs$log2_OREOMZEB_PACBIO, factor = 13), jitter(filt_multihit_hybe_v_ngs$log2_OREOMZEB_ILLUMINA, factor = 13), col = "gray", pch = 16, ylim=c(-6,6), xlim=c(-6,6))
#par(new=TRUE)
#plot(jitter(filt_hybe_v_ngs_ON_bias_aCGH$log2_OREOMZEB_PACBIO, factor = 13), jitter(filt_hybe_v_ngs_ON_bias_aCGH$log2_OREOMZEB_ILLUMINA, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=40, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
#par(new=TRUE)
#plot(jitter(filt_hybe_v_ngs_MZ_bias_aCGH$log2_OREOMZEB_PACBIO, factor = 13), jitter(filt_hybe_v_ngs_MZ_bias_aCGH$log2_OREOMZEB_ILLUMINA, factor = 13), col = rgb(red=100, green=25, blue=25, alpha=35, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
#par(new=TRUE)
#plot(jitter(filt_hybe_v_ngs_ON_bias_aCGH$log2_OREOMZEB_PACBIO, factor = 13), jitter(filt_hybe_v_ngs_ON_bias_aCGH$log2_OREOMZEB_ILLUMINA, factor = 13), col = rgb(red=35, green=35, blue=100, alpha=4, maxColorValue=100), pch = 16, ylim=c(-6,6), xlim=c(-6,6))
abline(h=0, v=0)
abline(lm(filt_multihit_hybe_v_ngs$log2_OREOMZEB_PACBIO ~ filt_multihit_hybe_v_ngs$log2_OREOMZEB_ILLUMINA), col = "red") # plot trendline for all probes with more than one in at least one assembly
#filt_hybe_v_ngs_biases <- filt_multihit_hybe_v_ngs[!(filt_multihit_hybe_v_ngs$log2_OREOMZEB_ILLUMINA==0 | filt_multihit_hybe_v_ngs$log2_OREOMZEB_PACBIO==0),]
filt_hybe_v_ngs_filt_gain3 <- filt_multihit_hybe_v_ngs
filt_hybe_v_ngs_filt_gain3 <- filt_hybe_v_ngs_filt_gain3[(filt_hybe_v_ngs_filt_gain3$log2_OREOMZEB_ILLUMINA != 0),]
filt_hybe_v_ngs_filt_gain3 <- filt_hybe_v_ngs_filt_gain3[(filt_hybe_v_ngs_filt_gain3$log2_OREOMZEB_PACBIO != 0),]
abline(lm(filt_hybe_v_ngs_filt_gain3$log2_OREOMZEB_PACBIO ~ filt_hybe_v_ngs_filt_gain3$log2_OREOMZEB_ILLUMINA), col = "green4") # plot trendline for all probes identified as gains for two plotted platforms
summary(lm(filt_multihit_hybe_v_ngs$log2_OREOMZEB_PACBIO ~ filt_multihit_hybe_v_ngs$log2_OREOMZEB_ILLUMINA))
summary(lm(filt_hybe_v_ngs_filt_gain3$log2_OREOMZEB_PACBIO ~ filt_hybe_v_ngs_filt_gain3$log2_OREOMZEB_ILLUMINA))
cor.test( ~ log2_OREOMZEB_PACBIO + log2_OREOMZEB_ILLUMINA,
          data=filt_multihit_hybe_v_ngs,
          method = "pearson")
cor.test( ~ log2_OREOMZEB_PACBIO + log2_OREOMZEB_ILLUMINA,
          data=filt_hybe_v_ngs_filt_gain3,
          method = "pearson")
par(mfrow=c(1,1))


## Check numbers of ON gains identified by all 3 technologies and their overlap
dim(filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median > 0.8),]) # Number of probes as ON gains in aCGH only
dim(filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0.8),]) # Number of probes as ON gains PacBio only
dim(filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0.8),]) # Number of probes as ON gains Illumina only
dim(filt_hybe_v_ngs_ON_gain_aCGH_PB) # Number of probes as ON gains in aCGH and PacBio
dim(filt_hybe_v_ngs_ON_gain_aCGH_IL) # Number of probes as ON gains in aCGH and Illumina
dim(filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0.8) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0.8),]) # Number of probes as ON gains in PacBio and Illumina
dim(filt_hybe_v_ngs_ON_gain_aCGH_IL[(filt_hybe_v_ngs_ON_gain_aCGH_IL$log2_OREOMZEB_PACBIO > 0.8),]) # Number of probes as ON gains in all 3 technologies
dim(filt_hybe_v_ngs_ON_gain_aCGH_PB[(filt_hybe_v_ngs_ON_gain_aCGH_PB$log2_OREOMZEB_ILLUMINA > 0.8),]) # This should match the number above, just a double check


## Check numbers of MZ gains identified by all 3 technologies and their overlap
dim(filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median < -0.8),]) # Number of probes as MZ gains in aCGH only
dim(filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < -0.8),])  # Number of probes as MZ gains in PacBio only
dim(filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < -0.8),]) # Number of probes as MZ gains in Illumina only
dim(filt_hybe_v_ngs_MZ_gain_aCGH_PB) # Number of probes as MZ gains in aCGH and PacBio
dim(filt_hybe_v_ngs_MZ_gain_aCGH_IL) # Number of probes as MZ gains in aCGH and Illumina
dim(filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < -0.8) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < -0.8),]) # Number of probes as MZ gains in PacBio and Illumina
dim(filt_hybe_v_ngs_MZ_gain_aCGH_IL[(filt_hybe_v_ngs_MZ_gain_aCGH_IL$log2_OREOMZEB_PACBIO < -0.8),]) # Number of probes as MZ gains in all 3 technologies
dim(filt_hybe_v_ngs_MZ_gain_aCGH_PB[(filt_hybe_v_ngs_MZ_gain_aCGH_PB$log2_OREOMZEB_ILLUMINA < -0.8),]) # This should match the number above, just a double check

## export probes from ON gains identified by all 3 technologies and their overlap
ON_gains_aCGH <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median > 0.8),] # Number of probes as ON gains in aCGH only
ON_gains_PB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0.8),] # Number of probes as ON gains PacBio only
ON_gains_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0.8),] # Number of probes as ON gains Illumina only
ON_gains_aCGH_PB <- filt_hybe_v_ngs_ON_gain_aCGH_PB # Number of probes as ON gains in aCGH and PacBio
ON_gains_aCGH_IL <- filt_hybe_v_ngs_ON_gain_aCGH_IL # Number of probes as ON gains in aCGH and Illumina
ON_gains_PB_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0.8) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0.8),] # Number of probes as ON gains in PacBio and Illumina
ON_gains_aCGH_PB_IL <- filt_hybe_v_ngs_ON_gain_aCGH_IL[(filt_hybe_v_ngs_ON_gain_aCGH_IL$log2_OREOMZEB_PACBIO > 0.8),] # Number of probes as ON gains in all 3 technologies
#write.table(ON_gains_aCGH, file = "./probe_lists/ON_gains_aCGH.txt", sep = "\t", row.names = FALSE)
#write.table(ON_gains_PB, file = "./probe_lists/ON_gains_PB.txt", sep = "\t", row.names = FALSE)
#write.table(ON_gains_IL, file = "./probe_lists/ON_gains_IL.txt", sep = "\t", row.names = FALSE)
#write.table(ON_gains_aCGH_PB, file = "./probe_lists/ON_gains_aCGH_PB.txt", sep = "\t", row.names = FALSE)
#write.table(ON_gains_aCGH_IL, file = "./probe_lists/ON_gains_aCGH_IL.txt", sep = "\t", row.names = FALSE)
#write.table(ON_gains_PB_IL, file = "./probe_lists/ON_gains_PB_IL.txt", sep = "\t", row.names = FALSE)
#write.table(ON_gains_aCGH_PB_IL, file = "./probe_lists/ON_gains_aCGH_PB_IL.txt", sep = "\t", row.names = FALSE)

## export probes from ON bias identified by all 3 technologies and their overlap
ON_bias_aCGH <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median > 0),] # Number of probes as ON bias in aCGH only
ON_bias_PB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0),] # Number of probes as ON bias PacBio only
ON_bias_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0),] # Number of probes as ON bias Illumina only
ON_bias_aCGH_PB <- filt_hybe_v_ngs_ON_bias_aCGH_PB # Number of probes as ON bias in aCGH and PacBio
ON_bias_aCGH_IL <- filt_hybe_v_ngs_ON_bias_aCGH_IL # Number of probes as ON bias in aCGH and Illumina
ON_bias_PB_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0),] # Number of probes as ON bias in PacBio and Illumina
ON_bias_aCGH_PB_IL <- filt_hybe_v_ngs_ON_gain_aCGH_IL[(filt_hybe_v_ngs_ON_gain_aCGH_IL$log2_OREOMZEB_PACBIO > 0),] # Number of probes as ON bias in all 3 technologies
#write.table(ON_bias_aCGH, file = "./probe_lists/ON_bias_aCGH.txt", sep = "\t", row.names = FALSE)
#write.table(ON_bias_PB, file = "./probe_lists/ON_bias_PB.txt", sep = "\t", row.names = FALSE)
#write.table(ON_bias_IL, file = "./probe_lists/ON_bias_IL.txt", sep = "\t", row.names = FALSE)
#write.table(ON_bias_aCGH_PB, file = "./probe_lists/ON_bias_aCGH_PB.txt", sep = "\t", row.names = FALSE)
#write.table(ON_bias_aCGH_IL, file = "./probe_lists/ON_bias_aCGH_IL.txt", sep = "\t", row.names = FALSE)
#write.table(ON_bias_PB_IL, file = "./probe_lists/ON_bias_PB_IL.txt", sep = "\t", row.names = FALSE)
#write.table(ON_bias_aCGH_PB_IL, file = "./probe_lists/ON_bias_aCGH_PB_IL.txt", sep = "\t", row.names = FALSE)

#export probes from MZ gains identified by all 3 technologies and their overlap
MZ_gains_aCGH <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median < -0.8),] # Number of probes as MZ gains in aCGH only
MZ_gains_PB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < -0.8),]  # Number of probes as MZ gains in PacBio only
MZ_gains_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < -0.8),] # Number of probes as MZ gains in Illumina only
MZ_gains_aCGH_PB <- filt_hybe_v_ngs_MZ_gain_aCGH_PB # Number of probes as MZ gains in aCGH and PacBio
MZ_gains_aCGH_IL <- filt_hybe_v_ngs_MZ_gain_aCGH_IL # Number of probes as MZ gains in aCGH and Illumina
MZ_gains_PB_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < -0.8) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < -0.8),] # Number of probes as MZ gains in PacBio and Illumina
MZ_gains_aCGH_PB_IL <- filt_hybe_v_ngs_MZ_gain_aCGH_IL[(filt_hybe_v_ngs_MZ_gain_aCGH_IL$log2_OREOMZEB_PACBIO < -0.8),] # Number of probes as MZ gains in all 3 technologies
#write.table(MZ_gains_aCGH, file = "./probe_lists/MZ_gains_aCGH.txt", sep = "\t", row.names = FALSE)
#write.table(MZ_gains_PB, file = "./probe_lists/MZ_gains_PB.txt", sep = "\t", row.names = FALSE)
#write.table(MZ_gains_IL, file = "./probe_lists/MZ_gains_IL.txt", sep = "\t", row.names = FALSE)
#write.table(MZ_gains_aCGH_PB, file = "./probe_lists/MZ_gains_aCGH_PB.txt", sep = "\t", row.names = FALSE)
#write.table(MZ_gains_aCGH_IL, file = "./probe_lists/MZ_gains_aCGH_IL.txt", sep = "\t", row.names = FALSE)
#write.table(MZ_gains_PB_IL, file = "./probe_lists/MZ_gains_PB_IL.txt", sep = "\t", row.names = FALSE)
#write.table(MZ_gains_aCGH_PB_IL, file = "./probe_lists/MZ_gains_aCGH_PB_IL.txt", sep = "\t", row.names = FALSE)

#export probes from MZ bias identified by all 3 technologies and their overlap
MZ_bias_aCGH <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$GCloess_weighted_median < 0),] # Number of probes as MZ bias in aCGH only
MZ_bias_PB <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < 0),]  # Number of probes as MZ bias in PacBio only
MZ_bias_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < 0),] # Number of probes as MZ bias in Illumina only
MZ_bias_aCGH_PB <- filt_hybe_v_ngs_MZ_bias_aCGH_PB # Number of probes as MZ bias in aCGH and PacBio
MZ_bias_aCGH_IL <- filt_hybe_v_ngs_MZ_bias_aCGH_IL # Number of probes as MZ bias in aCGH and Illumina
MZ_bias_PB_IL <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < 0) & (filt_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < 0),] # Number of probes as MZ bias in PacBio and Illumina
MZ_bias_aCGH_PB_IL <- filt_hybe_v_ngs_MZ_gain_aCGH_IL[(filt_hybe_v_ngs_MZ_gain_aCGH_IL$log2_OREOMZEB_PACBIO < 0),] # Number of probes as MZ bias in all 3 technologies
#write.table(MZ_bias_aCGH, file = "./probe_lists/MZ_bias_aCGH.txt", sep = "\t", row.names = FALSE)
#write.table(MZ_bias_PB, file = "./probe_lists/MZ_bias_PB.txt", sep = "\t", row.names = FALSE)
#write.table(MZ_bias_IL, file = "./probe_lists/MZ_bias_IL.txt", sep = "\t", row.names = FALSE)
#write.table(MZ_bias_aCGH_PB, file = "./probe_lists/MZ_bias_aCGH_PB.txt", sep = "\t", row.names = FALSE)
#write.table(MZ_bias_aCGH_IL, file = "./probe_lists/MZ_bias_aCGH_IL.txt", sep = "\t", row.names = FALSE)
#write.table(MZ_bias_PB_IL, file = "./probe_lists/MZ_bias_PB_IL.txt", sep = "\t", row.names = FALSE)
#write.table(MZ_bias_aCGH_PB_IL, file = "./probe_lists/MZ_bias_aCGH_PB_IL.txt", sep = "\t", row.names = FALSE)

#### Euler diagram time! ####
## First, create euler plots for probes
MZ_bias_PB_list <- MZ_bias_PB$PROBE_DESIGN_ID
MZ_bias_IL_list <- MZ_bias_IL$PROBE_DESIGN_ID
MZ_gains_aCGH_list <- MZ_gains_aCGH$PROBE_DESIGN_ID
MZ_gains_aCGH_list <- MZ_gains_aCGH_list[!is.na(MZ_gains_aCGH_list)]

ON_bias_PB_list <- ON_bias_PB$PROBE_DESIGN_ID
ON_bias_IL_list <- ON_bias_IL$PROBE_DESIGN_ID
ON_gains_aCGH_list <- ON_gains_aCGH$PROBE_DESIGN_ID
ON_gains_aCGH_list <- ON_gains_aCGH_list[!is.na(ON_gains_aCGH_list)]

##for CD-Hit clustering, export probe lists
#write.table(MZ_bias_PB_list, './probe_lists/MZ_bias_PB_list.txt', sep='\t', row.names = FALSE )
#write.table(MZ_bias_IL_list, './probe_lists/MZ_bias_IL_list.txt', sep='\t', row.names = FALSE )
#write.table(MZ_gains_aCGH_list, './probe_lists/MZ_gains_aCGH_list.txt', sep='\t', row.names = FALSE )
#write.table(ON_bias_PB_list, './probe_lists/ON_bias_PB_list.txt', sep='\t', row.names = FALSE )
#write.table(ON_bias_IL_list, './probe_lists/ON_bias_IL_list.txt', sep='\t', row.names = FALSE )
#write.table(ON_gains_aCGH_list, './probe_lists/ON_gains_aCGH_list.txt', sep='\t', row.names = FALSE )

MZ_venn <- list("MZ PacBio probes" = MZ_bias_PB_list,
                "MZ Illumina probes" = MZ_bias_IL_list,
                "MZ aCGH probes" = MZ_gains_aCGH_list)
MZ_euler_probe <- plot(euler(MZ_venn), quantities = TRUE)


ON_venn <- list("ON PacBio probes" = ON_bias_PB_list,
                "ON Illumina probes" = ON_bias_IL_list,
                "ON aCGH probes" = ON_gains_aCGH_list)
ON_euler_probe <- plot(euler(ON_venn), quantities = TRUE)


## Find how many genic probes in venn
MZ_bias_PB2 <- merge(MZ_bias_PB, annot, all.x=TRUE, by="PROBE_DESIGN_ID")
MZ_bias_PB_genic <- filter(MZ_bias_PB2, Priority_probe_type != "BACKBONE_PROBE")
MZ_bias_IL2 <- merge(MZ_bias_IL, annot, all.x=TRUE, by="PROBE_DESIGN_ID")
MZ_bias_IL_genic <- filter(MZ_bias_IL2, Priority_probe_type != "BACKBONE_PROBE")
MZ_gains_aCGH2 <- merge(MZ_gains_aCGH, annot, all.x=TRUE, by="PROBE_DESIGN_ID")
MZ_gains_aCGH_genic <- filter(MZ_gains_aCGH2, Priority_probe_type != "BACKBONE_PROBE")

MZ_bias_PB_list_genic <- unique(MZ_bias_PB_genic$PROBE_DESIGN_ID)
MZ_bias_IL_list_genic <- unique(MZ_bias_IL_genic$PROBE_DESIGN_ID)
MZ_gains_aCGH_list_genic <- MZ_gains_aCGH_genic$PROBE_DESIGN_ID
MZ_gains_aCGH_list_genic <- unique(MZ_gains_aCGH_list_genic[!is.na(MZ_gains_aCGH_list_genic)])
 
ON_bias_PB2 <- merge(ON_bias_PB, annot, all.x=TRUE, by="PROBE_DESIGN_ID")
ON_bias_PB_genic <- filter(ON_bias_PB2, Priority_probe_type != "BACKBONE_PROBE")
ON_bias_IL2 <- merge(ON_bias_IL, annot, all.x=TRUE, by="PROBE_DESIGN_ID")
ON_bias_IL_genic <- filter(ON_bias_IL2, Priority_probe_type != "BACKBONE_PROBE")
ON_gains_aCGH2 <- merge(ON_gains_aCGH, annot, all.x=TRUE, by="PROBE_DESIGN_ID")
ON_gains_aCGH_genic <- filter(ON_gains_aCGH2, Priority_probe_type != "BACKBONE_PROBE")

ON_bias_PB_list_genic <- unique(ON_bias_PB_genic$PROBE_DESIGN_ID)
ON_bias_IL_list_genic <- unique(ON_bias_IL_genic$PROBE_DESIGN_ID)
ON_gains_aCGH_list_genic <- ON_gains_aCGH_genic$PROBE_DESIGN_ID
ON_gains_aCGH_list_genic <- unique(ON_gains_aCGH_list_genic[!is.na(ON_gains_aCGH_list_genic)])

MZ_genic_venn <- list("MZ PacBio genic probes" = MZ_bias_PB_list_genic,
                "MZ Illumina genic probes" = MZ_bias_IL_list_genic,
                "MZ aCGH genic probes" = MZ_gains_aCGH_list_genic)
MZ_genic_euler_probe <- plot(euler(MZ_genic_venn), quantities = TRUE)


ON_genic_venn <- list("ON PacBio genic probes" = ON_bias_PB_list_genic,
                "ON Illumina genic probes" = ON_bias_IL_list_genic,
                "ON aCGH genic probes" = ON_gains_aCGH_list_genic)
ON_genic_euler_probe <- plot(euler(ON_genic_venn), quantities = TRUE)


### Now look at lists of unique annotations per set

MZ_gene_venn <- MZ_gene_venn <- list("MZ PacBio genes" = unique(MZ_bias_PB_genic$Priority_accession),
                                     "MZ Illumina genes" = unique(MZ_bias_IL_genic$Priority_accession),
                                     "MZ aCGH genes" = unique(MZ_gains_aCGH_genic$Priority_accession))
MZ_euler_gene <- plot(euler(MZ_gene_venn), quantities = TRUE)

ON_gene_venn <- list("ON PacBio genes" = unique(ON_bias_PB_genic$Priority_accession),
                     "ON Illumina genes" = unique(ON_bias_IL_genic$Priority_accession),
                     "ON aCGH genes" = unique(ON_gains_aCGH_genic$Priority_accession))
ON_euler_gene <- plot(euler(ON_gene_venn), quantities = TRUE)

## Plot all euler plots together
gridExtra::grid.arrange(MZ_euler_probe, MZ_euler_gene, ON_euler_probe, ON_euler_gene)

## Plot all euler plots together
gridExtra::grid.arrange(MZ_euler_probe, ON_euler_probe, MZ_genic_euler_probe, ON_genic_euler_probe, MZ_euler_gene, ON_euler_gene)

#write.table(MZ_bias_PB_list_genic, './probe_lists/MZ_genic_bias_PB_list.txt', sep='\t', row.names = FALSE )
#write.table(MZ_bias_IL_list_genic, './probe_lists/MZ_genic_bias_IL_list.txt', sep='\t', row.names = FALSE )
#write.table(MZ_gains_aCGH_list_genic, './probe_lists/MZ_genic_gains_aCGH_list.txt', sep='\t', row.names = FALSE )
#write.table(ON_bias_PB_list_genic, './probe_lists/ON_genic_bias_PB_list.txt', sep='\t', row.names = FALSE )
#write.table(ON_bias_IL_list_genic, './probe_lists/ON_genic_bias_IL_list.txt', sep='\t', row.names = FALSE )
#write.table(ON_gains_aCGH_list_genic, './probe_lists/ON_genic_gains_aCGH_list.txt', sep='\t', row.names = FALSE )

## PCA
acgh_PB.pca <- prcomp(na.omit(filt_hybe_v_ngs_v2[,c(23,25,26)]), scale. = TRUE)
acgh_PB.pca <- princomp(na.omit(filt_hybe_v_ngs_v2[,c(23,25,26)]), cor = TRUE)
summary(acgh_PB.pca)
ggbiplot(acgh_PB.pca)
acgh_PB.pca$loadings

## Import probe nucleotide compostion stats generated in Mesquite
## These stats are for probes identified as gains in either species through aCGH (+/-0.8 threshold) or biases in NGS assemblies (<>0) 
comp_freq <- read.table(file = 'DNA_comp_frequency_stats.txt', sep = '\t', header = TRUE)
## Now load lists of probes ID's as CNVs by 1, 2, or 3 platforms
gain_cat <- read.table(file = './probe_lists/gain_categories2.txt', sep = '\t', header = FALSE)
gain_cat <- cbind.data.frame(probeType = gain_cat$V1, probeID = gain_cat$V2)
gain_comp_freq <- join(gain_cat,comp_freq, by = "probeID")

gain_comp_freq$dup_check <- duplicated(gain_comp_freq$probeID) #ID second entry of duplicate probes
gain_comp_freq$dup_check2 <- duplicated(gain_comp_freq$probeID, fromLast = TRUE) #ID first entry of duplicate probes
gain_comp_freq_nr <- gain_comp_freq[!(gain_comp_freq$dup_check=="TRUE"),] #Remove duplicates
gain_comp_freq_nr <- gain_comp_freq_nr[!(gain_comp_freq_nr$dup_check2=="TRUE"),] #Remove duplicates

## format for PCA and MANOVA
gain_comp_freq2 <- gain_comp_freq_nr[,-2]
rownames(gain_comp_freq2) <- gain_comp_freq_nr[,2]
#gain_comp_freq2 <- gain_comp_freq2[,-1]
gain_comp_freq2 <- gain_comp_freq2[,-25]
gain_comp_freq2 <- gain_comp_freq2[,-24]
dim(gain_comp_freq2)

#comp_freq.pca <- prcomp(gain_comp_freq2[,c(2:23)], scale. = TRUE)
comp_freq.pca <- princomp(gain_comp_freq2[,c(2:23)], cor = TRUE)
summary(comp_freq.pca)
#PC1 and PC2
ggbiplot(comp_freq.pca, groups = gain_comp_freq2$probeType, obs.scale = 1, var.scale = 1, ellipse = TRUE, ellipse.prob = .75, alpha = 0, var.axes = TRUE)
#PC2 and PC3
ggbiplot(comp_freq.pca, choices = c(2,3), groups = gain_comp_freq2$probeType, obs.scale = 1, var.scale = 1, ellipse = TRUE, ellipse.prob = .75, alpha = 0, var.axes = TRUE)
comp_freq.pca$loadings
summary(comp_freq.pca)
#write.table(comp_freq.pca$loadings, file="probe_sets_PCA_loadings.txt", quote=FALSE, sep='\t')

## Filter out all probes found by multiple technologies
gain_comp_freq3 <- gain_comp_freq2[!grepl("PB_aCGH", gain_comp_freq2$probeType),]
gain_comp_freq3 <- gain_comp_freq3[!grepl("IL_aCGH", gain_comp_freq3$probeType),]
gain_comp_freq3 <- gain_comp_freq3[!grepl("IL_PB", gain_comp_freq3$probeType),]

## PCA for filtered probes
#comp_freq3.pca <- prcomp(gain_comp_freq3[,c(2:23)], scale. = TRUE)
comp_freq3.pca <- princomp(gain_comp_freq3[,c(2:23)], cor = TRUE)
summary(comp_freq3.pca)
#PC1 and PC2
ggbiplot(comp_freq3.pca, groups = gain_comp_freq3$probeType, obs.scale = 1, var.scale = 1, ellipse = TRUE, ellipse.prob = .95, alpha = 0, var.axes = TRUE, repel = TRUE)
#PC2 and PC3
ggbiplot(comp_freq3.pca, choices = c(2,3), groups = gain_comp_freq3$probeType, obs.scale = 1, var.scale = 1, ellipse = TRUE, ellipse.prob = .95, alpha = 0, var.axes = TRUE, repel = TRUE)
comp_freq3.pca$loadings
summary(comp_freq3.pca)
#write.table(comp_freq3.pca$loadings, file="platform_excl_PCA_loadings.txt", quote=FALSE, sep='\t')


## Check for pairwise differences in nucleotide composition between all probe categories
res.perm.freq <- pairwise.perm.manova(dist(gain_comp_freq[, 3:24],"euclidean"),gain_comp_freq$probeType, nperm = 100,
                                      progress = TRUE, p.method = "fdr")
res.perm.freq
summary(res.perm.freq)

## Check for pairwise differences in nucleotide composition between groups of probes found by only one technology
res.perm.freq_pltfm_excl <- pairwise.perm.manova(dist(gain_comp_freq3[, 2:23],"euclidean"),gain_comp_freq3$probeType, nperm = 100,
                                      progress = TRUE, p.method = "fdr")
res.perm.freq_pltfm_excl
summary(res.perm.freq_pltfm_excl)


## Post-hoc ANOVAs to further investigate specific nucleotide/dinucleotide frequency differences between platforms
gain_comp_freq_aCGH_IL <- gain_comp_freq3[gain_comp_freq3$probeType != 'PB_list',]
gain_comp_freq_aCGH_PB <- gain_comp_freq3[gain_comp_freq3$probeType != 'IL_list',]
gain_comp_freq_PB_IL <- gain_comp_freq3[gain_comp_freq3$probeType != 'aCGH_list',]

for(i in 2:ncol(gain_comp_freq_aCGH_IL))
{
  
  column <- names(gain_comp_freq_aCGH_IL[i])
  #tidy will summarise and return neat format
  aov_list <- broom::tidy(aov(gain_comp_freq_aCGH_IL[,i] ~ probeType, data = gain_comp_freq_aCGH_IL))
  
  # Add this condition if you only want aov with P < 0.05 printed
  if(aov_list$p.value[1] < 0.05) {
    
    print(column)
    print(aov_list)
  }
}

for(i in 2:ncol(gain_comp_freq_aCGH_PB))
{
  
  column <- names(gain_comp_freq_aCGH_PB[i])
  #tidy will summarise and return neat format
  aov_list <- broom::tidy(aov(gain_comp_freq_aCGH_PB[,i] ~ probeType, data = gain_comp_freq_aCGH_PB))
  
  # Add this condition if you only want aov with P < 0.05 printed
  if(aov_list$p.value[1] < 0.05) {
    
    print(column)
    print(aov_list)
  }
}

for(i in 2:ncol(gain_comp_freq_PB_IL))
{
  
  column <- names(gain_comp_freq_PB_IL[i])
  #tidy will summarise and return neat format
  aov_list <- broom::tidy(aov(gain_comp_freq_PB_IL[,i] ~ probeType, data = gain_comp_freq_PB_IL))
  
  # Add this condition if you only want aov with P < 0.05 printed
  if(aov_list$p.value[1] < 0.05) {
    
    print(column)
    print(aov_list)
  }
}


##########################################
## Boxplot to look at aCGH ratios of NGS-biased probes
par(mfrow=c(1,1))
full_aCGH <- as.data.frame(cbind(as.numeric(filt_hybe_v_ngs_v2$GCloess_weighted_median),"fulldata"))

ONPBbias_aCGH <- as.data.frame(cbind(as.numeric(filt_hybe_v_ngs_ONPB$GCloess_weighted_median),"ONPBbias"))
ONILbias_aCGH <- as.data.frame(cbind(as.numeric(filt_hybe_v_ngs_ONIL$GCloess_weighted_median),"ONILbias"))
ON_TechBias_aCGH <- rbind(ONPBbias_aCGH,ONILbias_aCGH)
ON_TechBias_aCGH <- cbind.data.frame(log_ratio=as.numeric(as.character(ON_TechBias_aCGH$V1)),tech_bias=as.factor(ON_TechBias_aCGH$V2))

MZPBbias_aCGH <- as.data.frame(cbind(as.numeric(filt_hybe_v_ngs_MZPB$GCloess_weighted_median),"MZPBbias"))
MZILbias_aCGH <- as.data.frame(cbind(as.numeric(filt_hybe_v_ngs_MZIL$GCloess_weighted_median),"MZILbias"))
MZ_TechBias_aCGH <- rbind(MZPBbias_aCGH,MZILbias_aCGH)
MZ_TechBias_aCGH <- cbind.data.frame(log_ratio=as.numeric(as.character(MZ_TechBias_aCGH$V1)),tech_bias=as.factor(MZ_TechBias_aCGH$V2))

TechBias_aCGH_boxplot <- rbind(ON_TechBias_aCGH,MZ_TechBias_aCGH)
boxplot(log_ratio~tech_bias,data=TechBias_aCGH_boxplot,outline=FALSE)
abline(h = 0, col = "gray")
abline(h = mean(as.numeric(as.character(full_aCGH$V1)), na.rm = TRUE), col = "red")

dim(ONPBbias_aCGH)
dim(ONILbias_aCGH)
dim(MZPBbias_aCGH)
dim(MZILbias_aCGH)

## Check for differences in aCGH ratios among NGS-biased probes
TechBias_aCGH_aov <- rbind(full_aCGH,ONPBbias_aCGH,ONILbias_aCGH,MZPBbias_aCGH,MZILbias_aCGH)
TechBias_aCGH_aov <- cbind.data.frame(log_ratio=as.numeric(as.character(TechBias_aCGH_aov$V1)),tech_bias=as.factor(TechBias_aCGH_aov$V2))
anova_res <- aov(formula = log_ratio ~ tech_bias, data = TechBias_aCGH_aov)
summary(anova_res)
## Which groups of probes differ
tuk<- TukeyHSD(anova_res)
tuk

## Violin plot to look at aCGH ratios of NGS-biased probes
vioplot(log_ratio~tech_bias,data=TechBias_aCGH_boxplot)
abline(h = 0, col = "gray")
abline(h = mean(as.numeric(as.character(full_aCGH$V1)), na.rm = TRUE), col = "red")



###########################################################################################################################

                      ##########################################################################
                      ### DO NOT RUN THIS SECTION UNTIL YOU ARE DONE WITH ALL ABOVE ANALYSES ###
                      ##########################################################################

## The following code is revised from the above R code and reimports new tables that overwrite variables for reanalysis ##
##            The goal of this code is to create Supplementary Figure 1 associated with the manuscript                  ##

## import ndf file, hybe ratio file, and hits per species_seqtech files
ndf <- read.table(file = '120718_Cichlid_EB_CGH_filt.ndf', sep = '\t', header = TRUE)
hybe <- read.table(file = 'combined_hybe_ratios2.txt', sep = '\t', header = TRUE)
## These versions of the hit-count files have counts at various BLASTn thresholds from 90% to 100% identity
mzi <- read.table(file = 'hits_Mzeb_illumina5.txt', sep = '\t', header = FALSE)
mzp <- read.table(file = 'hits_Mzeb_pacbio5.txt', sep = '\t', header = FALSE)
oni <- read.table(file = 'hits_Oreo_illumina5.txt', sep = '\t', header = FALSE)
onp <- read.table(file = 'hits_Oreo_pacbio5.txt', sep = '\t', header = FALSE)

## merge hits per species_seqtech files
mzeb_hits <- merge(mzi, mzp, all=TRUE, by=c("V1","V3"))
onil_hits <- merge(oni, onp, all=TRUE, by=c("V1","V3"))

## make final hits per species_seqtech files
ngs_hits <- merge (mzeb_hits, onil_hits, all=TRUE, by=c("V1","V3"))

## merge ndf and hybe ratio file
ndf_hybe <- merge(ndf, hybe, all=TRUE, by=c("CHROM", "POSITION"))

## merge previous ndf/hybe ratio file and hit
hybe_v_ngs <- merge (ndf_hybe, ngs_hits, all=TRUE, by.x="PROBE_DESIGN_ID", by.y="V1")

## remove rownames for potential downstream export back to command line
rownames(hybe_v_ngs) <- c()

## exported table as a test
#exptest <- write.table(hybe_v_ngs, file="exptest.tsv", quote=FALSE, sep='\t')

## set NA to 0 in hits per species_seqtech columns
#hybe_v_ngs$HITS_MZEB_ILLUMINA[is.na(hybe_v_ngs$HITS_MZEB_ILLUMINA)] <- 0
#hybe_v_ngs$HITS_MZEB_PACBIO[is.na(hybe_v_ngs$HITS_MZEB_PACBIO)] <- 0
#hybe_v_ngs$HITS_OREO_ILLUMINA[is.na(hybe_v_ngs$HITS_OREO_ILLUMINA)] <- 0
#hybe_v_ngs$HITS_OREO_PACBIO[is.na(hybe_v_ngs$HITS_OREO_PACBIO)] <- 0

## renaming columns
names(hybe_v_ngs)[names(hybe_v_ngs) == "X545383A05_CichlCGH02.WS1_635_qspline"] <- "GCloess_X545383A05_CichlCGH02.WS1_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X545385A05_CichlCGH03.wa05_635_qspline"] <- "GCloess_X545385A05_CichlCGH03.wa05_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X545409A05_CichlCGH05.wa04_635_qspline"] <- "GCloess_X545409A05_CichlCGH05.wa04_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X550971A05_CichlCGH09.wa01_635_qspline"] <- "GCloess_X550971A05_CichlCGH09.wa01_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X545383A05_CichlCGH02.WS1_635_qspline.1"] <- "MA2C_X545383A05_CichlCGH02.WS1_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X545385A05_CichlCGH03.wa05_635_qspline.1"] <- "MA2C_X545385A05_CichlCGH03.wa05_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X545409A05_CichlCGH05.wa04_635_qspline.1"] <- "MA2C_X545409A05_CichlCGH05.wa04_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "X550971A05_CichlCGH09.wa01_635_qspline.1"] <- "MA2C_X550971A05_CichlCGH09.wa01_635_qspline"
names(hybe_v_ngs)[names(hybe_v_ngs) == "chrom"] <- "chrom(dup)"
names(hybe_v_ngs)[names(hybe_v_ngs) == "pos"] <- "pos(dup)"

names(hybe_v_ngs)[names(hybe_v_ngs) == "V2.x.x"] <- "HITS_MZEB_ILLUMINA"
names(hybe_v_ngs)[names(hybe_v_ngs) == "V2.y.x"] <- "HITS_MZEB_PACBIO"
names(hybe_v_ngs)[names(hybe_v_ngs) == "V2.x.y"] <- "HITS_OREO_ILLUMINA"
names(hybe_v_ngs)[names(hybe_v_ngs) == "V2.y.y"] <- "HITS_OREO_PACBIO"

## Calculate weighted median hybe ratios across samples, using replicate samples as half the weight as singlicate
hybe_v_ngs$GCloess_rep_mean <- apply(hybe_v_ngs[,c("GCloess_X545383A05_CichlCGH02.WS1_635_qspline", "GCloess_X550971A05_CichlCGH09.wa01_635_qspline")],1,mean, na.rm = TRUE)
hybe_v_ngs$GCloess_weighted_median <- apply(hybe_v_ngs[,c("GCloess_rep_mean", "GCloess_X545385A05_CichlCGH03.wa05_635_qspline", "GCloess_X545409A05_CichlCGH05.wa04_635_qspline")],1,median, na.rm = TRUE)
### Meta-analysis lead us to use GCloess for primary normalization method in this study ###

## adding columns with Oreo/Mzeb log2 ratios for comparisons between sequencing techs 
## this means log2(HITS_OREO_ILLUMINA/HITS_MZEB_ILLUMINA) and the same for pacbio
hybe_v_ngs$log2_OREOMZEB_ILLUMINA <- log2(hybe_v_ngs$HITS_OREO_ILLUMINA/hybe_v_ngs$HITS_MZEB_ILLUMINA)
hybe_v_ngs$log2_OREOMZEB_PACBIO <- log2(hybe_v_ngs$HITS_OREO_PACBIO/hybe_v_ngs$HITS_MZEB_PACBIO)
## honestly not sure why I added these columns but maybe they'll be useful
hybe_v_ngs$log2_OREO_PACBIOILLUMINA <- log2(hybe_v_ngs$HITS_OREO_PACBIO/hybe_v_ngs$HITS_OREO_ILLUMINA)
hybe_v_ngs$log2_MZEB_PACBIOILLUMINA <- log2(hybe_v_ngs$HITS_MZEB_PACBIO/hybe_v_ngs$HITS_MZEB_ILLUMINA)


## filtering masterfile such that each probe will have at least one hit in each species_technology column
filt_hybe_v_ngs <- hybe_v_ngs[is.na(hybe_v_ngs$HITS_MZEB_ILLUMINA) == FALSE 
                              & is.na(hybe_v_ngs$HITS_MZEB_PACBIO) == FALSE 
                              & is.na(hybe_v_ngs$HITS_OREO_ILLUMINA) == FALSE 
                              & is.na(hybe_v_ngs$HITS_OREO_PACBIO) == FALSE, ]

## Select Probes that have multiple copies in at least one assembly
filt_hybe_v_ngs_v2 <- filt_hybe_v_ngs # copy master file to separate file and keep original as backup in case there are any formatting issues in future steps
filt_hybe_v_ngs_v2$dup_check <- duplicated(filt_hybe_v_ngs_v2[,c('PROBE_ID','V3')])
filt_hybe_v_ngs_v3 <- filt_hybe_v_ngs_v2[!(filt_hybe_v_ngs_v2$dup_check=="TRUE"),] 
dim(filt_hybe_v_ngs_v2) # check the number of rows before duplicate filtering
dim(filt_hybe_v_ngs_v3) # check the number of rows after duplicate filtering
## Should be the same number of rows, in which case it's fine to use the v2 masterfile for the following filter steps

## Select Probes that have multiple copies in any assembly
filt_multihit_hybe_v_ngs <- filter(filt_hybe_v_ngs_v2, (HITS_MZEB_PACBIO * HITS_MZEB_ILLUMINA * HITS_OREO_PACBIO * HITS_OREO_ILLUMINA) > 1)
## Select Probes that have multiple copies in at least one PacBio assembly
filt_multihitPB_hybe_v_ngs <- filter(filt_hybe_v_ngs_v2, (HITS_MZEB_PACBIO * HITS_OREO_PACBIO) > 1)
## Select Probes that have multiple copies in at least one Illumina assembly
filt_multihitIL_hybe_v_ngs <- filter(filt_hybe_v_ngs_v2, (HITS_MZEB_ILLUMINA * HITS_OREO_ILLUMINA) > 1)
## Select Probes that have at least one 100% hit in all assemblies
filt_hybe_v_ngs_v4 <- filt_hybe_v_ngs_v2[(filt_hybe_v_ngs_v2$V3==100),] 
filt_100hit_hybe_v_ngs <- filter(filt_hybe_v_ngs_v4, (HITS_MZEB_PACBIO * HITS_MZEB_ILLUMINA * HITS_OREO_PACBIO * HITS_OREO_ILLUMINA) > 0)

blasthit_list <- filt_100hit_hybe_v_ngs$PROBE_DESIGN_ID
filt_100hit_hybe_v_ngs_v2 <- subset(filt_hybe_v_ngs_v2, PROBE_DESIGN_ID %in% blasthit_list )

## Apply filters for running platform and species comparisons
## Use +/-0.3 and +/-0.8 filters for hybe ratios
filt_hybe_v_ngs_filt1 <- filt_100hit_hybe_v_ngs_v2
filt_hybe_v_ngs_filt2 <- filt_100hit_hybe_v_ngs_v2
filt_hybe_v_ngs_filt1[(filt_hybe_v_ngs_filt1 < 0.8) & (filt_hybe_v_ngs_filt1 > -0.8)] <- NA
filt_hybe_v_ngs_filt2[(filt_hybe_v_ngs_filt2 < 0.3) & (filt_hybe_v_ngs_filt2 > -0.3)] <- NA

filt_hybe_v_ngs_ONPB <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$log2_OREO_PACBIOILLUMINA > 0),] #PacBio biased in Onil
filt_hybe_v_ngs_ONIL <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$log2_OREO_PACBIOILLUMINA < 0),] #Illumina biased in Onil
filt_hybe_v_ngs_MZPB <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$log2_MZEB_PACBIOILLUMINA > 0),] #PacBio biased in Mzeb
filt_hybe_v_ngs_MZIL <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$log2_MZEB_PACBIOILLUMINA < 0),] #Illumina biased in Mzeb

filt_hybe_v_ngs_ONPB <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$HITS_OREO_PACBIO > 1),] #PacBio biased in Onil
filt_hybe_v_ngs_ONIL <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$HITS_OREO_ILLUMINA > 1),] #Illumina biased in Onil
filt_hybe_v_ngs_MZPB <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$HITS_MZEB_PACBIO > 1),] #PacBio biased in Mzeb
filt_hybe_v_ngs_MZIL <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$HITS_MZEB_ILLUMINA > 1),] #Illumina biased in Mzeb

filt_hybe_v_ngs_ON_bias_PB <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0),] #Onil biased in PacBio
filt_hybe_v_ngs_ON_bias_IL <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0),] #Onil biased in Illumina
filt_hybe_v_ngs_MZ_bias_PB <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < 0),] #Mzeb biased in PacBio
filt_hybe_v_ngs_MZ_bias_IL <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < 0),] #Mzeb biased in Illumina
filt_hybe_v_ngs_ON_bias_aCGH <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$GCloess_weighted_median > 0),] #Onil biased in aCGH
filt_hybe_v_ngs_MZ_bias_aCGH <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$GCloess_weighted_median < 0),] #Mzeb biased in aCGH

filt_hybe_v_ngs_ON_bias_aCGH_PB <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$GCloess_weighted_median > 0) & (filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0),] #Onil biased in aCGH and PacBio
filt_hybe_v_ngs_ON_bias_aCGH_IL <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$GCloess_weighted_median > 0) & (filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0),] #Onil biased in aCGH and Illumina
filt_hybe_v_ngs_MZ_bias_aCGH_PB <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$GCloess_weighted_median < 0) & (filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < 0),] #Mzeb biased in aCGH and PacBio
filt_hybe_v_ngs_MZ_bias_aCGH_IL <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$GCloess_weighted_median < 0) & (filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < 0),] #Mzeb biased in aCGH and Illumina

filt_hybe_v_ngs_ON_gain_PB <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0.8),] #Onil gain in PacBio (past 0.8 threshold)
filt_hybe_v_ngs_ON_gain_IL <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0.8),] #Onil gain in Illumina (past 0.8 threshold)
filt_hybe_v_ngs_MZ_gain_PB <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < -0.8),] #Mzeb gain in PacBio (past 0.8 threshold)
filt_hybe_v_ngs_MZ_gain_IL <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < -0.8),] #Mzeb gain in Illumina (past 0.8 threshold)
filt_hybe_v_ngs_ON_gain_aCGH <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$GCloess_weighted_median > 0.8),] #Onil gain in aCGH (past 0.8 threshold)
filt_hybe_v_ngs_MZ_gain_aCGH <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$GCloess_weighted_median < -0.8),] #Mzeb gain in aCGH (past 0.8 threshold)

filt_hybe_v_ngs_ON_gain_aCGH_PB <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$GCloess_weighted_median > 0.8) & (filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO > 0.8),] #Onil gain in aCGH and PacBio (past 0.8 threshold)
filt_hybe_v_ngs_ON_gain_aCGH_IL <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$GCloess_weighted_median > 0.8) & (filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA > 0.8),] #Onil gain in aCGH and Illumina (past 0.8 threshold)
filt_hybe_v_ngs_MZ_gain_aCGH_PB <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$GCloess_weighted_median < -0.8) & (filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_PACBIO < -0.8),] #Mzeb gain in aCGH and PacBio (past 0.8 threshold)
filt_hybe_v_ngs_MZ_gain_aCGH_IL <- filt_100hit_hybe_v_ngs_v2[(filt_100hit_hybe_v_ngs_v2$GCloess_weighted_median < -0.8) & (filt_100hit_hybe_v_ngs_v2$log2_OREOMZEB_ILLUMINA < -0.8),] #Mzeb gain in aCGH and Illumina (past 0.8 threshold)


## Initialize the empty files
MZ_threshold_IL <- data.frame()
MZ_threshold_PB <- data.frame()
ON_threshold_IL <- data.frame()
ON_threshold_PB <- data.frame()

##################################################################################################
###  Time to calculate reciprocal overlap in Euler plots for different assembly blast filters  ###
##################################################################################################
for(var in 90:100)
{
  filt_100hit_hybe_v_ngs_v3 <- filter(filt_100hit_hybe_v_ngs_v2, V3 == var)
  
  ## export probes from ON gains identified by all 3 technologies and their overlap
  ON_gains_aCGH <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$GCloess_weighted_median > 0.8),] # Number of probes as ON gains in aCGH only
  
  #ON_gains_aCGH <- filt_100hit_hybe_v_ngs_v3_bk[(filt_100hit_hybe_v_ngs_v3_bk$GCloess_weighted_median > 0.8),] # Number of probes as ON gains in aCGH only
  
  ON_gains_PB <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_PACBIO > 0.8),] # Number of probes as ON gains PacBio only
  ON_gains_IL <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_ILLUMINA > 0.8),] # Number of probes as ON gains Illumina only
  ON_gains_aCGH_PB <- filt_hybe_v_ngs_ON_gain_aCGH_PB # Number of probes as ON gains in aCGH and PacBio
  ON_gains_aCGH_IL <- filt_hybe_v_ngs_ON_gain_aCGH_IL # Number of probes as ON gains in aCGH and Illumina
  ON_gains_PB_IL <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_ILLUMINA > 0.8) & (filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_PACBIO > 0.8),] # Number of probes as ON gains in PacBio and Illumina
  ON_gains_aCGH_PB_IL <- filt_hybe_v_ngs_ON_gain_aCGH_IL[(filt_hybe_v_ngs_ON_gain_aCGH_IL$log2_OREOMZEB_PACBIO > 0.8),] # Number of probes as ON gains in all 3 technologies
  
  ## export probes from ON bias identified by all 3 technologies and their overlap
  ON_bias_aCGH <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$GCloess_weighted_median > 0),] # Number of probes as ON bias in aCGH only
  ON_bias_PB <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_PACBIO > 0),] # Number of probes as ON bias PacBio only
  ON_bias_IL <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_ILLUMINA > 0),] # Number of probes as ON bias Illumina only
  ON_bias_aCGH_PB <- filt_hybe_v_ngs_ON_bias_aCGH_PB # Number of probes as ON bias in aCGH and PacBio
  ON_bias_aCGH_IL <- filt_hybe_v_ngs_ON_bias_aCGH_IL # Number of probes as ON bias in aCGH and Illumina
  ON_bias_PB_IL <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_ILLUMINA > 0) & (filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_PACBIO > 0),] # Number of probes as ON bias in PacBio and Illumina
  ON_bias_aCGH_PB_IL <- filt_hybe_v_ngs_ON_gain_aCGH_IL[(filt_hybe_v_ngs_ON_gain_aCGH_IL$log2_OREOMZEB_PACBIO > 0),] # Number of probes as ON bias in all 3 technologies
  
  #export probes from MZ gains identified by all 3 technologies and their overlap
  MZ_gains_aCGH <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$GCloess_weighted_median < -0.8),] # Number of probes as MZ gains in aCGH only
  
  #MZ_gains_aCGH <- filt_100hit_hybe_v_ngs_v3_bk[(filt_100hit_hybe_v_ngs_v3_bk$GCloess_weighted_median < -0.8),] # Number of probes as MZ gains in aCGH only
  
  MZ_gains_PB <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_PACBIO < -0.8),]  # Number of probes as MZ gains in PacBio only
  MZ_gains_IL <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_ILLUMINA < -0.8),] # Number of probes as MZ gains in Illumina only
  MZ_gains_aCGH_PB <- filt_hybe_v_ngs_MZ_gain_aCGH_PB # Number of probes as MZ gains in aCGH and PacBio
  MZ_gains_aCGH_IL <- filt_hybe_v_ngs_MZ_gain_aCGH_IL # Number of probes as MZ gains in aCGH and Illumina
  MZ_gains_PB_IL <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_ILLUMINA < -0.8) & (filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_PACBIO < -0.8),] # Number of probes as MZ gains in PacBio and Illumina
  MZ_gains_aCGH_PB_IL <- filt_hybe_v_ngs_MZ_gain_aCGH_IL[(filt_hybe_v_ngs_MZ_gain_aCGH_IL$log2_OREOMZEB_PACBIO < -0.8),] # Number of probes as MZ gains in all 3 technologies
  
  #export probes from MZ bias identified by all 3 technologies and their overlap
  MZ_bias_aCGH <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$GCloess_weighted_median < 0),] # Number of probes as MZ bias in aCGH only
  MZ_bias_PB <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_PACBIO < 0),]  # Number of probes as MZ bias in PacBio only
  MZ_bias_IL <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_ILLUMINA < 0),] # Number of probes as MZ bias in Illumina only
  MZ_bias_aCGH_PB <- filt_hybe_v_ngs_MZ_bias_aCGH_PB # Number of probes as MZ bias in aCGH and PacBio
  MZ_bias_aCGH_IL <- filt_hybe_v_ngs_MZ_bias_aCGH_IL # Number of probes as MZ bias in aCGH and Illumina
  MZ_bias_PB_IL <- filt_100hit_hybe_v_ngs_v3[(filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_ILLUMINA < 0) & (filt_100hit_hybe_v_ngs_v3$log2_OREOMZEB_PACBIO < 0),] # Number of probes as MZ bias in PacBio and Illumina
  MZ_bias_aCGH_PB_IL <- filt_hybe_v_ngs_MZ_gain_aCGH_IL[(filt_hybe_v_ngs_MZ_gain_aCGH_IL$log2_OREOMZEB_PACBIO < 0),] # Number of probes as MZ bias in all 3 technologies
  
  ## First, find overlap in hypothetical euler plots for probes
  MZ_bias_PB_list <- unique(MZ_bias_PB$PROBE_DESIGN_ID)
  MZ_bias_IL_list <- unique(MZ_bias_IL$PROBE_DESIGN_ID)
  MZ_gains_aCGH_list <- MZ_gains_aCGH$PROBE_DESIGN_ID
  MZ_gains_aCGH_list <- unique(MZ_gains_aCGH_list[!is.na(MZ_gains_aCGH_list)])
  
  MZ_bias_PB_genic <- filter(MZ_bias_PB, CONTAINER != "BACKBONE")
  MZ_bias_IL_genic <- filter(MZ_bias_IL, CONTAINER != "BACKBONE")
  MZ_gains_aCGH_genic <- filter(MZ_gains_aCGH, CONTAINER != "BACKBONE")
  
  MZ_bias_PB_list_genic <- MZ_bias_PB_genic$PROBE_DESIGN_ID
  MZ_bias_IL_list_genic <- MZ_bias_IL_genic$PROBE_DESIGN_ID
  MZ_gains_aCGH_list_genic <- MZ_gains_aCGH_genic$PROBE_DESIGN_ID
  MZ_gains_aCGH_list_genic <- MZ_gains_aCGH_list_genic[!is.na(MZ_gains_aCGH_list_genic)]
  
  ON_bias_PB_list <- unique(ON_bias_PB$PROBE_DESIGN_ID)
  ON_bias_IL_list <- unique(ON_bias_IL$PROBE_DESIGN_ID)
  ON_gains_aCGH_list <- ON_gains_aCGH$PROBE_DESIGN_ID
  ON_gains_aCGH_list <- unique(ON_gains_aCGH_list[!is.na(ON_gains_aCGH_list)])
  
  ON_bias_PB_genic <- filter(ON_bias_PB, CONTAINER != "BACKBONE")
  ON_bias_IL_genic <- filter(ON_bias_IL, CONTAINER != "BACKBONE")
  ON_gains_aCGH_genic <- filter(ON_gains_aCGH, CONTAINER != "BACKBONE")
  
  ON_bias_PB_list_genic <- ON_bias_PB_genic$PROBE_DESIGN_ID
  ON_bias_IL_list_genic <- ON_bias_IL_genic$PROBE_DESIGN_ID
  ON_gains_aCGH_list_genic <- ON_gains_aCGH_genic$PROBE_DESIGN_ID
  ON_gains_aCGH_list_genic <- ON_gains_aCGH_list_genic[!is.na(ON_gains_aCGH_list_genic)]
  
  ## Build tables with counts in each section of hypothetical 2-way Euler plots
  MZ_threshold_IL <- rbind(MZ_threshold_IL,cbind(var, length(MZ_gains_aCGH_list),length(MZ_bias_IL_list),length(intersect(MZ_gains_aCGH_list,MZ_bias_IL_list))))
  MZ_threshold_PB <- rbind(MZ_threshold_PB,cbind(var, length(MZ_gains_aCGH_list),length(MZ_bias_PB_list),length(intersect(MZ_gains_aCGH_list,MZ_bias_PB_list))))
  ON_threshold_IL <- rbind(ON_threshold_IL,cbind(var, length(ON_gains_aCGH_list),length(ON_bias_IL_list),length(intersect(ON_gains_aCGH_list,ON_bias_IL_list))))
  ON_threshold_PB <- rbind(ON_threshold_PB,cbind(var, length(ON_gains_aCGH_list),length(ON_bias_PB_list),length(intersect(ON_gains_aCGH_list,ON_bias_PB_list))))
  
}

### Combine results from each pairwise Euler comparison and calculate proportions for plotting
##  These table values reflect the portion of probes identified as gains in by the first platform overlapping with the list of gains from the second platform
full_overlap_table <- data.frame(cbind(thresh = MZ_threshold_IL$var, MZ_aCGH_IL_overlap = MZ_threshold_IL$V4/MZ_threshold_IL$V2, MZ_IL_aCGH_overlap = MZ_threshold_IL$V4/MZ_threshold_IL$V3, MZ_aCGH_PB_overlap = MZ_threshold_PB$V4/MZ_threshold_PB$V2, MZ_PB_aCGH_overlap = MZ_threshold_PB$V4/MZ_threshold_PB$V3, ON_aCGH_IL_overlap = ON_threshold_IL$V4/ON_threshold_IL$V2, ON_IL_aCGH_overlap = ON_threshold_IL$V4/ON_threshold_IL$V3, ON_aCGH_PB_overlap = ON_threshold_PB$V4/ON_threshold_PB$V2, ON_PB_aCGH_overlap = ON_threshold_PB$V4/ON_threshold_PB$V3))

## Plot overlaps within MZ
MZ_SI_overlap <- ggplot(full_overlap_table, aes(x=thresh)) + 
  geom_line(aes(y = MZ_aCGH_IL_overlap), color = "darkred") + 
  geom_line(aes(y = MZ_IL_aCGH_overlap), color= "darkred", linetype="twodash") + 
  geom_line(aes(y = MZ_aCGH_PB_overlap), color = "steelblue") + 
  geom_line(aes(y = MZ_PB_aCGH_overlap), color= "steelblue", linetype="twodash") +
  ylim(0, 1) +
  scale_x_continuous(breaks = seq(90, 100, by = 1))

## Plot overlaps within ON
ON_SI_overlap <- ggplot(full_overlap_table, aes(x=thresh)) + 
  geom_line(aes(y = ON_aCGH_IL_overlap), color = "darkred") + 
  geom_line(aes(y = ON_IL_aCGH_overlap), color= "darkred", linetype="twodash") + 
  geom_line(aes(y = ON_aCGH_PB_overlap), color = "steelblue") + 
  geom_line(aes(y = ON_PB_aCGH_overlap), color= "steelblue", linetype="twodash") +
  ylim(0, 1) +
  scale_x_continuous(breaks = seq(90, 100, by = 1))

## Finally plot them together
grid.arrange(MZ_SI_overlap, ON_SI_overlap)
