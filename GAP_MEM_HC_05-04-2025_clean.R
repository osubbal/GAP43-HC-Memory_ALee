#Annie Lee
#ADNI GAP-43, HC, Memory Project
#Updated: 11/12/2024
#Cleaned, uploading to OneDrive

#load in packages
library(dplyr)
library(polycor)
library(psych)
library(ggplot2)
library(emmeans)

#set directory
#setwd("/Users/annielee/OneDrive/ADNI")
setwd("/Users/oliviahorn/Library/CloudStorage/OneDrive-TheOhioStateUniversity/BBAL/Projects/BADS_2021E0199/ADNI/ADNI_Annie_GAP43-HC-Mem/working")
#Read in datasets
adnimerge <- read.csv("data/ADNIMERGE.csv") #participant info and baseline HC
uwn <- read.csv("data/UWNPSYCHSUM_03_09_21.csv") #neuropsych composite scores
gap43 <- read.csv("data/BLENNOW_LAB_CSF_GAP_43_06_08_21.csv",header=T) #GAP-43
upenn <- read.csv("data/UPENNBIOMK9_04_19_17.csv") #CSF biomarkers
pet <- read.csv("data/UCBERKELEYFDG_8mm_02_17_23_22May2023.csv")
ucsf <- read.csv("data/UCSFFSX51_11_08_19_15Mar2024.csv") #has both accelerated and non-accelerated t1, https://adni.bitbucket.io/reference/ucsffsx51.html

#-----------------------------------------------------------------------------#
#-----------------------CLEANING AND MERGING FILES----------------------------#
#-----------------------------------------------------------------------------#

#Filter all datasets
adnimerge_filter <- adnimerge %>% filter(VISCODE=='bl')
gap43_filter <- gap43 %>% filter(VISCODE2=='bl')
upenn_filter <- upenn %>% filter(VISCODE2=='bl')
pet_filter <- pet %>% filter(VISCODE2=='bl')
ucsf_filter <- ucsf %>% filter(VISCODE2=='scmri', IMAGETYPE=='Non-Accelerated T1')

#Subset all column data that we want
adnimerge_filter <- adnimerge_filter %>% dplyr::select(RID,VISCODE,DX_bl,AGE,PTGENDER,PTEDUCAT,MMSE,APOE4,ADAS13,Hippocampus_bl, ICV_bl)
gap43_filter <- gap43_filter %>% dplyr::select(RID, VISCODE2, GAP_43)
uwn_filter <- uwn %>% dplyr::select(RID, VISCODE2, ADNI_MEM, ADNI_EF)
upenn_filter <- upenn_filter %>% dplyr::select(RID,VISCODE2,ABETA,PTAU,TAU)
pet_filter <- pet_filter %>% dplyr::select(RID,ROINAME,MEAN,VISCODE2)

#Change 'scmri' to 'bl' to merge
names(adnimerge_filter)[names(adnimerge_filter) == 'VISCODE'] <- "VISCODE2"
ucsf_filter$VISCODE2[ucsf_filter$VISCODE2=='scmri'] <- 'bl'   

#Clean PET dataset for AT(N)
pet_filter <- aggregate(MEAN ~RID,pet_filter, function(x) x[1]/x[-1])
pet_filter$VISCODE2 <- "bl"
names(pet_filter)[names(pet_filter) == 'MEAN'] <- "metaROI"

#Create baseline database with HC subfields ("all_data_hc")
all_data_hc <- Reduce(function(x,y) merge(x,y, by=c("RID","VISCODE2")), list(adnimerge_filter,uwn_filter,gap43_filter,upenn_filter, pet_filter, ucsf_filter))
length(all_data_hc$RID) #N=789


#----------------------- Clean AD biomarker data ------------------------------#
#1. Truncate biomarkers
#Total Tau: 80 ≤ TAU ≤ 1300 pg/ML
all_data_hc$TAU[as.numeric(all_data_hc$TAU) > 1300] <- 1300
all_data_hc$TAU[as.numeric(all_data_hc$TAU) < 80] <- 80
all_data_hc$TAU <- as.numeric(all_data_hc$TAU)

# PTAU: 8 ≤ PTAU ≤ 120
all_data_hc$PTAU[all_data_hc$PTAU=="<8"] <- 7
all_data_hc$PTAU[as.numeric(all_data_hc$PTAU) > 120] <- 120
all_data_hc$PTAU[as.numeric(all_data_hc$PTAU) < 8] <- 8
all_data_hc$PTAU <- as.numeric(all_data_hc$PTAU)

# ABETA: 200 ≤ ABETA ≤ 1700
all_data_hc$ABETA[all_data_hc$ABETA==">1700"] <- 1800
all_data_hc$ABETA[as.numeric(all_data_hc$ABETA) > 1700] <- 1700
all_data_hc$ABETA[as.numeric(all_data_hc$ABETA) < 200] <- 200
all_data_hc$ABETA <- as.numeric(all_data_hc$ABETA)

#2. Binarize AD biomarkers to control for AT(N) group
#Qing 2022 Reference:
#A = Pathological ABeta1-42 defined as < 977 pg/ml
#T = Pathological P-Tau181 defined as > 27 pg/ml
#N = Pathological FDG-PET defined as < 1.21 
# 1 = pathological ; 0 = not pathological

# ABETA binarize (0 = non-pathological Abeta ; 1 = pathological Abeta)
all_data_hc$ABETA <- as.numeric(all_data_hc$ABETA)
for (i in 1:length(all_data_hc$ABETA)) {
  if (all_data_hc$ABETA[i] < 977) {
    (all_data_hc$ABETA_atn[i] <- 1)
  } else if (all_data_hc$ABETA[i] >= 977) {
    (all_data_hc$ABETA_atn[i] <- 0)
  }
}

table(all_data_hc$DX_bl,all_data_hc$ABETA_atn)
table(all_data_hc$ABETA_atn)

# P-TAU binarize
all_data_hc$PTAU <- as.numeric(all_data_hc$PTAU)
for (i in 1:length(all_data_hc$PTAU)) {
  if (all_data_hc$PTAU[i] > 27) {
    (all_data_hc$PTAU_atn[i] <- 1)
  } else if (all_data_hc$PTAU[i] <= 27) {
    (all_data_hc$PTAU_atn[i] <- 0)
  }
}

table(all_data_hc$DX_bl,all_data_hc$PTAU_atn)
table(all_data_hc$PTAU_atn) 

# NEURODEGENERATION
all_data_hc$metaROI<- as.numeric(all_data_hc$metaROI)
for (i in 1:length(all_data_hc$metaROI)) {
  if (all_data_hc$metaROI[i] <= 1.21) {
    (all_data_hc$PET_atn[i] <- 1)
  } else if (all_data_hc$metaROI[i] > 1.21) {
    (all_data_hc$PET_atn[i] <- 0)
  }
}

table(all_data_hc$PET_atn) #0=309, 1=254


#2. Create ATN classification strings - "ATN_class" [A(+/-), T(+/-), N(+/-)]
# [A+/-] for ABETA
for (i in 1:length(all_data_hc$ABETA_atn)) {
  if (all_data_hc$ABETA_atn[i] == 1) {
    (all_data_hc$ATN_class[i] <- "A+")
  } else if (all_data_hc$ABETA_atn[i] == 0) {
    (all_data_hc$ATN_class[i] <- "A-")
  }
}

for (i in 1:length(all_data_hc$PTAU_atn)) {
  if (all_data_hc$PTAU_atn[i] == 1) {
    (all_data_hc$ATN_class[i] <- paste(all_data_hc$ATN_class[i],"T+"))
  } else if (all_data_hc$PTAU_atn[i] == 0) {
    (all_data_hc$ATN_class[i] <- paste(all_data_hc$ATN_class[i],"T-"))
  }
}

for (i in 1:length(all_data_hc$PET_atn)) {
  if (all_data_hc$PET_atn[i] == 1) {
    (all_data_hc$ATN_class[i] <- paste(all_data_hc$ATN_class[i],"N+"))
  } else if (all_data_hc$PET_atn[i] == 0) {
    (all_data_hc$ATN_class[i] <- paste(all_data_hc$ATN_class[i],"N-"))
  }
}

all_data_hc$ATN_class <- as.factor(all_data_hc$ATN_class)

#3. Create ATN classification groups groups -- "ATN_group" [CN, SNAP, AD]
for (i in 1:length(all_data_hc$ATN_class)) {
  if (all_data_hc$ATN_class[i] == "A- T- N-") {
    (all_data_hc$ATN_group[i] <- "CN")
  } else if (all_data_hc$ATN_class[i] == "A- T+ N+" | all_data_hc$ATN_class[i] == "A- T- N+" | all_data_hc$ATN_class[i] == "A- T+ N-") {
    (all_data_hc$ATN_group[i] <- "SNAP")
  }
  else if (all_data_hc$ATN_class[i] == "A+ T- N-" | all_data_hc$ATN_class[i] =="A+ T+ N-" | all_data_hc$ATN_class[i] == "A+ T+ N+" | all_data_hc$ATN_class[i] == "A+ T- N+" )
    (all_data_hc$ATN_group[i] <- "AD")
}

all_data_hc$ATN_group <- as.factor(all_data_hc$ATN_group)
table(all_data_hc$ATN_group) #AD n=419, SNAP n=123, CN n=247

#--------------------Create hippocampal subfields ROIs------------------------#
#Rename ROIs to CA1, Presubiculum, Subiculum, Dentate Gyrus (CA4)
names(all_data_hc)[names(all_data_hc) == 'ST131HS'] <- "CA1_Left"
names(all_data_hc)[names(all_data_hc) == 'ST139HS'] <- "CA1_Right"
names(all_data_hc)[names(all_data_hc) == 'ST136HS'] <- "Presubiculum_Left"
names(all_data_hc)[names(all_data_hc) == 'ST144HS'] <- "Presubiculum_Right"
names(all_data_hc)[names(all_data_hc) == 'ST137HS'] <- "Subiculum_Left"
names(all_data_hc)[names(all_data_hc) == 'ST145HS'] <- "Subiculum_Right"
names(all_data_hc)[names(all_data_hc) == 'ST141HS'] <- "DG_Right"
names(all_data_hc)[names(all_data_hc) == 'ST133HS'] <- "DG_Left"

#Create average (noted below as 'sum') of R+L subfields
all_data_hc$CA1_sum <- (all_data_hc$CA1_Left + all_data_hc$CA1_Right) / 2
all_data_hc$Presubiculum_sum <- (all_data_hc$Presubiculum_Left + all_data_hc$Presubiculum_Right) / 2
all_data_hc$Subiculum_sum <- (all_data_hc$Subiculum_Left + all_data_hc$Subiculum_Right) / 2
all_data_hc$DG_sum <- (all_data_hc$DG_Left + all_data_hc$DG_Right) / 2

#------------------------ Additional data cleaning ------------------------------------#
#Subset dataframe to select needed columns
all_data_hc <- all_data_hc %>% dplyr::select(RID, VISCODE2, DX_bl, AGE, PTGENDER, PTEDUCAT, MMSE, APOE4, Hippocampus_bl, ICV_bl,
                                             GAP_43, ABETA, PTAU, TAU, ADNI_MEM, ADNI_EF, CA1_sum, DG_sum, Subiculum_sum, Presubiculum_sum, ATN_group)

#Remove dupliates
all_data_hc[all_data_hc$RID %in% all_data_hc$RID[duplicated(all_data_hc$RID)],]
dups <- c(2238, 4037, 4164, 4209, 4271, 4308, 4376, 4387, 4427, 4741, 4765, 4809, 4877, 4910, 5005, 5023, 5023, 5037, 5269)

dups <- data.frame(dups)
all_data_hc <- all_data_hc[!(all_data_hc$RID %in% dups$dups),] #remove duplicates
length(all_data_hc$RID)

#Binarize APOE4 carrier status
all_data_hc$APOE4[all_data_hc$APOE4=="2"] <- 1

#Binarize PTGENDER
all_data_hc$PTGENDER <- ifelse(all_data_hc$PTGENDER=="Female", 1, 0) #Female=1
all_data_hc$PTGENDER <- as.factor(all_data_hc$PTGENDER)

#Missing value check
cbind(lapply(lapply(all_data_hc, is.na), sum)) #64 missing hippocampal volume

#Remove participants w/ missing HC values
all_data_hc <- all_data_hc[!is.na(all_data_hc$Hippocampus_bl), ]
length(all_data_hc$RID)

#***** N = 689  ******

#-------------- ICV correction for total/subfield hippocampal data -------------#
#Total HC stats
mean(all_data_hc$Hippocampus_bl)
sd(all_data_hc$Hippocampus_bl)

#ROI stats
mean(all_data_hc$CA1_sum)
mean(all_data_hc$Presubiculum_sum)
mean(all_data_hc$Subiculum_sum)
mean(all_data_hc$DG_sum)

#ICV correction for total HC data
all_data_hc$HC_ICV <- all_data_hc$Hippocampus_bl / all_data_hc$ICV_bl

#ICV correction for subfields
all_data_hc$CA1_sum_ICV <- all_data_hc$CA1_sum / all_data_hc$ICV_bl
all_data_hc$Presubiculum_sum_ICV <- all_data_hc$Presubiculum_sum / all_data_hc$ICV_bl
all_data_hc$Subiculum_sum_ICV <- all_data_hc$Subiculum_sum / all_data_hc$ICV_bl
all_data_hc$DG_sum_ICV <- all_data_hc$DG_sum / all_data_hc$ICV_bl

#Multiplying the ROIs if needed
# all_data_hc$Subiculum_sum_ICV_mul <- all_data_hc$Subiculum_sum_ICV * 100
# all_data_hc$Presubiculum_sum_ICV_mul <- all_data_hc$Presubiculum_sum_ICV * 100
# all_data_hc$CA1_sum_ICV_mul <- all_data_hc$CA1_sum_ICV * 100
# all_data_hc$HC_ICV_mul <- all_data_hc$HC_ICV * 100
# all_data_hc$DG_ICV_mul <- all_data_hc$DG_sum_ICV * 100

#---------------------------- Outlier Removal ---------------------------------#
scattOut <- function(VOI) {
  upper3 = mean(VOI,na.rm=TRUE) + 3*sd(VOI,na.rm=TRUE)
  lower3 = mean(VOI,na.rm=TRUE) - 3*sd(VOI,na.rm=TRUE)
  upper5 = mean(VOI,na.rm=TRUE) + 5*sd(VOI,na.rm=TRUE)
  lower5 = mean(VOI,na.rm=TRUE) - 5*sd(VOI,na.rm=TRUE)
  plot(VOI)
  abline(h=upper3,col="blue")
  abline(h=lower3,col="blue")
  abline(h=upper5,col="red")
  abline(h=lower5,col="red")
}

#Log transform GAP-43
all_data_hc$GAP43_log <- log(all_data_hc$GAP_43) #insert log GAP43 into all_data_hc_long

scattOut(all_data_hc$GAP43_log)
scattOut(all_data_hc$ADNI_MEM)
scattOut(all_data_hc$HC_ICV) #1
scattOut(all_data_hc$CA1_sum_ICV) #5
scattOut(all_data_hc$DG_sum_ICV) #1
scattOut(all_data_hc$Subiculum_sum_ICV)
scattOut(all_data_hc$Presubiculum_sum_ICV) #2

#Total HC volume
uppermax_hip <- mean(all_data_hc$HC_ICV) + 3*sd(all_data_hc$HC_ICV)
lowermax_hip <- mean(all_data_hc$HC_ICV) - 3*sd(all_data_hc$HC_ICV)
all_data_hc[which(all_data_hc$HC_ICV > uppermax_hip | all_data_hc$HC_ICV < lowermax_hip),] 
hc_outliers <- c(4784, 5241)

#CA1 volume
uppermax_ca1 <- mean(all_data_hc$CA1_sum_ICV) + 3*sd(all_data_hc$CA1_sum_ICV)
lowermax_ca1 <- mean(all_data_hc$CA1_sum_ICV) - 3*sd(all_data_hc$CA1_sum_ICV)
all_data_hc[which(all_data_hc$CA1_sum_ICV > uppermax_ca1 | all_data_hc$CA1_sum_ICV < lowermax_ca1),] 
ca1_outliers <- c(4071, 4226, 4755, 4784, 5197)

#DG volume
uppermax_dg <- mean(all_data_hc$DG_sum_ICV) + 3*sd(all_data_hc$DG_sum_ICV)
lowermax_dg <- mean(all_data_hc$DG_sum_ICV) - 3*sd(all_data_hc$DG_sum_ICV)
all_data_hc[which(all_data_hc$DG_sum_ICV > uppermax_dg | all_data_hc$DG_sum_ICV < lowermax_dg),] 
dg_outliers <- c(2278)

#Presubiculum
uppermax_pres <- mean(all_data_hc$Presubiculum_sum_ICV) + 3*sd(all_data_hc$Presubiculum_sum_ICV)
lowermax_pres <- mean(all_data_hc$Presubiculum_sum_ICV) - 3*sd(all_data_hc$Presubiculum_sum_ICV)
all_data_hc[which(all_data_hc$Presubiculum_sum_ICV > uppermax_pres | all_data_hc$Presubiculum_sum_ICV < lowermax_pres),] 
presub_outliers <- c(2234, 2245)

outliers <- as.numeric(c(hc_outliers, ca1_outliers, dg_outliers, presub_outliers))
outliers <- data.frame(outliers)
outliers <- outliers %>% distinct(outliers, .keep_all = TRUE)
length(outliers$outliers) #9 total outliers

all_data_hc <- all_data_hc[!(all_data_hc$RID %in% outliers$outliers),] #remove outliers
length(all_data_hc$RID)

#***** N = 680  ******

all_data_hc[all_data_hc$RID %in% all_data_hc$RID[duplicated(all_data_hc$RID)],]

#--------------------- Checking data distributions ----------------------------#
#checking continuous variables
qqnorm(all_data_hc$AGE)
qqnorm(all_data_hc$PTEDUCAT)
  hist(all_data_hc$PTEDUCAT) #some skew towards higher education
qqnorm(all_data_hc$ADNI_MEM)

#checking ROIs
qqnorm(all_data_hc$HC_ICV)
qqnorm(all_data_hc$CA1_sum_ICV)
qqnorm(all_data_hc$DG_sum_ICV)
qqnorm(all_data_hc$Presubiculum_sum_ICV)
qqnorm(all_data_hc$Subiculum_sum_ICV)

#checking biomarker
qqnorm(all_data_hc$GAP_43)
    shapiro.test(all_data_hc$GAP_43)


#--------------------------- Norm all data ------------------------------------#
#Z-score all demographic, biomarker data
all_data_hc$PTEDUCAT_norm <- scale(all_data_hc$PTEDUCAT)
all_data_hc$AGE_norm <- scale(all_data_hc$AGE)
all_data_hc$GAP43_log_norm <- scale(all_data_hc$GAP43_log)
all_data_hc$ABETA_norm <- scale(all_data_hc$ABETA)
all_data_hc$PTAU_norm <- scale(all_data_hc$PTAU)

length(unique(all_data_hc$RID)) #N=680, correct

#Z-score all HC volume data
all_data_hc$HC_ICV_norm <- scale(all_data_hc$HC_ICV)
all_data_hc$CA1_sum_ICV <- scale(all_data_hc$CA1_sum_ICV)
all_data_hc$Presubiculum_sum_ICV <- scale(all_data_hc$Presubiculum_sum_ICV)
all_data_hc$Subiculum_sum_ICV <- scale(all_data_hc$Subiculum_sum_ICV)
all_data_hc$DG_sum_ICV <- scale(all_data_hc$DG_sum_ICV)

#Z-score multiplied HC volumes
# all_data_hc$HC_ICV_norm_mul <- scale(all_data_hc$HC_ICV_mul)
# all_data_hc$Subiculum_sum_ICV_mul <- scale(all_data_hc$Subiculum_sum_ICV_mul)
# all_data_hc$Presubiculum_sum_ICV_mul <- scale(all_data_hc$Presubiculum_sum_ICV_mul)
# all_data_hc$CA1_sum_ICV_mul <- scale(all_data_hc$CA1_sum_ICV_mul)
# all_data_hc$DG_ICV_norm_mul <- scale(all_data_hc$DG_ICV_mul)


#----------------------------------------------------------------------------#
#------------------- Table 1. Participant Demographics ----------------------#
#----------------------------------------------------------------------------#
#demographics at baseline
#https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html 

#CREATING TABLE
library(table1)

rndr <- function(x, name, ...) {
  if (!is.numeric(x)) return(render.categorical.default(x))
  what <- switch(name,
                 AGE = "Mean (SD)",
                 PTEDUCAT="Mean (SD)",
                 # APOE4="Mean (SD)",
                 ABETA="Mean (SD)",
                 TAU="Mean (SD)",
                 PTAU="Mean (SD)",
                 MMSE="Mean (SD)",
                 ADNI_MEM="Mean (SD)",
                 GAP_43 = "Mean (SD)",
                 HC_ICV = "Mean (SD)",
                 CA1_sum_ICV = "Mean (SD)",
                 Subiculum_sum_ICV = "Mean (SD)",
                 Presubiculum_sum_ICV = "Mean(SD")
  parse.abbrev.render.code(c("", what))(x)
}

labels <- list(
  variables=list(AGE="Age (in years)",
                 PTGENDER="Gender",
                 PTEDUCAT="Education (in years)",
                 APOE4="APOE4 allele status",
                 ABETA="Amyloid Beta",
                 TAU="Tau",
                 PTAU="P-Tau",
                 MMSE="MMSE",
                 ADNI_MEM="Memory Composite Score",
                 GAP_43 = "CSF GAP-43",
                 HC_ICV = "Total hippocampal volume",
                 CA1_sum_ICV = "CA1 volume",
                 Subiculum_sum_ICV = "Subiculum volume",
                 Presubiculum_sum_ICV = "Presubiculum volume"))


strata <- all_data_hc 
strata$ATN_group <- factor(strata$ATN_group, levels = c("CN", "SNAP", "AD"))
strata$APOE4 <- factor(strata$APOE4)

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=5), c("", "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD))) }
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,sprintf("%d (%0.0f %%)", FREQ, PCT)))) }

# CA1_sum, Presubiculum_sum, Subiculum_sum, DG_sum, Hippocampus_bl
table1(~ AGE + PTGENDER + PTEDUCAT + APOE4 + MMSE + ADNI_MEM + ABETA + PTAU + GAP_43 + Hippocampus_bl + Subiculum_sum + Presubiculum_sum + DG_sum + CA1_sum | ATN_group,
  data = strata, render.continuous=my.render.cont, render.categorical=my.render.cat, topclass= "Rtable1-times")

all_data_hc2 <- all_data_hc %>% filter(DX_bl=="CN" | DX_bl=="SMC" | DX_bl=="MCI")

#------------------------------------------------------------------------------#
#-------------------------PRELIMINARY ANALYSES---------------------------------#
#------------------------------------------------------------------------------#
#Executive Function
#summary(lm(all_data_hc$ADNI_EF ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm*HC_ICV_norm, data=all_data_hc))
#Supplemental: GAP-43 does not predicting EF and there is no significant GAPxHC interaction on EF,
#suggesting that this relationship is specific to episodic memory performance

library(interactions)

#make sure CN (no pathology group) is the reference group
all_data_hc <- within(all_data_hc, ATN_group <- relevel(factor(ATN_group), ref="CN"))
contrasts(all_data_hc$ATN_group) 

#------------------------------------------------------------------------------#
#--------------------------- Analysis Aim 1------------------------------------#
#------------------------------------------------------------------------------#
#Aim 1. Main effect of GAP-43 on episodic memory performance

#Step 1. Covariates predicting memory
step1 <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data=all_data_hc)
summary(step1)

#Step 2. Main effect of GAP-43 predicting memory
step2 <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm, data=all_data_hc)
summary(step2) #significant main effect of gap43 -0.10333    0.02959  -3.493  0.00051 ***
anova(step1,step2) #sig
summary(step2)$r.squared - summary(step1)$r.squared #change in Rsq= 0.01209138


#------------------------------------------------------------------------------#
#-------------------------- Analysis Aim 2  -----------------------------------#
#------------------------------------------------------------------------------#

#Aim 2. HC volume as a moderator on the relationship between GAP-43 and episodic memory performance

#------------------------------------------------------------------------------#
#1. Total Hippocampal Volume
#Step 1. Covariates predicting memory [identical as aim 1]
step1_hc <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data=all_data_hc)
summary(step1_hc)

#Step 2. Main effect of GAP-43 predicting memory [identical as aim 1]
step2_hc <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm, data=all_data_hc)
summary(step2_hc)
anova(step1_hc, step2_hc)
summary(step2_hc)$r.squared - summary(step1_hc)$r.squared #change in Rsq= 0.01209138

#Step 3. Main effect of HC
step3_hc <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm + HC_ICV_norm, data=all_data_hc)
summary(step3_hc) #significant, p<.001
anova(step2_hc, step3_hc)
summary(step3_hc)$r.squared - summary(step2_hc)$r.squared #change in Rsq= 0.1281263

#Step 4. Interaction GAP x HC_ICV on memory
step4_hc <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + HC_ICV_norm*GAP43_log_norm, data=all_data_hc)
summary(step4_hc) #significant interaction
anova(step4_hc, step3_hc)
summary(step4_hc)$r.squared - summary(step3_hc)$r.squared #change in Rsq= 0.004223594

#Interaction Plot and Post-Hoc
intx_hc <- lm(all_data_hc$ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + HC_ICV_norm*GAP43_log_norm, data=all_data_hc)
intx_hc_plot <- interact_plot(model=intx_hc, pred=GAP43_log_norm, modx=HC_ICV_norm, plot.points=TRUE, data=all_data_hc)
intx_hc_plot

sim_slopes(intx_hc, pred=GAP43_log_norm, modx=HC_ICV_norm) #sig slope at -1SD and mean
#Greater GAP-43 is associated with worse episodic memory specifically in those with low HC volume and mean HC volume
#Contrastingly, those with higher HC volume do not exhibit this negative relationship between GAP-43 and memory

#Emmeans for post-hoc
effa <- mean(all_data_hc$HC_ICV_norm) + sd(all_data_hc$HC_ICV_norm)
eff <- mean(all_data_hc$HC_ICV_norm)
effb <- mean(all_data_hc$HC_ICV_norm) - sd(all_data_hc$HC_ICV_norm)

effar <- round(effa,1) #1
effr <- round(eff,1) #0 
effbr <- round(effb,1) #-1

mylist <- list(HC_ICV_norm=c(effbr,effr,effar))

#Test significance of the slope
emtrends(intx_hc , ~HC_ICV_norm, var="GAP43_log_norm", at=mylist) #shows CI for each group (significant for -1SD and mean)

#Test pairwise differences of slopes
emtrends(intx_hc, pairwise ~HC_ICV_norm, var="GAP43_log_norm", at=mylist, adjust="none")
#p value matches that shown in step3 of HLR
#Yes, all slopes are significantly different from each other


#------------------------------------------------------------------------------#
#2. CA1 subfield
#Step 1. Covariates predicting memory
step1_ca1 <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data=all_data_hc)
summary(step1_ca1)

#Step 2. Main effect of GAP-43 predicting memory
step2_ca1 <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm, data=all_data_hc)
summary(step2_ca1)
anova(step1_ca1, step2_ca1)

#Step 3. Main effect of HC
step3_ca1 <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm + CA1_sum_ICV, data=all_data_hc)
summary(step3_ca1) #significant
anova(step2_ca1, step3_ca1) #significant
summary(step3_ca1)$r.squared - summary(step2_ca1)$r.squared #change in Rsq= 0.02324049

#Step 4. Interaction GAP x HC_ICV ROI on memory
step4_ca1 <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm*CA1_sum_ICV, data=all_data_hc)
summary(step4_ca1)  #ns interaction
anova(step3_ca1,step4_ca1)
summary(step4_ca1)$r.squared - summary(step3_ca1)$r.squared #change in Rsq= 0.0214


#No post-hoc since interaction is not significant



#------------------------------------------------------------------------------#
#3. DG subfield
#Step 1. Covariates predicting memory
step1_dg <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data=all_data_hc)
summary(step1_dg)

#Step 2. Main effect of GAP-43 predicting memory
step2_dg <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm, data=all_data_hc)
summary(step2_dg)
anova(step1_dg, step2_dg)

#Step 3. Main effect of HC
step3_dg <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm + DG_sum_ICV, data=all_data_hc)
summary(step3_dg) #significant
anova(step2_dg, step3_dg) #significant
summary(step3_dg)$r.squared - summary(step2_dg)$r.squared #change in Rsq= 0.08201847

#Step 4. Interaction GAP x HC_ICV ROI on memory
step4_dg <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm*DG_sum_ICV, data=all_data_hc)
summary(step4_dg) 
anova(step3_dg, step4_dg) #significant interaction
summary(step4_dg)$r.squared - summary(step3_dg)$r.squared #change in Rsq= 0.00461105

#Interaction Plot and Post-Hoc
intx_dg <- lm(all_data_hc$ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + DG_sum_ICV*GAP43_log_norm, data=all_data_hc)
intx_dg_plot <- interact_plot(model=intx_dg, pred=GAP43_log_norm, modx=DG_sum_ICV, plot.points=TRUE, data=all_data_hc)
intx_dg_plot

sim_slopes(intx_dg, pred=GAP43_log_norm, modx=DG_sum_ICV) #significant slope for all groups

#Emmeans for post-hoc
effa <- mean(all_data_hc$DG_sum_ICV) + sd(all_data_hc$DG_sum_ICV)
eff <- mean(all_data_hc$DG_sum_ICV)
effb <- mean(all_data_hc$DG_sum_ICV) - sd(all_data_hc$DG_sum_ICV)

effar <- round(effa,1) #1
effr <- round(eff,1) #0 
effbr <- round(effb,1) #-1

mylist <- list(DG_sum_ICV=c(effbr,effr,effar))

#Test significancc of the slope
emtrends(intx_dg, ~DG_sum_ICV, var="GAP43_log_norm", at=mylist) #shows CI for each group (significant for all groups)

#Test pairwise differences of slopes
emtrends(intx_dg, pairwise ~DG_sum_ICV, var="GAP43_log_norm", at=mylist, adjust="none")
#Yes, all slopes are significantly different from each other



#------------------------------------------------------------------------------#
#4. Subiculum subfield
#Step 1. Covariates predicting memory
step1_sub <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data=all_data_hc)
summary(step1_sub)

#Step 2. Main effect of GAP-43 predicting memory
step2_sub <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm, data=all_data_hc)
summary(step2_sub)
anova(step1,step2_sub)

#Step 3. Main effect of HC
step3_sub <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm + Subiculum_sum_ICV, data=all_data_hc)
summary(step3_sub) #significant
anova(step2_sub, step3_sub) #significant
summary(step3_sub)$r.squared - summary(step2_sub)$r.squared #change in Rsq= 0.1139697

#Step 4. Interaction GAP x HC_ICV ROI on memory
step4_sub <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm*Subiculum_sum_ICV, data=all_data_hc)
summary(step4_sub) 
anova(step3_sub, step4_sub)
summary(step4_sub)$r.squared - summary(step3_sub)$r.squared #change in Rsq= 0.004325508

#Interaction plot and post-hoc
intx_sub <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm*Subiculum_sum_ICV, data=all_data_hc)
posthoc_sub <- interact_plot(model=intx_sub, pred=GAP43_log_norm, modx=Subiculum_sum_ICV, plot.points=TRUE, data=all_data_hc)
posthoc_sub
sim_slopes(intx_sub, pred=GAP43_log_norm, modx=Subiculum_sum_ICV)

#Emmeans for post-hoc
effa <- mean(all_data_hc$Subiculum_sum_ICV) + sd(all_data_hc$Subiculum_sum_ICV)
eff <- mean(all_data_hc$Subiculum_sum_ICV)
effb <- mean(all_data_hc$Subiculum_sum_ICV) - sd(all_data_hc$Subiculum_sum_ICV)

effar <- round(effa,1) #1
effr <- round(eff,1) #0 
effbr <- round(effb,1) #-1

mylist <- list(Subiculum_sum_ICV=c(effbr,effr,effar))

#Test significance of the slope
emtrends(intx_sub, ~Subiculum_sum_ICV, var="GAP43_log_norm", at=mylist) #shows CI for each group (significant for all groups)

#Test pairwise differences of slopes
emtrends(intx_sub, pairwise ~Subiculum_sum_ICV, var="GAP43_log_norm", at=mylist, adjust="none")
#Yes, all slopes are significantly different from each other



#------------------------------------------------------------------------------#
#5. Presubiculum subfield
#Step 1. Covariates predicting memory
step1_presub <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data=all_data_hc)
summary(step1_presub)

#Step 2. Main effect of GAP-43 predicting memory
step2_presub <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm, data=all_data_hc)
summary(step2_presub)
anova(step1_presub, step2_presub)

#Step 3. Main effect of HC
step3_presub <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm + Presubiculum_sum_ICV, data=all_data_hc)
summary(step3_presub) #significant
anova(step2_presub, step3_presub) #significant
summary(step3_presub)$r.squared - summary(step2_presub)$r.squared #change in Rsq= 0.1395107

#Step 4. Interaction GAP x HC_ICV ROI on memory
step4_presub <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm*Presubiculum_sum_ICV, data=all_data_hc)
summary(step4_presub) 
anova(step3_presub, step4_presub)
summary(step4_presub)$r.squared - summary(step3_presub)$r.squared #change in Rsq= 0.005541383

#post hoc for presubiculum interaction
presub_intx <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm*Presubiculum_sum_ICV, data=all_data_hc)
posthoc_presub <- interact_plot(model=presub_intx, pred=GAP43_log_norm, modx=Presubiculum_sum_ICV, plot.points=TRUE, data=all_data_hc)
posthoc_presub
sim_slopes(presub_intx, pred=GAP43_log_norm, modx=Presubiculum_sum_ICV) 

#Emmeans for post-hoc
effa <- mean(all_data_hc$Presubiculum_sum_ICV) + sd(all_data_hc$Presubiculum_sum_ICV)
eff <- mean(all_data_hc$Presubiculum_sum_ICV)
effb <- mean(all_data_hc$Presubiculum_sum_ICV) - sd(all_data_hc$Presubiculum_sum_ICV)

effar <- round(effa,1) #1
effr <- round(eff,1) #0 
effbr <- round(effb,1) #-1

mylist <- list(Presubiculum_sum_ICV=c(effbr,effr,effar))

#Test significance of the slope
emtrends(presub_intx, ~Presubiculum_sum_ICV, var="GAP43_log_norm", at=mylist) #shows CI for each group (significant for all groups)

#Test pairwise differences of slopes
emtrends(presub_intx, pairwise ~Presubiculum_sum_ICV, var="GAP43_log_norm", at=mylist, adjust="none")
#Yes, all slopes are significantly different from each other

#------------------------------------------------------------------------------#
#------------------ Figures: Significant Effects ------------------------------#
#------------------------------------------------------------------------------#

library(ggpubr)
library(jtools) #theme_apa()

# FIGURE 1: MAIN EFFECTS
## Main effect of GAP43 (panel A)

# create residuals
resid_memory <- residuals(lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data = all_data_hc))
resid_gap43 <- residuals(lm(GAP43_log_norm ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data = all_data_hc))

# create panel A
a <- ggplot(data = all_data_hc, aes(x = resid_gap43, y = resid_memory)) +
  geom_point(color = "skyblue", alpha=.3) +
  geom_smooth(method = "lm", se = TRUE, color = "#1A73B5") +
  labs(x = "CSF GAP-43 (pg/mL)", y = "Composite Episodic Memory\n(z-score)", title = "") +
  theme_apa() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 12),
        # plot.title = element_text(size=11),
        panel.background = element_rect(fill='transparent'),
        plot.title.position = "plot")
a

## Main effect of hippocampal ROIs (panel B)

# create residuals
resid_memory <- residuals(lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data = all_data_hc))
resid_hc <- residuals(lm(HC_ICV_norm ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data = all_data_hc))
resid_sub <- residuals(lm(Subiculum_sum_ICV ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data = all_data_hc))
resid_presub <- residuals(lm(Presubiculum_sum_ICV ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data = all_data_hc))
resid_dg <- residuals(lm(DG_sum_ICV ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data = all_data_hc))
resid_ca1 <- residuals(lm(CA1_sum_ICV ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group, data = all_data_hc))
partial_df <- data.frame(resid_memory, resid_hc, resid_sub, resid_presub, resid_dg, resid_ca1)

# create plots for panel B
b1 <- ggplot(data = all_data_hc, aes(x = resid_hc, y = resid_memory)) +
  geom_point(color = "orange", alpha=.3) +
  geom_smooth(method = "lm", se = TRUE, color = "#E59F14") +
  labs(x = "Total Hippocampal Volume", y = "Composite Episodic Memory\n(z-score)", title = "") +
  theme_apa() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill='transparent'),
        plot.title.position = "panel"
  )

b2 <- ggplot(data = all_data_hc, aes(x = resid_sub, y = resid_memory)) +
  geom_point(color = "orange", alpha=.3) +
  geom_smooth(method = "lm", se = TRUE, color = "#E59F14") +
  labs(x = "Subiculum Volume", y = "Composite Episodic Memory\n(z-score)", title = "") +
  theme_apa() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 12)
  )

b3 <- ggplot(data = all_data_hc, aes(x = resid_presub, y = resid_memory)) +
  geom_point(color = "orange", alpha=.3) +
  geom_smooth(method = "lm", se = TRUE, color = "#E59F14") +
  labs(x = "Presubiculum Volume", y = "Composite Episodic Memory\n(z-score)", title = "") +
  theme_apa() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 12)
  )

b4 <- ggplot(data = all_data_hc, aes(x = resid_dg, y = resid_memory)) +
  geom_point(color = "orange", alpha=.3) +
  geom_smooth(method = "lm", se = TRUE, color = "#E59F14") +
  labs(x = "Dentate Gyrus Volume", y = "Composite Episodic Memory\n(z-score)", title = "") +
  theme_apa() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 12)
  )

b5 <- ggplot(data = all_data_hc, aes(x = resid_ca1, y = resid_memory)) +
  geom_point(color = "orange", alpha=.3) +
  geom_smooth(method = "lm", se = TRUE, color = "#E59F14") +
  labs(x = "CA1 Volume", y = "Composite Episodic Memory\n(z-score)", title = "") +
  theme_apa() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 12)
  )

## Combine panels to create Figure 1
library(gridExtra)

fig1 <- ggarrange(a, b1, b2, b3, b4, b5,
                  nrow = 3, ncol = 2,
                  labels = c("A. Step 2", "B. Step 3"),
                  label.x = -0.11,
                  font.label = list(size = 13, face = "bold"))
fig1
ggsave("data/Figure1.jpg",
       fig1, bg = "transparent", dpi=300,
       width = 5.5, 
       height = 9.3, 
       units = "in")

# FIGURE 2: INTERACTIONS
# CATEGORICAL MODERATOR FOR VISUALIZATION PURPOSES

# create new categorical variables so that above 0.5 = +1SD group, -0.5 - 0.5 SD = mean group, and below -0.5 = -1SD group
all_data_hc$cat_totalhc <- ifelse(all_data_hc$HC_ICV_norm > 0.5, "+1 SD", ifelse(
                                    all_data_hc$HC_ICV_norm < -0.5, "-1 SD", "Mean"))
all_data_hc$cat_totalhc <- factor(all_data_hc$cat_totalhc, levels = c("+1 SD", "Mean", "-1 SD"))

all_data_hc$cat_sub <- ifelse(all_data_hc$Subiculum_sum_ICV > 0.5, "+1 SD", ifelse(
                                all_data_hc$Subiculum_sum_ICV < -0.5, "-1 SD", "Mean"))
all_data_hc$cat_sub <- factor(all_data_hc$cat_sub, levels = c("+1 SD", "Mean", "-1 SD"))

all_data_hc$cat_presub <- ifelse(all_data_hc$Presubiculum_sum_ICV > 0.5, "+1 SD", ifelse(
                                all_data_hc$Presubiculum_sum_ICV < -0.5, "-1 SD", "Mean"))
all_data_hc$cat_presub <- factor(all_data_hc$cat_presub, levels = c("+1 SD", "Mean", "-1 SD"))

all_data_hc$cat_dg <- ifelse(all_data_hc$DG_sum_ICV > 0.5, "+1 SD", ifelse(
                              all_data_hc$DG_sum_ICV < -0.5, "-1 SD", "Mean"))
all_data_hc$cat_dg <- factor(all_data_hc$cat_dg, levels = c("+1 SD", "Mean", "-1 SD"))

# run linear models (need to do with variables as numeric)
all_data_hc <- all_data_hc %>%
  mutate(
    AGE_norm = as.numeric(AGE_norm),
    PTEDUCAT_norm = as.numeric(PTEDUCAT_norm),
    GAP43_log_norm = as.numeric(GAP43_log_norm),
    HC_ICV_norm = as.numeric(HC_ICV_norm),
    Subiculum_sum_ICV = as.numeric(Subiculum_sum_ICV),
    Presubiculum_sum_ICV = as.numeric(Presubiculum_sum_ICV),
    DG_sum_ICV = as.numeric(DG_sum_ICV)
  )

intx_hc <- lm(all_data_hc$ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + HC_ICV_norm*GAP43_log_norm, data=all_data_hc)
intx_sub <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm*Subiculum_sum_ICV, data=all_data_hc)
presub_intx <- lm(ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + GAP43_log_norm*Presubiculum_sum_ICV, data=all_data_hc)
intx_dg <- lm(all_data_hc$ADNI_MEM ~ AGE_norm + PTEDUCAT_norm + PTGENDER + APOE4 + ATN_group + DG_sum_ICV*GAP43_log_norm, data=all_data_hc)

# get coefficients from model
coef.hc <- coef(intx_hc)
coef.sub <- coef(intx_sub)
coef.presub <- coef(presub_intx)
coef.dg <- coef(intx_dg)

# define levels of hc volume (+/- 1SD, mean)
mean_hc <- mean(all_data_hc$HC_ICV_norm, na.rm = TRUE)
sd_hc <- sd(all_data_hc$HC_ICV_norm, na.rm = TRUE)
hc_levels <- c(mean_hc - sd_hc, mean_hc, mean_hc + sd_hc)
regression_data.hc <- data.frame(
  cat_totalhc = factor(c("-1 SD", "Mean", "+1 SD"), levels = c("+1 SD", "Mean", "-1 SD")),
  intercept = coef.hc["(Intercept)"] + coef.hc["HC_ICV_norm"] * hc_levels,
  slope = coef.hc["GAP43_log_norm"] + coef.hc["HC_ICV_norm:GAP43_log_norm"] * hc_levels
)

mean_sub <- mean(all_data_hc$Subiculum_sum_ICV, na.rm = TRUE)
sd_sub <- sd(all_data_hc$Subiculum_sum_ICV, na.rm = TRUE)
sub_levels <- c(mean_sub - sd_sub, mean_sub, mean_sub + sd_sub)
regression_data.sub <- data.frame(
  cat_sub = factor(c("-1 SD", "Mean", "+1 SD"), levels = c("+1 SD", "Mean", "-1 SD")),
  intercept = coef.sub["(Intercept)"] + coef.sub["Subiculum_sum_ICV"] * sub_levels,
  slope = coef.sub["GAP43_log_norm"] + coef.sub["GAP43_log_norm:Subiculum_sum_ICV"] * sub_levels
)

mean_presub <- mean(all_data_hc$Presubiculum_sum_ICV, na.rm = TRUE)
sd_presub <- sd(all_data_hc$Presubiculum_sum_ICV, na.rm = TRUE)
presub_levels <- c(mean_presub - sd_presub, mean_presub, mean_presub + sd_presub)
regression_data.presub <- data.frame(
  cat_presub = factor(c("-1 SD", "Mean", "+1 SD"), levels = c("+1 SD", "Mean", "-1 SD")),
  intercept = coef.presub["(Intercept)"] + coef.presub["Presubiculum_sum_ICV"] * presub_levels,
  slope = coef.presub["GAP43_log_norm"] + coef.presub["GAP43_log_norm:Presubiculum_sum_ICV"] * presub_levels
)

mean_dg <- mean(all_data_hc$DG_sum_ICV, na.rm = TRUE)
sd_dg <- sd(all_data_hc$DG_sum_ICV, na.rm = TRUE)
dg_levels <- c(mean_dg - sd_dg, mean_dg, mean_dg + sd_dg)
regression_data.dg <- data.frame(
  cat_dg = factor(c("-1 SD", "Mean", "+1 SD"), levels = c("+1 SD", "Mean", "-1 SD")),
  intercept = coef.dg["(Intercept)"] + coef.dg["DG_sum_ICV"] * dg_levels,
  slope = coef.dg["GAP43_log_norm"] + coef.dg["DG_sum_ICV:GAP43_log_norm"] * dg_levels
)

# create plot
p1 <- ggplot(data = all_data_hc, aes(x = GAP43_log_norm, y = ADNI_MEM, color = cat_totalhc, shape = cat_totalhc, linetype = cat_totalhc)) +
  geom_point(alpha = .3, size = 2.5) +
  geom_abline(data = regression_data.hc, aes(intercept = intercept, slope = slope, color = cat_totalhc, linetype = cat_totalhc), 
              size = 1.25, show.legend = FALSE) +
  geom_segment(data = regression_data.hc, 
               aes(x = -Inf, xend = -Inf, y = -Inf, yend = -Inf, color = cat_totalhc, linetype = cat_totalhc), 
               size = 1.5, inherit.aes = FALSE, show.legend = TRUE) +
  scale_color_manual(values = c("black", "#CC79A7", "#019E73")) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed")) +
  labs(title = "A. Total Hippocampal Volume",
       x = "CSF GAP-43 (pg/mL)", y = "Composite Episodic Memory\n(z-score)",
       color = "Hippocampal ROI Volume", shape = "Hippocampal ROI Volume", linetype = "Hippocampal ROI Volume") +
  theme_apa(legend.use.title = TRUE) +
  theme(legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title=element_text(size=14)) + 
  guides(shape = guide_legend(override.aes = list(size = 5)))

p2 <- ggplot(data = all_data_hc, aes(x = GAP43_log_norm, y = ADNI_MEM, color = cat_sub, shape = cat_sub, linetype = cat_sub)) +
  geom_point(alpha = .3, size = 2.5) +
  geom_abline(data = regression_data.sub, aes(intercept = intercept, slope = slope, color = cat_sub, linetype = cat_sub), 
              size = 1.25, show.legend = FALSE) +
  geom_segment(data = regression_data.sub, 
               aes(x = -Inf, xend = -Inf, y = -Inf, yend = -Inf, color = cat_sub, linetype = cat_sub), 
               size = 1.5, inherit.aes = FALSE, show.legend = TRUE) +
  scale_color_manual(values = c("black", "#CC79A7", "#019E73")) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed")) +
  labs(title = "B. Subiculum Volume",
       x = "CSF GAP-43 (pg/mL)", y = "Composite Episodic Memory\n(z-score)",
       color = "Hippocampal ROI Volume", shape = "Hippocampal ROI Volume", linetype = "Hippocampal ROI Volume") +
  theme_apa(legend.use.title = TRUE) +
  theme(legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title=element_text(size=14)) + 
  guides(shape = guide_legend(override.aes = list(size = 5)))

p3 <- ggplot(data = all_data_hc, aes(x = GAP43_log_norm, y = ADNI_MEM, color = cat_presub, shape = cat_presub, linetype = cat_presub)) +
  geom_point(alpha = .3, size = 2.5) +
  geom_abline(data = regression_data.presub, aes(intercept = intercept, slope = slope, color = cat_presub, linetype = cat_presub), 
              size = 1.25, show.legend = FALSE) +
  geom_segment(data = regression_data.presub, 
               aes(x = -Inf, xend = -Inf, y = -Inf, yend = -Inf, color = cat_presub, linetype = cat_presub), 
               size = 1.5, inherit.aes = FALSE, show.legend = TRUE) +
  scale_color_manual(values = c("black", "#CC79A7", "#019E73")) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed")) +
  labs(title = "C. Presubiculum Volume",
       x = "CSF GAP-43 (pg/mL)", y = "Composite Episodic Memory\n(z-score)",
       color = "Hippocampal ROI Volume", shape = "Hippocampal ROI Volume", linetype = "Hippocampal ROI Volume") +
  theme_apa(legend.use.title = TRUE) +
  theme(legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title=element_text(size=14)) + 
  guides(shape = guide_legend(override.aes = list(size = 5)))

p4 <- ggplot(data = all_data_hc, aes(x = GAP43_log_norm, y = ADNI_MEM, color = cat_dg, shape = cat_dg, linetype = cat_dg)) +
  geom_point(alpha = .3, size = 2.5) +
  geom_abline(data = regression_data.dg, aes(intercept = intercept, slope = slope, color = cat_dg, linetype = cat_dg), 
              size = 1.25, show.legend = FALSE) +
  geom_segment(data = regression_data.dg, 
               aes(x = -Inf, xend = -Inf, y = -Inf, yend = -Inf, color = cat_dg, linetype = cat_dg), 
               size = 1.5, inherit.aes = FALSE, show.legend = TRUE) +
  scale_color_manual(values = c("black", "#CC79A7", "#019E73")) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed")) +
  labs(title = "D. Dentate Gyrus Volume",
       x = "CSF GAP-43 (pg/mL)", y = "Composite Episodic Memory\n(z-score)",
       color = "Hippocampal ROI Volume", shape = "Hippocampal ROI Volume", linetype = "Hippocampal ROI Volume") +
  theme_apa(legend.use.title = TRUE) +
  theme(legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title=element_text(size=14)) + 
  guides(shape = guide_legend(override.aes = list(size = 5)))

fig2 <- ggpubr::ggarrange(
  p1, p2, p3, p4,
  common.legend = TRUE,
  legend = "top",
  align = "hv",
  nrow = 2, ncol = 2
)
fig2

ggsave("data/Figure2.jpg", fig2, width = 190, units = "mm", dpi = 300, bg = "white")


#------------------------------ End of Code------------------------------------#

#For notes from previous abstract runs (INS) including triple interaction, see code 10-13-2024.R

