##########################################################################################################
#### SEPT 2016 10yr HF risk prediction from 80 proteins + ARIC risk factors   ############################
#### Christoph Nowak - christoph.nowak@medsci.uu.se   ####################################################
#### cohorts: ULSAM-77 & PIVUS-70   ######################################################################

dir <- ("~/Proteomics_HF")
setwd(dir)

library(readstata13)
library(foreign)
library(psych)
library(plyr)
library(splines)
library(survival)
library(glmnet)
library(rms)
library(boot)

######## 1) Pre-Processing: Upload, data-cleaning, variable re-scaling, merging #######

HF <- data.frame(read.dta13("ULSAM77 PIVUS70 all heart failure events.dta"))
ogtt <- data.frame(read.dta("age70-clamp-ogtt.dta"))
kont <- data.frame(read.dta("ULSAM_kont.dta"))
u <- merge(ogtt,kont,by="pat")

u$dat70_numeric <- as.numeric(as.Date(as.character(u$dat70),"%Y%m%d"))  
d <- data.frame(read.dta("d140916.dta"))
u <- merge(u,d,by="pat")

hr <- data.frame(read.dta("u77_v642_5.dta")) ## Ulsam heart rate
ami <- data.frame(read.dta13("ULSAM77 prevalent AMI.dta")) ## Ulsam acute MI

u <- merge(u,hr,by="pat")
u <- merge(u,ami,by="pat")

ingel70 <- data.frame(read.table("ingel70.txt"))
ingel70$id <- ingel70$pat
Ins <- data.frame(read.dta("PIVUS-Insulin-CD40-LPK.dta"))
biokem <- data.frame(read.dta("PIVUS data klin kemi.dta"))
biokem$id <- biokem$lpnr
basal <- data.frame(read.dta("PIVUS-basalfil_70.dta"))

p <- merge(basal,ingel70,by="id") 
p <- merge(p,Ins,by="id") 
p <- merge(p,biokem,by="id")

######## 2) 10-yr ARIC-HF risk score prep: Sex, ""Race"", CurrentSmoker, FormerSmoker, Age, HR, SBP, BMI, Hx of CHD, aHT-Tx, currentDiabetes
  
p$AGE <- p$age70.x ## Age
u$AGE <- u$age77
u$SEX <- 0 ## Female "1"
p$SEX <- p$kn
u$CurrentSmo <- 0; u$CurrentSmo[u$a117 == 1] <- 1 ## Tobacco smoking
u$formerSmo <- 0; u$formerSmo[u$a118 == 1 | u$a119 == 1] <- 1 
p$CurrentSmo <- 0; p$CurrentSmo[p$rkarenu==1] <- 1
p$formerSmo <- 0; p$formerSmo[p$tidigarerkare==1] <- 1
u$SBP <- u$v013 ## SBP
p$SBP <- p$manuelltsbp # manually-measured SBP because ARIC is community physician-focused
p$HR <- p$hr ## HR
u$HR <- u$V645 
u$BMI <- u$v290 ## BMI
p$BMI <- p$bmi
p$DM <- (p$diabetesmell)  # Diabetes mellitus: clinically registered 
                          # not based on blood glucose because ARIC asks persons "Do you have DM?"
u$DM <- u$v378
u$aHT <- u$v101 ## antihypertensive Medx
p$aHT <- p$medicinmothgtbt
u$AMI <- u$amiprev ## prevalent acute MI
p$AMI <- p$validmi70

  ## so-called "race" which is ARIC is irrelevant as homogeneous European cohorts

u$id2 <- u$pat+10000
p$id2 <- p$id
l <- list(u,p)
u_p <- do.call(rbind.fill, l)
HF$id2 <- HF$id
HF$id2[HF$study == "ULSAM"] <- HF$id2[HF$study == "ULSAM"]+10000
a <- merge(u_p,HF,by="id2")
a$Cohort <- 0; a$Cohort[a$study == "ULSAM"] <- 1 # cohort dummy (ULSAM coded)

  ## table(a$SEX) 
  ## table(a$aHT) 
  ## table(a$DM)
  ## table(a$CurrentSmo) 
  ## table(a$formerSmo) 
  ## table(a$AMI) 
  ## describe(a$AGE) 
  ## describe(a$BMI) 
  ## describe(a$SBP) 

######## 3) establish the random 2/3 training, 1/3 validation split 

  ## custom-made split-function ##

split_data <- function(dat, props = c(.6666666, .33333), which.adjust = 1){
  stopifnot(all(props >= 0), which.adjust <= length(props))
  props <- props/sum(props)
  n <- nrow(dat)
  ns <- round(n * props)
  ns[which.adjust] <- n - sum(ns[-which.adjust])
  ids <- rep(1:length(props), ns)
  which.group <- sample(ids)
  split(dat, which.group)
}

set.seed(93468)
split <- split_data(a,c(0.666666,0.3333)) 
HF_train <- split[[1]]
HF_val <- split[[2]]
HF_combined <- a

  ## Mean-imputation of NAs => heart rate has 78 and 32 missing in Train and Val sets; NT-pro-bnp has 1 and 1 missing in Train and Val sets
HF_train[is.na(HF_train[,"HR"]),"HR"] <- mean(HF_train[,"HR"],  na.rm = TRUE) 
HF_train[is.na(HF_train[,"stdntprobnp"]),"stdntprobnp"] <- mean(HF_train[,"stdntprobnp"],  na.rm = TRUE)
HF_val[is.na(HF_val[,"HR"]),"HR"] <- mean(HF_val[,"HR"],  na.rm = TRUE) 
HF_val[is.na(HF_val[,"stdntprobnp"]),"stdntprobnp"] <- mean(HF_val[,"stdntprobnp"],  na.rm = TRUE)


######## 4) Penalized LASSO-Cox in n = 915 training set ### ### ###

s = HF_train
ARIC <- c("AGE","BMI","SBP","Cohort","SEX","aHT","HR","DM","CurrentSmo","formerSmo","AMI") 
all_std_prot <- (names(s)[(grep("^std",names(s)))])

  ## exclude from these proteins those that didn't pass quality control (see manuscript)
Proteins <- all_std_prot[-which(all_std_prot %in% c("stdntprobnp","stdprot_beta_ngf","stdprot_mmp_7","stdprot_sirt2","stdprot_nemo" ))]
y = Surv(as.numeric(as.Date(s$allsviktdatestata)),s$allsvikt)
x = as.matrix(s[,c(ARIC,Proteins)])  # matrix for wo.-NTproBNP
x_probnp = as.matrix(s[,c(ARIC,"stdntprobnp",Proteins)])  # with 

  ## penalties for ARIC +- ntproBNP  
pen <- rep(1,ncol(x)) 
pen_probnp1 <- rep(1,ncol(x_probnp)) 
pen_probnp2 <- rep(1,ncol(x_probnp)) 

pen[which(colnames(x) %in% c(ARIC,"Cohort"))] <- 0 # Model 1: no ntprobnp
pen_probnp1[which(colnames(x_probnp) %in% c(ARIC,"Cohort"))] <- 0 # Model 2: ntprobnp included but not forced inside 
pen_probnp2[which(colnames(x_probnp) %in% c(ARIC,"Cohort","stdntprobnp"))] <- 0 # Model 2: ntprobnp forced inside

set.seed(1234)
lasso_wo_probnp <- cv.glmnet(x, y, family="cox",alpha=1, type.measure="deviance", nfolds=10,  penalty.factor = pen, nlambda=100) 
set.seed(1234)
lasso_with_probnp_free <- cv.glmnet(x_probnp, y, family="cox",alpha=1, type.measure="deviance", nfolds=10,  penalty.factor = pen_probnp1, nlambda=100) 
set.seed(1234)
lasso_with_probnp_forced <- cv.glmnet(x_probnp, y, family="cox",alpha=1, type.measure="deviance", nfolds=10,  penalty.factor = pen_probnp2, nlambda=100) 

  ## find lambda with mean squared error being minimum and any number of proteins
ind1  <- which(lasso_wo_probnp$glmnet.fit$df < max(lasso_wo_probnp$glmnet.fit$df))
ind2  <- which(lasso_with_probnp_free$glmnet.fit$df < max(lasso_with_probnp_free$glmnet.fit$df))
ind3  <- which(lasso_with_probnp_forced$glmnet.fit$df < max(lasso_with_probnp_forced$glmnet.fit$df))

  ## or: select max. 10 proteins
  # ind1  <- which(  lasso_wo_probnp$glmnet.fit$df < 20) 
  # ind2  <- which(  lasso_with_probnp_free$glmnet.fit$df < 20)
  # ind3  <- which(  lasso_with_probnp_forced$glmnet.fit$df < 21)

optimal.lambda.index1 <- which(lasso_wo_probnp$cvm == min(lasso_wo_probnp$cvm[ind1]))
optimal.lambda1  <- lasso_wo_probnp$lambda[optimal.lambda.index1]
optimal.lambda.index2 <- which(lasso_with_probnp_free$cvm == min(lasso_with_probnp_free$cvm[ind2]))
optimal.lambda2  <- lasso_with_probnp_free$lambda[optimal.lambda.index2]
optimal.lambda.index3 <- which(lasso_with_probnp_forced$cvm == min(lasso_with_probnp_forced$cvm[ind3]))
optimal.lambd3  <- lasso_with_probnp_forced$lambda[optimal.lambda.index3]

  ## selected variables
optimal.beta1  <- lasso_wo_probnp$glmnet.fit$beta[,optimal.lambda.index1] 
selectedBeta1 <- optimal.beta1[abs(optimal.beta1)>0 ] 
optimal.beta2  <- lasso_with_probnp_free$glmnet.fit$beta[,optimal.lambda.index2] 
selectedBeta2 <- optimal.beta2[abs(optimal.beta2)>0 ] 
optimal.beta3  <- lasso_with_probnp_forced$glmnet.fit$beta[,optimal.lambda.index3] 
selectedBeta3 <- optimal.beta3[abs(optimal.beta3)>0 ] 

## View((sort(names(selectedBeta1))))
## View((sort(names(selectedBeta2))))
## View((sort(names(selectedBeta3))))

  ## C-index in the validation sample
s2 = HF_val
y2 = Surv(as.numeric(as.Date(s2$allsviktdatestata)),s2$allsvikt)
x2 = as.matrix(s2[,c(ARIC,Proteins)])  ## cohort, ARIC, 80 proteins
x2pro = as.matrix(s2[,c(ARIC,"stdntprobnp",Proteins)])  ## cohort, ARIC, ntprobnp, 80 proteins

Cv1 <- survConcordance (y2 ~ predict(lasso_wo_probnp,newx=x2,s="lambda.min"))$concordance
Cvse1 <- survConcordance (y2 ~ predict(lasso_wo_probnp,newx=x2,s="lambda.min"))$std.err
Cv2 <- survConcordance (y2 ~ predict(lasso_with_probnp_free,newx=x2pro,s="lambda.min"))$concordance
Cvse2 <- survConcordance (y2 ~ predict(lasso_with_probnp_free,newx=x2pro,s="lambda.min"))$std.err
Cv3 <- survConcordance (y2 ~ predict(lasso_with_probnp_forced,newx=x2pro,s="lambda.min"))$concordance
Cvse3 <- survConcordance (y2 ~ predict(lasso_with_probnp_forced,newx=x2pro,s="lambda.min"))$std.err

## paste(Cv1,Cvse1,Cv1-(1.96*Cv1),Cv1+(1.96*Cvse1))
## paste(Cv2,Cvse2,Cv2-(1.96*Cv2),Cv2+(1.96*Cvse2))
## paste(Cv3,Cvse3,Cv3-(1.96*Cv3),Cv3+(1.96*Cvse3))

  #### ARIC factors-only C-index: without and with ntprobnp

cox1 <- (coxph(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ AGE+BMI+SBP+Cohort+SEX+aHT+HR+DM+CurrentSmo+formerSmo+AMI,data=HF_val))
cox2 <- (coxph(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ AGE+BMI+SBP+Cohort+SEX+aHT+HR+DM+CurrentSmo+formerSmo+AMI+
                 stdprot_agrp+stdprot_casp_8+stdprot_ccl20+stdprot_egf+stdprot_esm_1+stdprot_fgf_23+stdprot_fs+stdprot_gh+stdprot_hk11+
                 stdprot_il_18+stdprot_il_1ra+stdprot_il_6+stdprot_il27_a+stdprot_mmp_10+stdprot_mpo+stdprot_rage+stdprot_scf+stdprot_sele+
                 stdprot_spon1+stdprot_st2+stdprot_tie2+stdprot_tim+stdprot_trail+stdprot_u_par,data=HF_val))
anova(cox1,cox2)

C1 <- survConcordance(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ predict(cox1),data=HF_val)$concordance
Cse1 <- survConcordance(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ predict(cox1),data=HF_val)$std.err
C2 <- survConcordance(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ predict(cox2),data=HF_val)$concordance
Cse2 <- survConcordance(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ predict(cox2),data=HF_val)$std.err

## paste("C-index ARIC-only, ",C1,Cse1,(C1-1.96*Cse1),(C1+1.96*Cse1))     
## paste("C-index ARIC-only, ",C2,Cse2,(C2-1.96*Cse2),(C2+1.96*Cse2)) 

cox3 <- (coxph(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ stdntprobnp+AGE+BMI+SBP+Cohort+SEX+aHT+HR+DM+CurrentSmo+formerSmo+AMI,data=HF_val))
cox4 <- (coxph(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ stdntprobnp+AGE+BMI+SBP+Cohort+SEX+aHT+HR+DM+CurrentSmo+formerSmo+AMI+stdprot_agrp+
                 stdprot_casp_8+stdprot_ccl20+stdprot_il_18+stdprot_il_1ra+stdprot_opg+stdprot_scf+stdprot_sele+stdprot_st2+stdprot_tim+stdprot_trail,data=HF_val))
anova(cox3,cox4)

C3 <- survConcordance(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ predict(cox3),data=HF_val)$concordance
Cse3 <- survConcordance(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ predict(cox3),data=HF_val)$std.err
C4 <- survConcordance(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ predict(cox4),data=HF_val)$concordance
Cse4 <- survConcordance(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ predict(cox4),data=HF_val)$std.err

## paste("C-index ARIC-only ntprobnp, ",C3,Cse3,(C3-1.96*Cse3),(C3+1.96*Cse3))     
## paste("C-index ARIC-only ntprobnp, ",C4,Cse4,(C4-1.96*Cse4),(C4+1.96*Cse4)) 


######## 5) Boostrapping for C-index 95%-CI

Cindexdiff <- function(data,indices){
  data <- data[indices,] # select obs. in bootstrap sample
  # C-statistic without the biomarker:
  C1   <- survConcordance(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ predict(cox1),data=data)$concordance
  # C-statistic with the biomarker:
  C2   <- survConcordance(Surv(as.numeric(as.Date(allsviktdatestata)),allsvikt) ~ predict(cox2),data=data)$concordance
  as.numeric(C2-C1) # returns the difference
}

set.seed(123)
bootthis <- boot(HF_val, Cindexdiff, R=10000)

#########################################################################################################################