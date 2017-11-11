library(survival)
library(data.table)
library(foreach)
library(doSNOW)
library(aod)
library(RJSONIO)

setwd("E:/RaulEpicProject")
clsb <- readLines("cl_sbcl.txt")
ind_ph <- readLines("mainvars.txt")
sub_vars <- readLines("sub_vars.txt")
varlist <- c(clsb,ind_ph)
outputvars <- c("Caseclrt_Crs_Frst","Casecl_Crs_Frst","Casert_Crs_Frst","Colon_Prox_Frst","Colon_Dist_Frst","Colon_Nos_Frst")

setwd("E:/RaulEpicProject/male/bysmoking/current")

Age_Recr <- scan("age_recr",sep="\t")
Agexit_Frst <-scan("Agexit_Frst",sep="\t")
L_School <- scan("L_School",sep="\t",what="character")
Smoke_Intensity <- scan("Smoke_Intensity ",sep="\t",what="character")
Pa_Index <- scan("Pa_Index",sep="\t",what="character")
Bmi_C <- scan("Bmi_C",sep="\t")
meat_red_pr <- scan("meat_red_processed",sep="\t")
QE_CA <- scan("QE_CA",sep="\t")
QE_Energy <- scan("QE_ENERGY",sep="\t")
Alc_Re <- scan("Alc_Re",sep="\t")
QE_Fibt <- scan("QE_FIBT",sep="\t")
Len1 <- scan("len1",sep="\t")
Centr_C <- scan("Cntr_C",sep="\t",what="character")

## Non Interacting models
### log2 model ###
#varlist <- clsb # input the varList here. 
#varlist <- ind_ph[101:475] # input the varList here.


donevar <- scan("E:/RaulEpicProject/donevar.txt",sep="\t",what="character")
varlist <- varlist[which(varlist%in%donevar==FALSE)]


cl<-makeCluster(6)
registerDoSNOW(cl)
foreach (i = 386:length(varlist)) %dopar% {  
  library(RJSONIO)
  library(survival)
  library(aod)
  outputvars <- c("Caseclrt_Crs_Frst","Casecl_Crs_Frst","Casert_Crs_Frst","Colon_Prox_Frst","Colon_Dist_Frst","Colon_Nos_Frst")
  var1 <- scan(varlist[i],sep="\t")
  var1[is.na(var1)] <- 0
  var2 <- c() # holder for median as factor
  var3 <- c() # holder for median as median value

  q0ind <- which(var1 ==0)
  quants <- c()
  
  if (length(q0ind)>=length(var1)/5) {  # in case of more than 20% zeros, subset the vector, and then do quantile on child. and all zeros belong to Q1. 
    var1_ft <- var1[-q0ind]
    quants <- quantile(var1_ft,probs=c(0.20,0.40,0.60,0.80,1.00))  
  } else {
    quants <- quantile(var1,probs=c(0.20,0.40,0.60,0.80,1.00))  
  }
  
  q1ind <- which(var1 <= quants[1])
  q2ind <- which(var1 > (quants[1])& var1 < quants[2])
  q3ind <- which(var1 >= quants[2] & var1 < quants[3])
  q4ind <- which(var1 >= quants[3] & var1 < quants[4])
  q5ind <- which(var1 >= quants[4] & var1 <= quants[5])
  
  var2[q1ind] <- 0
  var2[q2ind] <- 1
  var2[q3ind] <- 2
  var2[q4ind] <- 3
  var2[q5ind] <- 4
  var2[q0ind] <- 0
  
  var3[q1ind] <- median(var1[q1ind])
  var3[q2ind] <- median(var1[q2ind])
  var3[q3ind] <- median(var1[q3ind])
  var3[q4ind] <- median(var1[q4ind])
  var3[q5ind] <- median(var1[q5ind])
  var3[q0ind] <- median(var1[q0ind])
  
  varmed <- list()
  varmed$q1$min <- min(var1[q1ind])
  varmed$q1$max <- max(var1[q1ind])
  
  varmed$q2$min <- min(var1[q2ind])
  varmed$q2$max <- max(var1[q2ind])
  
  varmed$q3$min <- min(var1[q3ind])
  varmed$q3$max <- max(var1[q3ind])
  
  varmed$q4$min <- min(var1[q4ind])
  varmed$q4$max <- max(var1[q4ind])
  
  varmed$q5$min <- min(var1[q5ind])
  varmed$q5$max <- max(var1[q5ind])
  
  qvar <- round(log((var1+0.0000001),base=2),digit=10)
  
  for (j in 1:6) { 
    outputvar <- scan(outputvars[j],sep="\t")
    
    outputcount <- list()
    outputcount$q1count <- table(outputvar[(which(var2==0))])[2]
    outputcount$q2count <- table(outputvar[(which(var2==1))])[2]
    outputcount$q3count <- table(outputvar[(which(var2==2))])[2]
    outputcount$q4count <- table(outputvar[(which(var2==3))])[2]
    outputcount$q5count <- table(outputvar[(which(var2==4))])[2]
    
    n5_cr <- summary(coxph(Surv(Age_Recr,Agexit_Frst,outputvar)~factor(var2)+strata(Len1,Centr_C),na.action=na.exclude))
    n5_mv <- summary(coxph(Surv(Age_Recr,Agexit_Frst,outputvar)~factor(var2)+factor(L_School)+factor(Smoke_Intensity)+factor(Pa_Index)+Bmi_C+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Len1,Centr_C),na.action=na.exclude))
    zz_cr <- summary(coxph(Surv(Age_Recr,Agexit_Frst,outputvar)~var3+strata(Len1,Centr_C),na.action=na.exclude))
    zz_mv <- summary(coxph(Surv(Age_Recr,Agexit_Frst,outputvar)~var3+factor(L_School)+factor(Smoke_Intensity)+factor(Pa_Index)+Bmi_C+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Len1,Centr_C),na.action=na.exclude))
    log2_cr <- summary(coxph(Surv(Age_Recr,Agexit_Frst,outputvar)~qvar+strata(Len1,Centr_C),na.action=na.exclude))
    log2_mv <- summary(coxph(Surv(Age_Recr,Agexit_Frst,outputvar)~qvar+factor(L_School)+factor(Smoke_Intensity)+factor(Pa_Index)+Bmi_C+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Len1,Centr_C),na.action=na.exclude))
    
    n5_cr$varmed <- varmed
    n5_cr$outputcount <- outputcount
    n5_mv$varmed <- varmed
    n5_mv$outputcount <- outputcount
    
    zz_cr$varmed <- varmed
    zz_cr$outputcount <- outputcount
    zz_mv$varmed <- varmed
    zz_mv$outputcount <- outputcount
    
    log2_cr$varmed <- varmed
    log2_cr$outputcount <- outputcount
    log2_mv$varmed <- varmed
    log2_mv$outputcount <- outputcount
    
        
    writeLines(toJSON(n5_cr),paste(c("male_smoking_n5_cr__",varlist[i],"__",outputvars[j]),collapse=""))  
    writeLines(toJSON(n5_mv),paste(c("male_smoking_n5_mv__",varlist[i],"__",outputvars[j]),collapse="")) 
    writeLines(toJSON(zz_cr),paste(c("male_smoking_zz_cr__",varlist[i],"__",outputvars[j]),collapse="")) 
    writeLines(toJSON(zz_mv),paste(c("male_smoking_zz_mv__",varlist[i],"__",outputvars[j]),collapse="")) 
    writeLines(toJSON(log2_cr),paste(c("male_smoking_log2_cr__",varlist[i],"__",outputvars[j]),collapse="")) 
    writeLines(toJSON(log2_mv),paste(c("male_smoking_log2_mv__",varlist[i],"__",outputvars[j]),collapse="")) 
  }
}
stopCluster(cl)
gc()                         

#### Interaction Models with Age#####
Age_Cat <- Age_Recr
Age_Cat[which(Age_Recr<65&Age_Recr>55)] <- 1
Age_Cat[which(Age_Recr>65)] <- 2
Age_Cat[which(Age_Recr<55)] <- 0

cl<-makeCluster(6)
registerDoSNOW(cl)
foreach (i = 1:length(clsb)) %dopar% {  
  library(RJSONIO)
  library(survival)
  library(aod)
  outputvars <- c("Caseclrt_Crs_Frst","Casecl_Crs_Frst","Casert_Crs_Frst","Colon_Prox_Frst","Colon_Dist_Frst","Colon_Nos_Frst")
  qvar <- scan(clsb[i],sep="\t")
  qvar[is.na(qvar)] <- 0
  qvar <- round(log((qvar+0.0000001),base=2),digit=10)

  j <- 1
  intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Age_Cat)+factor(L_School)+factor(Pa_Index)+Bmi_C+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Centr_C),na.action=na.exclude)
  waldres <- wald.test(b = coef(intmod_mv)[18:19], Sigma = vcov(intmod_mv)[c(18:19),c(18:19)], Terms=1:2)  
  res1 <- summary(intmod_mv)
  res1$aodwald <- waldres$result
  writeLines(toJSON(res1),paste(c("male__smoking__int_mod__age__",clsb[i],"__",outputvars[j]),collapse=""))
  
  j <- 2
  intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Age_Cat)+factor(L_School)+factor(Pa_Index)+Bmi_C+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Centr_C),na.action=na.exclude)
  waldres <- wald.test(b = coef(intmod_mv)[18:19], Sigma = vcov(intmod_mv)[c(18:19),c(18:19)], Terms=1:2)  
  res1 <- summary(intmod_mv)
  res1$aodwald <- waldres$result
  writeLines(toJSON(res1),paste(c("male__smoking__int_mod__age__",clsb[i],"__",outputvars[j]),collapse=""))
  
  j <- 3
  intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Age_Cat)+factor(L_School)+factor(Pa_Index)+Bmi_C+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Centr_C),na.action=na.exclude)
  waldres <- wald.test(b = coef(intmod_mv)[18:19], Sigma = vcov(intmod_mv)[c(18:19),c(18:19)], Terms=1:2)  
  res1 <- summary(intmod_mv)
  res1$aodwald <- waldres$result
  writeLines(toJSON(res1),paste(c("male__smoking__int_mod__age__",clsb[i],"__",outputvars[j]),collapse=""))
  
  j <- 4
  intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Age_Cat)+factor(L_School)+factor(Pa_Index)+Bmi_C+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Centr_C),na.action=na.exclude)
  waldres <- wald.test(b = coef(intmod_mv)[18:19], Sigma = vcov(intmod_mv)[c(18:19),c(18:19)], Terms=1:2)  
  res1 <- summary(intmod_mv)
  res1$aodwald <- waldres$result
  writeLines(toJSON(res1),paste(c("male__smoking__int_mod__age__",clsb[i],"__",outputvars[j]),collapse=""))
  
  j <- 5
  intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Age_Cat)+factor(L_School)+factor(Pa_Index)+Bmi_C+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Centr_C),na.action=na.exclude)
  waldres <- wald.test(b = coef(intmod_mv)[18:19], Sigma = vcov(intmod_mv)[c(18:19),c(18:19)], Terms=1:2)  
  res1 <- summary(intmod_mv)
  res1$aodwald <- waldres$result
  writeLines(toJSON(res1),paste(c("male__smoking__int_mod__age__",clsb[i],"__",outputvars[j]),collapse=""))
  
  j <- 6
  intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Age_Cat)+factor(L_School)+factor(Pa_Index)+Bmi_C+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Centr_C),na.action=na.exclude)
  waldres <- wald.test(b = coef(intmod_mv)[18:19], Sigma = vcov(intmod_mv)[c(18:19),c(18:19)], Terms=1:2)  
  res1 <- summary(intmod_mv)
  res1$aodwald <- waldres$result
  writeLines(toJSON(res1),paste(c("male__smoking__int_mod__age__",clsb[i],"__",outputvars[j]),collapse=""))
  
}
stopCluster(cl)
gc()

#### Interaction Models with BMI#####
Bmi_cat <- Bmi_C
Bmi_cat[which(Bmi_C<30&Bmi_C >25)] <- 1
Bmi_cat[which(Bmi_C>=30)] <- 2
Bmi_cat[which(Bmi_C<=25)] <- 0

cl<-makeCluster(6)
registerDoSNOW(cl)
foreach (i = 1:length(clsb)) %dopar% {  
  library(RJSONIO)
  library(survival)
  library(aod)
  outputvars <- c("Caseclrt_Crs_Frst","Casecl_Crs_Frst","Casert_Crs_Frst","Colon_Prox_Frst","Colon_Dist_Frst","Colon_Nos_Frst")
  qvar <- scan(clsb[i],sep="\t")
  qvar[is.na(qvar)] <- 0
  qvar <- round(log((qvar+0.0000001),base=2),digit=10)
  
    j <- 1
    intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Bmi_cat)+factor(L_School)+factor(Pa_Index)+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Len1,Centr_C),na.action=na.exclude)
    waldres <- wald.test(b = coef(intmod_mv)[17:18], Sigma = vcov(intmod_mv)[c(17:18),c(17:18)], Terms=1:2)  
    res1 <- summary(intmod_mv)
    res1$aodwald <- waldres$result
    writeLines(toJSON(res1),paste(c("male__smoking__int_mod__bmi__",clsb[i],"__",outputvars[j]),collapse=""))
  
  j <- 2
  intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Bmi_cat)+factor(L_School)+factor(Pa_Index)+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Len1,Centr_C),na.action=na.exclude)
  waldres <- wald.test(b = coef(intmod_mv)[17:18], Sigma = vcov(intmod_mv)[c(17:18),c(17:18)], Terms=1:2)  
  res1 <- summary(intmod_mv)
  res1$aodwald <- waldres$result
  writeLines(toJSON(res1),paste(c("male__smoking__int_mod__bmi__",clsb[i],"__",outputvars[j]),collapse=""))
  
  j <- 3
  intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Bmi_cat)+factor(L_School)+factor(Pa_Index)+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Len1,Centr_C),na.action=na.exclude)
  waldres <- wald.test(b = coef(intmod_mv)[17:18], Sigma = vcov(intmod_mv)[c(17:18),c(17:18)], Terms=1:2)  
  res1 <- summary(intmod_mv)
  res1$aodwald <- waldres$result
  writeLines(toJSON(res1),paste(c("male__smoking__int_mod__bmi__",clsb[i],"__",outputvars[j]),collapse=""))
  
  j <- 4
  intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Bmi_cat)+factor(L_School)+factor(Pa_Index)+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Len1,Centr_C),na.action=na.exclude)
  waldres <- wald.test(b = coef(intmod_mv)[17:18], Sigma = vcov(intmod_mv)[c(17:18),c(17:18)], Terms=1:2)  
  res1 <- summary(intmod_mv)
  res1$aodwald <- waldres$result
  writeLines(toJSON(res1),paste(c("male__smoking__int_mod__bmi__",clsb[i],"__",outputvars[j]),collapse=""))
  
  j <- 5
  intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Bmi_cat)+factor(L_School)+factor(Pa_Index)+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Len1,Centr_C),na.action=na.exclude)
  waldres <- wald.test(b = coef(intmod_mv)[17:18], Sigma = vcov(intmod_mv)[c(17:18),c(17:18)], Terms=1:2)  
  res1 <- summary(intmod_mv)
  res1$aodwald <- waldres$result
  writeLines(toJSON(res1),paste(c("male__smoking__int_mod__bmi__",clsb[i],"__",outputvars[j]),collapse=""))
  
  j <- 6
  intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Bmi_cat)+factor(L_School)+factor(Pa_Index)+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Len1,Centr_C),na.action=na.exclude)
  waldres <- wald.test(b = coef(intmod_mv)[17:18], Sigma = vcov(intmod_mv)[c(17:18),c(17:18)], Terms=1:2)  
  res1 <- summary(intmod_mv)
  res1$aodwald <- waldres$result
  writeLines(toJSON(res1),paste(c("male__smoking__int_mod__bmi__",clsb[i],"__",outputvars[j]),collapse=""))

}
stopCluster(cl)
gc()

#### Interaction Models with Alcohol#####  
Alc_cat <- Alc_Re
Alc_cat[which(Alc_Re<30&Alc_Re >5)] <- 1
Alc_cat[which(Alc_Re>=30)] <- 2
Alc_cat[which(Alc_Re<5)] <- 0


cl<-makeCluster(6)
registerDoSNOW(cl)
foreach (i = 1:length(clsb)) %dopar% {  
  library(RJSONIO)
  library(survival)
  library(aod)
  outputvars <- c("Caseclrt_Crs_Frst","Casecl_Crs_Frst","Casert_Crs_Frst","Colon_Prox_Frst","Colon_Dist_Frst","Colon_Nos_Frst")
  
  qvar <- scan(clsb[i],sep="\t")
  qvar[is.na(qvar)] <- 0
  qvar <- round(log((qvar+0.0000001),base=2),digit=10)
  
  for (j in 1:6) { ## this loop is causing big memory problem
    intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Alc_cat)+factor(L_School)+factor(Pa_Index)+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Len1,Centr_C),na.action=na.exclude)
    waldres <- wald.test(b = coef(intmod_mv)[20:21], Sigma = vcov(intmod_mv)[c(20:21),c(20:21)], Terms=1:2)  
    res1 <- summary(intmod_mv)
    res1$aodwald <- waldres$result
    writeLines(toJSON(res1),paste(c("male__smoking__int_mod_alc__",clsb[i],"__",outputvars[j]),collapse=""))
  }
}
stopCluster(cl)
gc()

#### Interaction Models with Fib_Intake#####
Fibt_Cat <- c()
var1 <- QE_Fibt
q0ind <- which(var1 ==0)
quants <- c()
if (length(q0ind)>=(length(var1)/5)) {  # in case of more than 20% zeros, subset the vector, and then do quantile on child. and all zeros belong to Q1. 
  var1_ft <- var1[-q0ind]
  quants <- quantile(var1_ft,probs=c(0.20,0.40,0.60,0.80,1.00))  
} else {
  quants <- quantile(var1,probs=c(0.20,0.40,0.60,0.80,1.00))  
}
q1ind <- which(var1 <= quants[1])
q2ind <- which(var1 > (quants[1])& var1 < quants[2])
q3ind <- which(var1 >= quants[2] & var1 < quants[3])
q4ind <- which(var1 >= quants[3] & var1 < quants[4])
q5ind <- which(var1 >= quants[4] & var1 <= quants[5])
Fibt_Cat[q1ind] <- 0
Fibt_Cat[q2ind] <- 1
Fibt_Cat[q3ind] <- 2
Fibt_Cat[q4ind] <- 3
Fibt_Cat[q5ind] <- 4
Fibt_Cat[q0ind] <- 0

cl<-makeCluster(6)
registerDoSNOW(cl)
foreach (i = 1:length(clsb)) %dopar% {  
  library(RJSONIO)
  library(survival)
  library(aod)
  outputvars <- c("Caseclrt_Crs_Frst","Casecl_Crs_Frst","Casert_Crs_Frst","Colon_Prox_Frst","Colon_Dist_Frst","Colon_Nos_Frst")
  qvar <- scan(clsb[i],sep="\t")
  qvar[is.na(qvar)] <- 0
  qvar <- round(log((qvar+0.0000001),base=2),digit=10)
  
  for (j in 1:6) { 
    intmod_mv <- coxph(Surv(Age_Recr,Agexit_Frst,scan(outputvars[j],sep="\t"))~qvar*factor(Fibt_Cat)+factor(L_School)+factor(Pa_Index)+Bmi_C+meat_red_pr+QE_Energy+QE_CA+Alc_Re+QE_Fibt+strata(Len1,Centr_C),na.action=na.exclude)
    waldres <- wald.test(b = coef(intmod_mv)[23:26], Sigma = vcov(intmod_mv)[c(23:26),c(23:26)], Terms=1:4)  
    res1 <- summary(intmod_mv)
    res1$aodwald <- waldres$result
    writeLines(toJSON(res1),paste(c("male__smoking__int_mod__fibt__",clsb[i],"_",outputvars[j]),collapse=""))
  }
}
stopCluster(cl)
gc()
