setwd("/projects/users/xxx")
library(R6)
library(data.table)
library(SuperLearner)
library(glmnet)
library(data.table)
library(tidyverse)
library(BBmisc)
library(uuid)
library(origami)
library(assertthat)
library(digest)
library(delayed)
library(imputeMissings)
install.packages("sl3-1.3.7.tar.gz", repos=NULL, type="source")
library(sl3)

source("utils.R") 
source("fit_mechanisms.R") 
source("estimators.R")
source("funs_piie_cv.R")

## Load data. Used MTO final analytic dataset

implist <- lapply(Sys.glob("xxxx*.csv"), read.csv)

## data cleaning, organizing

  datimp<-implist[[1]]
  datimp$re0v1<-ifelse(datimp$raceethn==0,1, 0)
  datimp$re2v1<-ifelse(datimp$raceethn==2,1, 0)
  datimp$re3v1<-ifelse(datimp$raceethn==3,1, 0)
  Wnames=c( "RA_SITE" , "re0v1", "re2v1", "re3v1", "F_SVY_AGE_BL_IMP",  "X_F_C1_BEHPRB" , "X_F_C1_GIFTED","X_F_AD_EDGRADHS",  "X_F_AD_NEVMARR","X_F_AD_PARENTU18",
            "X_F_AD_WORKING", "X_F_HH_AFDC" , "X_F_HH_DISABL", "X_F_HH_SIZE2", "X_F_HH_SIZE3", "X_F_HH_SIZE4",
            "X_F_HOOD_UNSAFENIT",  "X_F_HOOD_VERYDISSAT", "X_F_HOUS_MOV3TM","X_F_HOUS_MOVSCHL", "X_F_HOUS_SEC8BEF" , "F_C9010T_PERPOV_BL"  )
  #Wnames=c( "re0v1", "re2v1", "re3v1", "F_SVY_AGE_BL_IMP",  "F_C9010T_PERPOV_BL"  )
  A<-"RA_POOLGRP_EXPS8"
  Z<-"F_SVY_CMOVE.x"
  weights<-"F_WT_TOTSVY"

 ##Neighborhood
  Mnhoodpov<-"F_C9010T_PERPOV_DW"
  
  ##school environment
  Mse1<-"F_SC_S_FRRDLUNCH_DUR"
  Mse2<-"F_SC_TITLE1_ALL_DUR"
  Mse3<-"F_SC_RANK_DUR"
  Mse4<-"F_SC_RATIO_PUPTCH_DUR"
  Mse<-c(Mse1, Mse2, Mse3, Mse4)
  
  ##school/home instability
  Minstab1<-"F_SH_CNT_SCH"
  Minstab2<-"F_SPL_MOVES_N"
  Minstab3<-"F_SH_BLDIST_RCNT"
  Minstab5<-"F_SH_CNT_YR_SCHCHG"
  Minstab<-c(Minstab1, Minstab2, Minstab3,  Minstab5)
  
  allmediators<-c(Mnhoodpov, Mse, Minstab)
    
    A<-"RA_POOLGRP_EXPS8"
    Z<-"F_SVY_CMOVE.x"
    weights<-"F_WT_TOTSVY"

    Y<-"F_RB_ALC_DRNK_M"
    alldat<-datimp[datimp$F_SVY_GENDER==0 & datimp$RA_SITE!=1, c(Wnames,A, Z,  allmediators, Y, weights) ]
    
    setnames(alldat, old = c(A,Z,Y, weights), new = c( "A", "Z", "Y", "weights"))
    
    newWnames<-c()
    for(g in 1:length(Wnames)){
      newWnames[g]<-paste0("W",g)
    }
    
    newMnames<-c()
    for(m in 1:length(allmediators)){
      newMnames[m]<-paste0("M",m)
    }
    setnames(alldat, old=c(Wnames,allmediators), new=c(newWnames, newMnames))
    
    
    A<-alldat[,"A"]
    M<-alldat[,substr(names(alldat), 1,1)=="M"]
    Z<-alldat[,"Z"]
    Y<-alldat[,"Y"]
    W<-alldat[substr(names(alldat),1,1)=="W"]

## Make learners
mean_lrnr<-Lrnr_mean$new()
fglm_lrnr<-Lrnr_glm_fast$new(family=binomial())
fglm_lrnr_uv<-Lrnr_glm_fast$new(family=gaussian())
mars_lrnr<-Lrnr_earth$new()
glmnet_lrnr<-Lrnr_glmnet$new(family=binomial())
xgboost_lrnr<-Lrnr_xgboost$new()

slbin<-Lrnr_sl$new(learners=list(fglm_lrnr,mean_lrnr,xgboost_lrnr, mars_lrnr, glmnet_lrnr), metalearner=Lrnr_nnls$new())
sluv<-Lrnr_sl$new(learners=list(fglm_lrnr_uv, mean_lrnr, xgboost_lrnr, mars_lrnr, glmnet_lrnr), metalearner=Lrnr_nnls$new())

## estimate d(v) and the PIIE, PITE using d(v)
#6 min
esttest<-EPIIEd.tmle(data=alldat,learners = slbin, uvlearners = sluv,V=c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8","W9", "W10", "W11", "W12","W13", "W14", "W15", "W16","W17", "W18", "W19", "W20","W21", "W22"),nfolds=5, family.outcome="binomial")

est<-esttest

#make results table using SuperLearner
piieest<-c(est$EPIIEd.est, est$EPIIE.est, est$Econtrast.est)
piielci<-c(est$EPIIEd.est - est$ci.widths, est$EPIIE.est-est$ci.widths.typ, est$Econtrast.est-est$contrast.ci.widths)
piieuci<-c(est$EPIIEd.est + est$ci.widths, est$EPIIE.est+est$ci.widths.typ, est$Econtrast.est+est$contrast.ci.widths)

piteest<-c(est$totalddest, est$totaldest, est$totaldcontrastest)
pitelci<-c(est$totalddest-est$ci.widths.dtot, est$totaldest - est$ci.widths.tot, est$totaldcontrastest-est$ci.widths.contrasttot)
piteuci<-c(est$totalddest+est$ci.widths.dtot, est$totaldest + est$ci.widths.tot, est$totaldcontrastest+est$ci.widths.contrasttot)

restable<-as.data.frame(rbind(cbind(piieest, piielci, piieuci), cbind(piteest, pitelci, piteuci)))
colnames(restable)<-c("est", "lci", "uci")
restable$effect<-c(rep("indirect", 3), rep("total", 3))
restable$rule<-rep(c("rule", "no rule", "contrast"),2)
restable

#estimation using the adaptive lasso
set.seed(4309)

#lasso approach
alldat$d<-as.numeric(I(est$eif11-est$eif10)<0)
datamatrix<-model.matrix( ~ (W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12+W13+W14+W15+W16+W17+W18+W19+W20+W21+W22)^2, alldat)[,-1]
newdat<-as.data.frame(datamatrix)
newdat$d<-alldat$d
res<-lm(d ~ ., data=newdat)
pfacm<-(1/abs(res$coefficients))
pfacm<-ifelse(is.na(pfacm), 2000, pfacm)

mstarfit<-cv.glmnet(datamatrix, alldat$d, family="binomial", penalty.factor=pfacm, nfolds=5)

coef(mstarfit, s=mstarfit$lambda.min, gamma=c(1), relax=TRUE)
vsecondorderint<-rownames(coef(mstarfit, s="lambda.min"))[coef(mstarfit, s="lambda.min")[,1]!=0]

suboptrulelasso<-predict(mstarfit, s="lambda.min", gamma=c(1), relax=TRUE, newx=datamatrix, type="class")
learners = slbin 
uvlearners = sluv
g_learners<-h_learners<-b_learners<-q_learners<-r_learners<-learners
u_learners<-v_learners<-uvlearners
data<-alldat
A <- data[, "A"]
M <- data[, substr(names(data), 1, 1) == "M"]
Z <- data[, substr(names(data), 1, 1) == "Z"]
Y <- data[, "Y"]
W <- data[, substr(names(data), 1, 1) == "W"]

obs_weights = rep(1, length(Y))
weights<-data[, "weights"]

data <- data.table::as.data.table(cbind(Y, M, Z, A, W, obs_weights))
w_names <- paste("W", seq_len(dim(data.table::as.data.table(W))[2]),
                 sep = "_"
)
m_names <- paste("M", seq_len(dim(data.table::as.data.table(M))[2]),
                 sep = "_"
)
data.table::setnames(data, c("Y", m_names, "Z", "A", w_names, "obs_weights"))

d.est.subopt<-as.numeric(suboptrulelasso)

eifdd<-est_onestep(data=data, contrast=list(d.est.subopt,d.est.subopt), g_learners=g_learners, h_learners=h_learners, b_learners=b_learners, q_learners=q_learners, r_learners=r_learners, u_learners=u_learners, v_learners=v_learners, w_names=w_names, m_names=m_names, y_bounds=c(min_y, max_y), ext_weights=weights, cv_folds=5)
eifd0<-est_onestep(data=data, contrast=list(d.est.subopt ,rep(0,nrow(data))), g_learners=g_learners, h_learners=h_learners, b_learners=b_learners, q_learners=q_learners, r_learners=r_learners, u_learners=u_learners, v_learners=v_learners, w_names=w_names, m_names=m_names, y_bounds=c(min_y, max_y), ext_weights=weights, cv_folds=5)

eifpiiecontrastdest<- eifdd$eif- eifd0$eif
eifoptusualcontrast<-eifpiiecontrastdest - (esttest$eif11 - esttest$eif10)
  
# Estimate EPIIEd
os.est = stats::weighted.mean(eifpiiecontrastdest, weights)
ci.widths = c(qnorm(0.975)*sd(eifpiiecontrastdest*weights)/sqrt(nrow(data)))
  
#estimate contrast
contrast.est<-stats::weighted.mean(eifoptusualcontrast, weights)
ci.widths.cont = c(qnorm(0.975)*sd(eifoptusualcontrast*weights)/sqrt(nrow(data)))
  
#total
totalddest<-stats::weighted.mean(eifdd$eif- esttest$eif00, weights)
ci.widths.dtot<-c(qnorm(0.975)*sd((eifdd$eif- esttest$eif00)* weights)/sqrt(nrow(data)))
  
totaldcontrastest<-stats::weighted.mean(eifdd$eif- esttest$eif11, weights)
ci.widths.contrasttot<-c(qnorm(0.975)*sd((eifdd$eif- esttest$eif11)*weights)/sqrt(nrow(data)))
  
suboptruleest<-list(EPIIEd.est=os.est,ci.widths=ci.widths,Econtrast.est=contrast.est,contrast.ci.widths=ci.widths.cont, totalddest=totalddest, ci.widths.dtot=ci.widths.dtot,totaldcontrastest=totaldcontrastest, ci.widths.contrasttot=ci.widths.contrasttot )
