
EPIIEd.tmle = function( data,
                        learners = sl3::Lrnr_glm_fast$new(),
                        uvlearners = sl3::Lrnr_hal9001$new(max_degree=5),
                        V=c('W1','W2'),
                        nfolds=5, 
                        family.outcome){

  g_learners<-h_learners<-b_learners<-q_learners<-r_learners<-learners
  u_learners<-v_learners<-uvlearners

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

# bound outcome Y in unit interval
  min_y <- min(data[["Y"]])
  max_y <- max(data[["Y"]])
  if(min_y<0 | max_y>1){
  data.table::set(data, j = "Y", value = scale_to_unit(data[["Y"]]))
}

#get CV eifpiie
eif11<-est_onestep(data=data, contrast=list(rep(1,nrow(data)),rep(1, nrow(data))), g_learners=g_learners, h_learners=h_learners, b_learners=b_learners, q_learners=q_learners, r_learners=r_learners, u_learners=u_learners, v_learners=v_learners, w_names=w_names, m_names=m_names, y_bounds=c(min_y, max_y), ext_weights=weights, cv_folds=5)
eif10<-est_onestep(data=data, contrast=list(rep(1,nrow(data)),rep(0, nrow(data))), g_learners=g_learners, h_learners=h_learners, b_learners=b_learners, q_learners=q_learners, r_learners=r_learners, u_learners=u_learners, v_learners=v_learners, w_names=w_names, m_names=m_names, y_bounds=c(min_y, max_y), ext_weights=weights, cv_folds=5)
eif00<-est_onestep(data=data, contrast=list(rep(0,nrow(data)),rep(0, nrow(data))), g_learners=g_learners, h_learners=h_learners, b_learners=b_learners, q_learners=q_learners, r_learners=r_learners, u_learners=u_learners, v_learners=v_learners, w_names=w_names, m_names=m_names, y_bounds=c(min_y, max_y), ext_weights=weights, cv_folds=5)

eifpiie<-eif11$eif-eif10$eif

eifpite<-eif11$eif-eif00$eif

#add to data table as new variable named D1
data<-data %>% 
  mutate(D1 = eifpiie)

SL.fit.D1<-D_eif(data=data, uvlearners=v_learners, w_names=w_names, cv_folds=5)

d.est.Dpi = as.numeric(SL.fit.D1<0)
  
#get CV eifdd - CV eifd0
eifdd<-est_onestep(data=data, contrast=list(d.est.Dpi ,d.est.Dpi ), g_learners=g_learners, h_learners=h_learners, b_learners=b_learners, q_learners=q_learners, r_learners=r_learners, u_learners=u_learners, v_learners=v_learners, w_names=w_names, m_names=m_names, y_bounds=c(min_y, max_y), ext_weights=weights, cv_folds=5)
eifd0<-est_onestep(data=data, contrast=list(d.est.Dpi ,rep(0,nrow(data))), g_learners=g_learners, h_learners=h_learners, b_learners=b_learners, q_learners=q_learners, r_learners=r_learners, u_learners=u_learners, v_learners=v_learners, w_names=w_names, m_names=m_names, y_bounds=c(min_y, max_y), ext_weights=weights, cv_folds=5)

eifpiiecontrastdest<- eifdd$eif- eifd0$eif
  	#difference between indirect effect with optimal rule and normal indirect effect
eifoptusualcontrast<-eifpiiecontrastdest - eifpiie

	os.est = stats::weighted.mean(eifpiiecontrastdest, weights)
	ci.widths = c(qnorm(0.975)*sd(eifpiiecontrastdest*weights)/sqrt(nrow(data)))
	#estimate EPIIE
	os.est.typ = stats::weighted.mean(eifpiie, weights)
	ci.widths.typ = c(qnorm(0.975)*sd(eifpiie*weights)/sqrt(nrow(data)))

	#estimate contrast
	contrast.est<-stats::weighted.mean(eifoptusualcontrast, weights)
	ci.widths.cont = c(qnorm(0.975)*sd(eifoptusualcontrast*weights)/sqrt(nrow(data)))

#total
totaldest<-stats::weighted.mean(eifpite, weights)
ci.widths.tot<-c(qnorm(0.975)*sd(eifpite*weights)/sqrt(nrow(data)))

totalddest<-stats::weighted.mean(eifdd$eif- eif00$eif, weights)
ci.widths.dtot<-c(qnorm(0.975)*sd((eifdd$eif- eif00$eif)*weights)/sqrt(nrow(data)))

totaldcontrastest<-stats::weighted.mean(eifdd$eif- eif11$eif, weights)
ci.widths.contrasttot<-c(qnorm(0.975)*sd((eifdd$eif- eif11$eif)*weights)/sqrt(nrow(data)))


return(list(EPIIEd.est=os.est,ci.widths=ci.widths,EPIIE.est=os.est.typ,ci.widths.typ=ci.widths.typ,Econtrast.est=contrast.est,contrast.ci.widths=ci.widths.cont, totalddest=totalddest, ci.widths.dtot=ci.widths.dtot, totaldest=totaldest, ci.widths.tot=ci.widths.tot, totaldcontrastest=totaldcontrastest, ci.widths.contrasttot=ci.widths.contrasttot, eif11= eif11$eif,  eif10=eif10$eif,  eif00=eif00$eif,  eifdd=eifdd$eif,  eifd0=eifd0$eif ))

	
#return(list(EPIIEd.est=os.est,ci.widths=ci.widths,EPIIE.est=os.est.typ,ci.widths.typ=ci.widths.typ,Econtrast.est=contrast.est,contrast.ci.widths=ci.widths.cont, PITE = eifpite))

}

EPIIE.tmle = function( data,
                        learners = sl3::Lrnr_glm_fast$new(),
                        uvlearners = sl3::Lrnr_hal9001$new(max_degree=5),
                        nfolds=5, 
                        family.outcome){
  
  g_learners<-h_learners<-b_learners<-q_learners<-r_learners<-learners
  u_learners<-v_learners<-uvlearners
  
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
  
  # bound outcome Y in unit interval
  min_y <- min(data[["Y"]])
  max_y <- max(data[["Y"]])
  if(min_y<0 | max_y>1){
    data.table::set(data, j = "Y", value = scale_to_unit(data[["Y"]]))
  }
  
  #get CV eifpiie
  eif11<-est_onestep(data=data, contrast=list(rep(1,nrow(data)),rep(1, nrow(data))), g_learners=g_learners, h_learners=h_learners, b_learners=b_learners, q_learners=q_learners, r_learners=r_learners, u_learners=u_learners, v_learners=v_learners, w_names=w_names, m_names=m_names, y_bounds=c(min_y, max_y), ext_weights=weights, cv_folds=5)
  eif10<-est_onestep(data=data, contrast=list(rep(1,nrow(data)),rep(0, nrow(data))), g_learners=g_learners, h_learners=h_learners, b_learners=b_learners, q_learners=q_learners, r_learners=r_learners, u_learners=u_learners, v_learners=v_learners, w_names=w_names, m_names=m_names, y_bounds=c(min_y, max_y), ext_weights=weights, cv_folds=5)
  eif00<-est_onestep(data=data, contrast=list(rep(0,nrow(data)),rep(0, nrow(data))), g_learners=g_learners, h_learners=h_learners, b_learners=b_learners, q_learners=q_learners, r_learners=r_learners, u_learners=u_learners, v_learners=v_learners, w_names=w_names, m_names=m_names, y_bounds=c(min_y, max_y), ext_weights=weights, cv_folds=5)
  
  eifpiie<-eif11$eif-eif10$eif
  
  eifpite<-eif11$eif-eif00$eif
  
  os.est.typ = stats::weighted.mean(eifpiie, weights)
  ci.widths.typ = c(qnorm(0.975)*sd(eifpiie*weights)/sqrt(nrow(data)))
    
  #total
  totaldest<-stats::weighted.mean(eifpite, weights)
  ci.widths.tot<-c(qnorm(0.975)*sd(eifpite*weights)/sqrt(nrow(data)))
  
  return(list(EPIIE.est=os.est.typ,ci.widths.typ=ci.widths.typ, totaldest=totaldest, ci.widths.tot=ci.widths.tot))
  
  
  #return(list(EPIIEd.est=os.est,ci.widths=ci.widths,EPIIE.est=os.est.typ,ci.widths.typ=ci.widths.typ,Econtrast.est=contrast.est,contrast.ci.widths=ci.widths.cont, PITE = eifpite))
  
}
