cv_eif <- function(fold,
                   data_in,
                   #contrast,
                   g_learners,
                   h_learners,
                   b_learners,
                   q_learners,
                   r_learners,
                   u_learners,
                   v_learners,
                   w_names,
                   m_names) {
  # make training and validation data
  train_data <- origami::training(data_in)
  valid_data <- origami::validation(data_in)

  # 1) fit regression for propensity score regression
  g_out <- fit_treat_mech(
    train_data = train_data,
    valid_data = valid_data,
    #contrast = contrast,
    learners = g_learners,
    w_names = w_names,
    m_names = m_names,
    type = "g"
  )

  # 2) fit clever regression for treatment, conditional on mediators
  h_out <- fit_treat_mech(
    train_data = train_data,
    valid_data = valid_data,
    #contrast = contrast,
    learners = h_learners,
    w_names = w_names,
    m_names = m_names,
    type = "h"
  )

  # 3) fit outcome regression
  b_out <- fit_out_mech(
    train_data = train_data,
    valid_data = valid_data,
    #contrast = contrast,
    learners = b_learners,
    m_names = m_names,
    w_names = w_names
  )

  # 4) fit mediator-outcome confounder regression, excluding mediator(s)
  q_out <- fit_moc_mech(
    train_data = train_data,
    valid_data = valid_data,
    #contrast = contrast,
    learners = q_learners,
    m_names = m_names,
    w_names = w_names,
    type = "q"
  )

  # 5) fit mediator-outcome confounder regression, conditioning on mediator(s)
  r_out <- fit_moc_mech(
    train_data = train_data,
    valid_data = valid_data,
    #contrast = contrast,
    learners = r_learners,
    m_names = m_names,
    w_names = w_names,
    type = "r"
  )

  # extract components; NOTE: only do this for observations in validation set
  b_prime <- b_out$b_est_valid$b_pred_A_prime
  h_star <- h_out$treat_est_valid$treat_pred_A_star
  g_star <- g_out$treat_est_valid$treat_pred_A_star
  h_prime <- h_out$treat_est_valid$treat_pred_A_prime
  g_prime <- g_out$treat_est_valid$treat_pred_A_prime
  q_prime_Z_one <- q_out$moc_est_valid_Z_one$moc_pred_A_prime
  r_prime_Z_one <- r_out$moc_est_valid_Z_one$moc_pred_A_prime
  q_prime_Z_natural <- q_out$moc_est_valid_Z_natural$moc_pred_A_prime
  r_prime_Z_natural <- r_out$moc_est_valid_Z_natural$moc_pred_A_prime

  # need pseudo-outcome regressions with intervention set to a contrast
  # NOTE: training fits of these nuisance functions must be performed using the
  #       data corresponding to the natural intervention value but predictions
  #       are only needed for u(z,a',w) and v(a*,w) as per form of the EIF
  valid_data_a_prime <- data.table::copy(valid_data)[, A := C1]
  valid_data_a_star <- data.table::copy(valid_data)[, A := C2]
  u_out <- fit_nuisance_u(
    train_data = train_data,
    valid_data = valid_data_a_prime,
    learners = u_learners,
    b_out = b_out,
    q_out = q_out,
    r_out = r_out,
    g_out = g_out,
    h_out = h_out,
    w_names = w_names
  )
  u_prime <- u_out$u_pred

  v_out <- fit_nuisance_v(
    train_data = train_data,
    valid_data = valid_data_a_star,
    #contrast = contrast,
    learners = v_learners,
    b_out = b_out,
    q_out = q_out,
    m_names = m_names,
    w_names = w_names
  )
  v_star <- v_out$v_pred

  # NOTE: assuming Z in {0,1}, other cases not supported yet
  u_int_eif <- lapply(c(1, 0), function(z_val) {
    # intervene on training and validation data sets
    valid_data_z_interv <- data.table::copy(valid_data)
    C1<-valid_data_z_interv[,"C1"]
    valid_data_z_interv[, `:=`(
      Z = z_val,
      A = C1,
      U_pseudo = u_prime
    )]

    # predict u(z, a', w) using intervened data with treatment set A = a'
    u_task_valid_z_interv <- sl3::sl3_Task$new(
      data = valid_data_z_interv,
      weights = "obs_weights",
      covariates = c("Z", "A", w_names),
      outcome = "U_pseudo",
      outcome_type = "continuous"
    )

    # return partial pseudo-outcome for v nuisance regression
    out_valid <- u_out[["u_fit"]]$predict(u_task_valid_z_interv)
    return(out_valid)
  })
  u_int_eif <- do.call(`-`, u_int_eif)

  # create inverse probability weights
  C1<-valid_data$C1
  C2<-valid_data$C2
  ipw_a_prime <- as.numeric(valid_data$A == C1) / g_prime
  ipw_a_star <- as.numeric(valid_data$A == C2) / g_star

  # residual term for outcome component of EIF
  c_star <- (g_prime / g_star) * (q_prime_Z_natural / r_prime_Z_natural) *
    (h_star / h_prime)

  # compute uncentered efficient influence function components
  eif_y <- ipw_a_prime * c_star / mean(ipw_a_prime * c_star) *
    (valid_data$Y - b_prime)
  eif_u <- ipw_a_prime / mean(ipw_a_prime) * u_int_eif *
    (valid_data$Z - q_prime_Z_one)
  eif_v <- ipw_a_star / mean(ipw_a_star) * (v_out$v_pseudo - v_star)

  # un-centered efficient influence function
  eif <- eif_y + eif_u + eif_v + v_star

  # output list
  out <- list(data.table::data.table(
    # components necessary for fluctuation step of TMLE
    g_prime = g_prime, g_star = g_star, h_prime = h_prime, h_star = h_star,
    q_prime_Z_natural = q_prime_Z_natural, q_prime_Z_one = q_prime_Z_one,
    r_prime_Z_natural = r_prime_Z_natural, r_prime_Z_one = r_prime_Z_one,
    v_star = v_star, u_int_diff = u_int_eif,
    b_prime = b_prime, b_prime_Z_zero = v_out$b_A_prime_Z_zero,
    b_prime_Z_one = v_out$b_A_prime_Z_one,
    # efficient influence function and fold IDs
    D_star = eif, fold = origami::fold_index()
  ))
  return(out)
}

D_eif<-function(data, 
                uvlearners, 
                w_names, 
                cv_folds=5){

  # make sure that more than one fold is specified
  assertthat::assert_that(cv_folds > 1)

  # create cross-validation folds
  folds <- origami::make_folds(data,
    fold_fun = origami::folds_vfold,
    V = cv_folds
  )

# estimate the EIF on a per-fold basis
  cv_d_eif_results <- origami::cross_validate(
    cv_fun = cv_d_eif,
    folds = folds,
    data_in = data,
    uvlearners = uvlearners,
    w_names = w_names,
    use_future = FALSE,
    .combine = FALSE
  )

  # get estimated efficient influence function
  cv_d_eif_est <- do.call(c, lapply(cv_d_eif_results[[1]], `[[`, "D_star"))
  obs_valid_idx <- do.call(c, lapply(folds, `[[`, "validation_set"))
  cv_d_eif_est <- cv_d_eif_est[order(obs_valid_idx)]

  d_eif_est_out<-cv_d_eif_est

  return(d_eif_est_out)


}

cv_d_eif <- function(fold,
                   data_in,
                   uvlearners,
                   w_names) {
  
  # make training and validation data
  train_data <- origami::training(data_in)
  valid_data <- origami::validation(data_in)


 ## construct task for treatment mechanism fit

  #ker nov2020 the bug is in the treat task portion of hte fit_rule_mech
 # treat_task <- sl3::sl3_Task$new(
 #   data = train_data,
 #   weights = "obs_weights",
 #   covariates = cov_names,
 #   outcome = "D1",
 #   outcome_type = "continuous"
 # )

d_out <- fit_rule_mech(
    train_data = train_data,
    valid_data = valid_data,
    learners = uvlearners,
    w_names = w_names
  )

  # extract components; NOTE: only do this for observations in validation set
  d_est <- d_out$treat_est_valid$treat_pred

  # output list
  out <- list(data.table::data.table(
    # components necessary for fluctuation step of TMLE
        D_star = d_est, fold = origami::fold_index()
  ))
  return(out)
}

est_onestep <- function(data,
                        contrast,
                        g_learners,
                        h_learners,
                        b_learners,
                        q_learners,
                        r_learners,
                        u_learners,
                        v_learners,
                        w_names,
                        m_names,
                        y_bounds,
                        ext_weights = NULL,
                        cv_folds = 5) {

  # make sure that more than one fold is specified
  assertthat::assert_that(cv_folds > 1)

  # create cross-validation folds
  folds <- origami::make_folds(data,
    fold_fun = origami::folds_vfold,
    V = cv_folds
  )

  data<-as.data.table(data %>% 
  mutate(C1 = contrast[[1]], C2=contrast[[2]]))

  # estimate the EIF on a per-fold basis
  cv_eif_results <- origami::cross_validate(
    cv_fun = cv_eif,
    folds = folds,
    data_in = data,
    #contrast = contrast,
    g_learners = g_learners,
    h_learners = h_learners,
    b_learners = b_learners,
    q_learners = q_learners,
    r_learners = r_learners,
    u_learners = u_learners,
    v_learners = v_learners,
    w_names = w_names,
    m_names = m_names,
    use_future = FALSE,
    .combine = FALSE
  )

  # get estimated efficient influence function
  cv_eif_est <- do.call(c, lapply(cv_eif_results[[1]], `[[`, "D_star"))
  obs_valid_idx <- do.call(c, lapply(folds, `[[`, "validation_set"))
  cv_eif_est <- cv_eif_est[order(obs_valid_idx)]

  # re-scale efficient influence function
  eif_est_rescaled <- cv_eif_est 
  #%>%
  #  scale_from_unit(y_bounds[2], y_bounds[1])

  # compute one-step estimate and variance from efficient influence function
  if (is.null(ext_weights)) {
    os_est <- mean(eif_est_rescaled)
    eif_est_out <- eif_est_rescaled
  } else {
    # compute a re-weighted one-step, with re-weighted influence function
    os_est <- stats::weighted.mean(eif_est_rescaled, ext_weights)
    eif_est_out <- eif_est_rescaled * ext_weights
  }
  os_var <- stats::var(eif_est_out) / length(eif_est_out)

  # output
  os_est_out <- list(
    theta = os_est,
    var = os_var,
    eif = eif_est_out,
    type = "onestep"
  )
  return(os_est_out)
}
