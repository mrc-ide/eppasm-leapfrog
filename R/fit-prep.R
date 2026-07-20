#' Prepare fp object for fitting and prepare likelihood data
#'
#' @param obj specfp or eppfp object
#' @param ... Additional parameters to override on `obj`
#' @param epp If TRUE prep EPP model fit
#'
#' @returns List containing `fp` object ready for fitting and likelihood
#'   data `likdat`
#' @export
prep_fp_fitmod <- function(obj, ..., epp = FALSE) {
  ## ... : updates to fixed parameters (fp) object to specify fitting options

  if (epp) {
    fp <- stats::update(attr(obj, 'eppfp'), ...)
  } else {
    fp <- stats::update(attr(obj, 'specfp'), ...)
  }

  ## Prepare likelihood data
  eppd <- attr(obj, "eppd")

  has_ancrtsite <- exists("ancsitedat", eppd) && any(eppd$ancsitedat$type == "ancrt")
  has_ancrtcens <- !is.null(eppd$ancrtcens) && nrow(eppd$ancrtcens)

  if (!has_ancrtsite) {
    fp$ancrtsite.beta <- 0
  }

  if (has_ancrtsite & has_ancrtcens) {
    fp$ancrt <- "both"
  } else if (has_ancrtsite & !has_ancrtcens) {
    fp$ancrt <- "site"
  } else if (!has_ancrtsite & has_ancrtcens) {
    fp$ancrt <- "census"
  } else {
    fp$ancrt <- "none"
  }

  likdat <- prepare_likdat(eppd, fp)
  fp$ancsitedata <- as.logical(nrow(likdat$ancsite.dat$df))

  if(fp$eppmod %in% c("logrw", "rhybrid")) { # THIS IS REALLY MESSY, NEED TO REFACTOR CODE

    fp$SIM_YEARS <- as.integer(max(likdat$ancsite.dat$df$yidx,
                                   likdat$hhs.dat$yidx,
                                   likdat$ancrtcens.dat$yidx,
                                   likdat$hhsincid.dat$idx))

    fp$proj.steps <- seq(fp$ss$proj_start+0.5, fp$ss$proj_start-1+fp$SIM_YEARS+0.5, by=1/fp$ss$hiv_steps_per_year)
  } else {
    fp$SIM_YEARS <- fp$ss$PROJ_YEARS
  }

  ## Prepare the EPP model
  tsEpidemicStart <- if(epp) fp$tsEpidemicStart else fp$ss$time_epi_start+0.5
  if(fp$eppmod == "rspline")
    fp <- prepare_rspline_model(fp, tsEpidemicStart=tsEpidemicStart)
  else if(fp$eppmod == "rtrend")
    fp <- prepare_rtrend_model(fp)
  else if(fp$eppmod == "logrw")
    fp <- prepare_logrw(fp)
  else if(fp$eppmod == "rhybrid")
    fp <- prepare_rhybrid(fp)
  else if(fp$eppmod == "rlogistic")
    fp$tsEpidemicStart <- fp$proj.steps[which.min(abs(fp$proj.steps - fp$ss$time_epi_start+0.5))]

  fp$logitiota <- TRUE

  ## Prepare the incidence model
  fp$incidmod <- "eppspectrum"

  list(
    fp = fp,
    likdat = likdat
  )
}

#' Prepare leapfrog inputs object for fitting and prepare likelihood data
#'
#' @param params Leapfrog input data
#' @param eppd EPP input data
#' @param ... Additional parameters to set/override on `params`
#'
#' @returns List containing `fp` object ready for fitting and likelihood
#'   data `likdat`
#' @export
prep_fp_fitmod_lf <- function(params, eppd, epp_t0, ...) {

  params <- modifyList(params, list(...))

  # Get number of years of projection from input data. Using survival rate
  # for this but it could be any of the inputs.
  # proj_years is how many years of input data we have
  # sim_years is how long the simulation will run, when fitting with
  # transmission input (i.e. not direct incidence) sim_years may be
  # less than number of proj_years. In this case, we simulate for sim_years
  # and then project forward to fill in up to proj_years.
  ss <- leapfrog::get_leapfrog_ss(LEAPFROG_MODEL_CONFIG)
  params$proj_years <- as.integer(dim(params$Sx)[length(dim(params$Sx))])
  params$sim_years <- params$proj_years
  params$projection_end_year <- params$projection_start_year + params$proj_years - 1

  params$ss$hiv_steps_per_year <- ss$HIV_STEPS_PER_YEAR
  params$proj_steps <- get_projection_steps(params)
  params$ss$AGE_START <- ss$p_idx_hiv_first_adult
  params$ss$f_idx <- 2 # Index to use for women
  params$ss$p_fert_idx <- seq_len(ss$p_fertility_age_groups)
  params$ss$h_ag_span <- ss$hAG_span
  params$ss$hAG <- ss$hAG
  params$ss$ag_idx <- rep(1:params$ss$hAG, params$ss$h_ag_span)
  params$ss$h_fert_idx <- which((params$ss$AGE_START - 1 + cumsum(params$ss$h_ag_span)) %in% 15:49)

  # leapfrog::process_pjnz returns fert_rat and frr_art6mos already mapped onto
  # the coarse HIV fertility age groups (one entry per h_fert_idx), so unlike the
  # Spectrum inputs no age-group remapping is needed here.
  n_h_fert <- length(params$ss$h_fert_idx)
  fert_rat <- params$fert_rat[, seq_len(params$proj_years), drop = FALSE]

  params$frr_cd4 <- array(1, c(ss$hDS, n_h_fert, params$proj_years))
  params$frr_cd4[,,] <- rep(fert_rat, each=ss$hDS)
  params$frr_cd4 <- sweep(params$frr_cd4, 1, params$cd4fert_rat, "*")
  params$frr_cd4 <- params$frr_cd4 * params$frr_scalar

  params$frr_art <- array(1.0, c(ss$hTS, ss$hDS, n_h_fert, params$proj_years))
  params$frr_art[1,,,] <- params$frr_cd4 # 0-6 months
  params$frr_art[2:ss$hTS, , , ] <- sweep(params$frr_art[2:ss$hTS, , , ], 3, params$frr_art6mos * params$frr_scalar, "*") # 6-12mos, >1 years

  has_ancrtsite <- exists("ancsitedat", eppd) && any(eppd$ancsitedat$type == "ancrt")
  has_ancrtcens <- !is.null(eppd$ancrtcens) && (nrow(eppd$ancrtcens) > 0)

  if (!has_ancrtsite) {
    params$ancrtsite_beta <- 0
  }

  if (has_ancrtsite & has_ancrtcens) {
    params$ancrt <- "both"
  } else if (has_ancrtsite & !has_ancrtcens) {
    params$ancrt <- "site"
  } else if (!has_ancrtsite & has_ancrtcens) {
    params$ancrt <- "census"
  } else {
    params$ancrt <- "none"
  }

  likdat <- prepare_likdat_lf(params, eppd)
  params$ancsitedata <- as.logical(nrow(likdat$ancsite_dat$df))

  if (params$eppmod %in% c("logrw", "rhybrid")) {

    params$sim_years <- as.integer(max(likdat$ancsite_dat$df$yidx,
                                       likdat$hhs_dat$yidx,
                                       likdat$ancrtcens_dat$yidx,
                                       likdat$hhsincid_dat$idx))

    params$proj_steps <- get_projection_steps(params)
  }

  ## Prepare the EPP model
  params$ts_epidemic_start <- epp_t0 + 0.5
  if (params$eppmod == "rspline") {
    fp <- prepare_rspline_model_lf(params,
                                   ts_epidemic_start = params$ts_epidemic_start)
  } else if (params$eppmod == "logrw") {
    fp <- prepare_logrw_lf(params,
                           ts_epidemic_start = params$ts_epidemic_start)
  } else if (params$eppmod == "rhybrid") {
    fp <- prepare_rhybrid_lf(params)
  } else if (params$eppmod == "rlogistic") {
    fp <- params
    fp$ts_epidemic_start <-
      params$proj_steps[
        which.min(abs(params$proj_steps - params$ts_epidemic_start))
      ]
  }

  fp$logitiota <- TRUE

  ## Prepare the incidence model
  fp$incidmod <- "eppspectrum"

  list(
    fp = fp,
    likdat = likdat
  )
}

get_projection_steps <- function(params) {
  params$projection_start_year + 0.5 + 0:(params$ss$hiv_steps_per_year * (params$sim_years - 1)) / params$ss$hiv_steps_per_year
}
