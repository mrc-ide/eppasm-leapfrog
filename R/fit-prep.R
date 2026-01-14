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
prep_fp_fitmod_lf <- function(params, eppd, ...) {

  params <- modifyList(params, list(...))

  ## TODO: Get hiv steps per year and age start from leapfrog state space
  hiv_steps_per_year <- 10
  params$ss$hiv_steps_per_year <- hiv_steps_per_year
  params$ss$AGE_START <- 15
  params$sim_years <- params$proj_years
  params$proj_steps <- get_projection_steps(params)

  has_ancrtsite <- exists("ancsitedat", eppd) && any(eppd$ancsitedat$type == "ancrt")
  has_ancrtcens <- !is.null(eppd$ancrtcens) && (nrow(eppd$ancrtcens) > 0)

  ## TODO: rename
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

  if(params$eppmod %in% c("logrw", "rhybrid")) {

    params$sim_years <- as.integer(max(likdat$ancsite_dat$df$yidx,
                                       likdat$hhs_dat$yidx,
                                       likdat$ancrtcens_dat$yidx,
                                       likdat$hhsincid_dat$idx))

    params$proj_steps <- get_projection_steps(params)
  }

  ## Prepare the EPP model
  ## TODO: what should we do with ts_epidemic_start, this isn't really in leapfrog
  ## just hardcoding for now
  ts_epidemic_start <- 1975.5
  if(params$eppmod == "rspline") {
    fp <- prepare_rspline_model_lf(params, ts_epidemic_start = ts_epidemic_start)
  } else if(params$eppmod == "logrw") {
    stop("TODO impl logrw for leapfrog")
    fp <- prepare_logrw(params)
  } else if(params$eppmod == "rhybrid") {
    fp <- prepare_rhybrid_lf(params)
  } else if(params$eppmod == "rlogistic") {
    fp <- params
    fp$ts_epidemic_start <- params$proj_steps[which.min(abs(params$proj_steps - ts_epidemic_start + 0.5))]
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
