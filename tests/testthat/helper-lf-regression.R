# Helpers for the regression tests which check that the leapfrog-structured
# likelihood / fitting code (ll_lf, fitmod_lf, ...) gives the same answers as
# the original code (ll, fitmod, ...) when both are handed the *same model
# output*.
#
# Both code paths ultimately call a model-running function, simmod() for the
# original path and simmod_lf() for the leapfrog path. We mock both so that
# they return the same numbers from the EPPASM C model, just laid out the way
# each path expects:
#
#   original path:  simmod(fp)     -> "spec" array with attributes
#   leapfrog path:  simmod_lf(fp)  -> list(pop, h_hivpop, h_artpop,
#                                          prev_15to49, incid_15to49)
#
# Any difference in the output of the two paths is then down to the
# likelihood / prior / fitting code or to the input data conversion
# (prep_fp_fitmod_lf, prepare_likdat_lf etc.), not to the model engine.

# Run the original EPPASM C engine, this is the reference "model".
run_eppasm_model <- function(fp) {
  # The C code writes to fp$birthslag in place. This does not change the model
  # output (every element is written before it is read within a run) but the
  # vector is shared with the cached inputs and every fp derived from them by
  # update() or modifyList(), so it would be silently modifying test fixtures.
  # `+ 0` forces a copy.
  fp$birthslag <- fp$birthslag + 0
  simmod.specfp(fp, VERSION = "C")
}

# Convert the output of EPPASM into the structure returned by simmod_lf(), which
# is what the leapfrog likelihood functions read. Only the elements that the
# likelihood code uses are populated.
convert_to_leapfrog_data <- function(mod) {
  pop <- mod
  attributes(pop) <- list(dim = dim(mod))
  list(
    pop = pop,
    h_hivpop = attr(mod, "hivpop"),
    h_artpop = attr(mod, "artpop"),
    prev_15to49 = attr(mod, "prev15to49"),
    incid_15to49 = attr(mod, "incid15to49")
  )
}

# The parameters set by fnCreateParam_lf() that change the model output (as
# opposed to only changing the likelihood). simmod_lf() is given a leapfrog fp
# so, to run EPPASM on the matching input, these are copied into an
# EPPASM-structured fp. If a new theta-dependent parameter starts to affect the
# model it needs adding here, otherwise the mocked models will silently diverge.
MODEL_PARAMS <- c("rvec", "iota", "frr_cd4", "frr_art")

# Mock simmod() and simmod_lf() for the duration of the calling test.
#
# fp_eppasm is the EPPASM-structured fp (from prep_fp_fitmod()) the model is
# run with, it is a template which gets the theta-dependent parameters from
# the leapfrog fp every time simmod_lf() is called.
local_mock_models <- function(fp_eppasm, .env = parent.frame()) {
  testthat::local_mocked_bindings(
    simmod = function(fp, ...) run_eppasm_model(fp),
    simmod_lf = function(fp) {
      model_params <- unclass(fp)[MODEL_PARAMS]
      # Check here as the C model segfaults, rather than errors, on missing inputs
      missing <- setdiff(MODEL_PARAMS, names(Filter(length, model_params)))
      if (length(missing) > 0) {
        stop("Leapfrog fp is missing or has empty model parameter(s): ",
             paste(missing, collapse = ", "), call. = FALSE)
      }
      fp_matching <- utils::modifyList(fp_eppasm, model_params)
      convert_to_leapfrog_data(run_eppasm_model(fp_matching))
    },
    .package = "eppasm.lf",
    .env = .env
  )
}

# Build the inputs to fitting for the original code path (ll, fitmod) and for
# the leapfrog code path (ll_lf, fitmod_lf) from the same PJNZ.
#
# Returns the output of prep_fp_fitmod() as `eppasm` and of
# prep_fp_fitmod_lf() as `lf`, each containing the `fp` and `likdat`.
# Cached as reading the PJNZ and preparing the splines is slow-ish.
.regression_inputs_cache <- new.env()

build_regression_inputs <- function(eppmod, proj_end = 2025.5, ...) {
  key <- paste(eppmod, proj_end, paste(names(list(...)), list(...), collapse = ","))
  if (!is.null(.regression_inputs_cache[[key]])) {
    return(.regression_inputs_cache[[key]])
  }

  pjnz <- system_file("extdata", "testpjnz", "Botswana2018.PJNZ")

  eppasm_inputs <- suppressMessages(prepare_spec_fit(pjnz, proj_end))
  # do.call as update.specfp() evaluates its arguments in the fp, so a variable
  # called eppmod could not be passed through by name
  eppasm <- do.call(prep_fp_fitmod,
                    c(list(eppasm_inputs$Urban, eppmod = eppmod), list(...)))

  lf_inputs <- suppressMessages(
    leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE))
  eppd <- prep_epp_data(pjnz)
  epp_t0 <- read_epp_t0(pjnz)
  lf <- prep_fp_fitmod_lf(lf_inputs, eppd$Urban, epp_t0["Urban"],
                          eppmod = eppmod, ...)

  out <- list(
    pjnz = pjnz, eppasm = eppasm, lf = lf,
    # for passing to fitmod() and fitmod_lf()
    eppasm_obj = eppasm_inputs$Urban,
    lf_params = lf_inputs, eppd = eppd$Urban, epp_t0 = epp_t0["Urban"]
  )
  .regression_inputs_cache[[key]] <- out
  out
}

# Epidemic models supported by both code paths, with any extra options
# needed to prepare them. rtrend is not supported by leapfrog.
regression_models <- list(
  rhybrid = list(rw_start = 2005),
  rspline = list(),
  logrw = list(),
  rlogistic = list()
)

# Draw thetas from the prior; used to try the likelihood out at some different,
# but plausible, parameter values.
sample_thetas <- function(inputs, n = 5, seed = 1) {
  set.seed(seed)
  sample.prior(n, inputs$eppasm$fp)
}

# The names of parameters differ between the two code paths, e.g. v.infl and
# v_infl. Standardise so they can be compared.
standardise_names <- function(x) {
  names(x) <- gsub(".", "_", names(x), fixed = TRUE)
  x[order(names(x))]
}

# Run an expression with a time limit. fitmod() and fitmod_lf() retry forever
# if the IMIS step errors, so a regression could otherwise hang the test run.
with_time_limit <- function(expr, seconds = 600) {
  setTimeLimit(elapsed = seconds, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf))
  expr
}
