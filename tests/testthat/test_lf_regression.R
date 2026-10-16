# Regression tests: the leapfrog-structured likelihood and fitting code
# (ll_lf, lprior_lf, fitmod_lf ...) should give the same answers as the
# original code (ll, lprior, fitmod ...) when run on the same model output.
#
# The model is mocked in both code paths (see helper-lf-regression.R) so that
# only the code around the model is compared.
#
# Thetas are drawn from the prior; there is one test per epidemic model so that
# a failure in one does not hide the others.

for (eppmod in names(regression_models)) {
  local({
    eppmod <- eppmod
    get_inputs <- function() {
      do.call(build_regression_inputs, c(list(eppmod), regression_models[[eppmod]]))
    }

    test_that(paste0("fnCreateParam_lf gives the same parameters as fnCreateParam: ", eppmod), {
      inputs <- get_inputs()
      thetas <- sample_thetas(inputs, 3)

      for (i in seq_len(nrow(thetas))) {
        expected <- fnCreateParam(thetas[i, ], inputs$eppasm$fp)
        actual <- fnCreateParam_lf(thetas[i, ], inputs$lf$fp)
        expect_equal(standardise_names(actual), standardise_names(expected),
                     info = paste("theta", i))
      }
    })

    test_that(paste0("lprior_lf gives the same value as lprior: ", eppmod), {
      inputs <- get_inputs()
      thetas <- sample_thetas(inputs, 5)

      for (i in seq_len(nrow(thetas))) {
        expect_equal(lprior_lf(thetas[i, ], inputs$lf$fp),
                     lprior(thetas[i, ], inputs$eppasm$fp),
                     info = paste("theta", i))
      }
    })

    test_that(paste0("sample_prior_lf draws the same samples as sample.prior: ", eppmod), {
      inputs <- get_inputs()
      set.seed(2)
      expected <- sample.prior(10, inputs$eppasm$fp)
      set.seed(2)
      actual <- sample_prior_lf(10, inputs$lf$fp)
      expect_equal(actual, expected)
    })

    test_that(paste0("dsamp_lf gives the same density as dsamp: ", eppmod), {
      inputs <- get_inputs()
      thetas <- sample_thetas(inputs, 5)

      expect_equal(dsamp_lf(thetas, inputs$lf$fp, log = TRUE),
                   dsamp(thetas, inputs$eppasm$fp, log = TRUE))
    })

    test_that(paste0("simmod_lf returns what the model mock provides: ", eppmod), {
      # The mocked models are only a valid stand in if the real simmod_lf output
      # has the same structure as the mock, which is checked here by running the
      # real leapfrog model.
      inputs <- get_inputs()
      theta <- sample_thetas(inputs, 1)[1, ]

      fp_lf <- utils::modifyList(inputs$lf$fp, fnCreateParam_lf(theta, inputs$lf$fp))
      fp_eppasm <- update(inputs$eppasm$fp,
                          list = fnCreateParam(theta, inputs$eppasm$fp))

      actual <- simmod_lf(fp_lf)
      mock <- convert_to_leapfrog_data(run_eppasm_model(fp_eppasm))

      expect_true(all(names(mock) %in% names(actual)))
      for (name in c("pop", "h_hivpop", "h_artpop")) {
        expect_equal(dim(actual[[name]]), dim(mock[[name]]), info = name)
      }
      # Unlike the arrays above, simmod_lf does not pad the 15 to 49 vectors out
      # to proj_years
      for (name in c("prev_15to49", "incid_15to49")) {
        expect_equal(length(actual[[name]]), fp_lf$sim_years, info = name)
        expect_true(length(mock[[name]]) >= fp_lf$sim_years, info = name)
      }
    })

    test_that(paste0("ll_lf gives the same likelihood as ll: ", eppmod), {
      inputs <- get_inputs()
      local_mock_models(inputs$eppasm$fp)
      thetas <- sample_thetas(inputs, 5)

      for (i in seq_len(nrow(thetas))) {
        expected <- ll(thetas[i, ], inputs$eppasm$fp, inputs$eppasm$likdat)
        actual <- ll_lf(thetas[i, ], inputs$lf$fp, inputs$lf$likdat)
        expect_equal(actual, expected, info = paste("theta", i))
      }
    })
  })
}

test_that("likelihood_lf and prior_lf give the same values as likelihood and prior", {
  inputs <- build_regression_inputs("rhybrid", rw_start = 2005)
  local_mock_models(inputs$eppasm$fp)
  thetas <- sample_thetas(inputs, 5)

  expect_equal(
    likelihood_lf(thetas, inputs$lf$fp, inputs$lf$likdat, log = TRUE),
    likelihood(thetas, inputs$eppasm$fp, inputs$eppasm$likdat, log = TRUE))
  expect_equal(prior_lf(thetas, inputs$lf$fp, log = TRUE),
               prior(thetas, inputs$eppasm$fp, log = TRUE))
})

test_that("fitmod_lf with optfit gives the same fit as fitmod", {
  inputs <- build_regression_inputs("rhybrid", rw_start = 2005)
  local_mock_models(inputs$eppasm$fp)

  fit_args <- list(optfit = TRUE, opthess = FALSE, B0 = 200, opt_maxit = 20)

  set.seed(10)
  expected <- with_time_limit(suppressWarnings(
    do.call(fitmod, c(list(inputs$eppasm_obj, eppmod = "rhybrid", rw_start = 2005),
                      fit_args))))
  set.seed(10)
  actual <- with_time_limit(suppressWarnings(
    do.call(fitmod_lf, c(list(inputs$lf_params, inputs$eppd, inputs$epp_t0,
                              eppmod = "rhybrid", rw_start = 2005),
                         fit_args))))

  # ll and ll_lf differ by floating point error (~1e-12), which the optimiser
  # amplifies over its iterations by an amount that varies by platform (e.g.
  # ancbias differed by ~3e-8 on CI), so use a looser tolerance than the default.
  expect_equal(actual$par, expected$par, tolerance = 1e-6)
  expect_equal(actual$value, expected$value, tolerance = 1e-6)
  expect_equal(actual$counts, expected$counts)
  expect_equal(standardise_names(actual$param), standardise_names(expected$param),
               tolerance = 1e-6)
})

# IMIS is a function of the likelihood evaluated at random draws, so with the
# same seed and the same likelihood values the two fits should be identical.
run_fitmod_pair <- function(inputs, eppmod, fit_args) {
  model_args <- c(list(eppmod = eppmod), regression_models[[eppmod]])
  eppasm_args <- c(list(inputs$eppasm_obj), model_args, fit_args)
  lf_args <- c(list(inputs$lf_params, inputs$eppd, inputs$epp_t0), model_args, fit_args)

  # IMIS prints progress for every stage
  utils::capture.output({
    set.seed(10)
    expected <- with_time_limit(suppressWarnings(do.call(fitmod, eppasm_args)))
    set.seed(10)
    actual <- with_time_limit(suppressWarnings(do.call(fitmod_lf, lf_args)))
  })
  list(expected = expected, actual = actual)
}

# Number of IMIS stages to compare. With these small sample sizes the rspline
# posterior is almost degenerate (effective sample size ~1 and a proposal
# covariance with condition number ~1e20), so floating point differences between
# ll and ll_lf (~1e-12) get amplified at later stages, e.g. at stage 6 theta
# differs by ~1e-8 and the likelihood by ~1e-4. Stages 1-5 are identical.
imis_number_k <- c(rhybrid = 6, rspline = 5, logrw = 6, rlogistic = 6)

for (eppmod in names(regression_models)) {
  local({
    eppmod <- eppmod

    test_that(paste0("fitmod_lf with IMIS gives the same fit as fitmod: ", eppmod), {
      inputs <- do.call(build_regression_inputs,
                        c(list(eppmod), regression_models[[eppmod]]))
      local_mock_models(inputs$eppasm$fp)

      fits <- run_fitmod_pair(
        inputs, eppmod,
        list(B0 = 500, B = 100, B.re = 200,
             number_k = imis_number_k[[eppmod]], opt_iter = 0))

      expect_s3_class(fits$actual, "specfit")
      expect_s3_class(fits$expected, "specfit")
      # Last column of stat is the time taken
      expect_equal(fits$actual$stat[, 1:6], fits$expected$stat[, 1:6])
      expect_equal(fits$actual$center, fits$expected$center)
      expect_equal(fits$actual$resample, fits$expected$resample)
    })
  })
}

test_that("fitmod_lf with IMIS and optimisation step gives a similar fit to fitmod", {
  # The optimisation step uses Nelder-Mead, which is chaotic: ll and ll_lf agree
  # to ~1e-15 relative but not bit for bit, and that is enough for the later
  # stages (and so the resample) to drift apart. Stage 1 happens before the
  # optimisation step so is still exact; after it only check the quantities
  # that should be stable: the log marginal likelihood at each stage and the
  # mode found by the optimiser.
  inputs <- build_regression_inputs("rhybrid", rw_start = 2005)
  local_mock_models(inputs$eppasm$fp)

  fits <- run_fitmod_pair(
    inputs, "rhybrid",
    list(B0 = 500, B = 100, B.re = 200, number_k = 6, opt_iter = 1))

  expect_equal(fits$actual$stat[1, 1:6], fits$expected$stat[1, 1:6])
  expect_equal(fits$actual$stat[, 1], fits$expected$stat[, 1], tolerance = 1e-2)
  expect_equal(fits$actual$center[[1]], fits$expected$center[[1]], tolerance = 1e-3)
})
