test_that("fit with imis works with leapfrog data", {
  # TODO: Maybe we want to turn this off on CI as it is quite slow
  # is there a way we can make it run faster
  pjnz <- system_file("extdata", "testpjnz", "Botswana2018.PJNZ")
  inputs <- leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE)
  # When running with transmission input i.e. not direct incidence EPPASM
  # will run the model for the years for which there is data and then project
  # forward to fill in up to PROJ_YEARS
  eppd <- prep_epp_data(pjnz)
  epp_t0 <- read_epp_t0(pjnz)

  fit <- fitmod_lf(inputs, eppd$Urban,
                   epp_t0["Urban"],
                   eppmod = "rhybrid",
                   rw_start = 2005,
                   B0=1e3,
                   B=1e2,
                   opt_iter = 1,
                   number_k=50)

  # TODO: What to test about this?
  expect_true(!is.null(fit))
  expect_s3_class(fit, "specfit")
})

test_that("fit with optfit works with leapfrog data", {
  pjnz <- system_file("extdata", "testpjnz", "Botswana2018.PJNZ")
  inputs <- leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE)
  eppd <- prep_epp_data(pjnz)
  epp_t0 <- read_epp_t0(pjnz)

  fit <- fitmod_lf(inputs, eppd$Urban,
                   epp_t0["Urban"],
                   eppmod = "rhybrid",
                   optfit = TRUE,
                   opthess = FALSE,
                   rw_start = 2005,
                   B0=1e3,
                   B=1e2,
                   opt_iter = 1,
                   number_k=50)

  # TODO: What to test about this?
  expect_true(!is.null(fit))
  expect_true(all(c("fp", "likdat", "par", "mod") %in% names(fit)))
})

# TODO: add a test with opthess = TRUE (do we have any appropriate data for this?)
