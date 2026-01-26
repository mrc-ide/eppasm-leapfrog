test_that("prior works as expected", {
  pjnz <- system_file("extdata", "testpjnz", "Botswana2018.PJNZ")
  inputs <- leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE)
  # When running with transmission input i.e. not direct incidence EPPASM
  # will run the model for the years for which there is data and then project
  # forward to fill in up to PROJ_YEARS
  eppd <- prep_epp_data(pjnz)
  epp_t0 <- read_epp_t0(pjnz)
  prep <- prep_fp_fitmod_lf(inputs, eppd$Urban, epp_t0["Urban"],
                            eppmod = "rhybrid")

  fp <- modifyList(prep$fp, fnCreateParam_lf(theta_rhybrid, prep$fp))

  lval <- prior_lf(theta_rhybrid, fp)

  # TODO: better tests, what can we actually check here
  expect_true(is.numeric(lval))
  expect_length(lval, 1)
})

test_that("sample_prior_lf works as expected", {
  pjnz <- system_file("extdata/testpjnz", "Botswana2018.PJNZ")
  inputs <- leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE)
  # When running with transmission input i.e. not direct incidence EPPASM
  # will run the model for the years for which there is data and then project
  # forward to fill in up to PROJ_YEARS
  eppd <- prep_epp_data(pjnz)
  epp_t0 <- read_epp_t0(pjnz)
  prep <- prep_fp_fitmod_lf(inputs, eppd$Urban, epp_t0["Urban"],
                            eppmod = "rhybrid")

  fp <- modifyList(prep$fp, fnCreateParam_lf(theta_rhybrid, prep$fp))

  mat <- sample_prior_lf(1e5, fp)

  # TODO: better tests, what can we actually check here
  expect_equal(dim(mat), c(100000, 10))
})
