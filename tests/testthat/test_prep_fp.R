test_that("can prep leapfrog inputs for rhybrid model", {
  pjnz <- system_file("extdata/testpjnz", "Mozambique_Maputo_Cidade2018.PJNZ")
  inputs <- leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE)
  inputs$proj_years <- 49
  eppd <- prep_epp_data(pjnz)
  prep <- prep_fp_fitmod_lf(inputs, eppd[["Maputo Cidade"]], eppmod = "rhybrid")

  # TODO: What to test here?
  expect_equal(prep$fp$eppmod, "rhybrid")
  expect_equal(prep$fp$proj_steps[1] , 1970.5)
  expect_equal(length(prep$fp$proj_steps), 471)
  # TODO: We need all the below for leapfrog?
  expect_equal(prep$fp$ancrt, "both")
  expect_true(is.null(prep$fp$ancrtsite_beta))
  expect_true(prep$fp$ancsitedata)
  expect_equal(prep$fp$incidmod, "eppspectrum")

  # TODO: Test something about the returned items
  expect_equal(names(prep$likdat), c("hhs_dat", "ancsite_dat", "ancrtcens_dat"))

  # Test fit params
  expect_equal(prep$fp$ts_epidemic_start, 1975.5)
  expect_equal(prep$fp$rt$eppmod, "rhybrid")
  expect_equal(prep$fp$rt$dt, 0.1)
  expect_equal(length(prep$fp$rt$rw_transition), 145)
  expect_equal(prep$fp$rt$rw_start, 2003)
  expect_equal(prep$fp$rt$proj_steps, prep$fp$proj_steps)
})

test_that("can prep leapfrog inputs for rspline model", {
  pjnz <- system_file("extdata/testpjnz", "Mozambique_Maputo_Cidade2018.PJNZ")
  inputs <- leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE)
  inputs$proj_years <- 49
  eppd <- prep_epp_data(pjnz)
  prep <- prep_fp_fitmod_lf(inputs, eppd[["Maputo Cidade"]], eppmod = "rspline")

  # TODO: What to test here?
  expect_equal(prep$fp$eppmod, "rspline")
  expect_equal(prep$fp$proj_steps[1] , 1970.5)
  expect_equal(length(prep$fp$proj_steps), 481)
  # TODO: We need all the below for leapfrog?
  expect_equal(prep$fp$ancrt, "both")
  expect_true(is.null(prep$fp$ancrtsite_beta))
  expect_true(prep$fp$ancsitedata)
  expect_equal(prep$fp$incidmod, "eppspectrum")

  # TODO: Test something about the returned items
  expect_equal(names(prep$likdat), c("hhs_dat", "ancsite_dat", "ancrtcens_dat"))

  # Test fit params
  expect_equal(prep$fp$ts_epidemic_start, 1975.5)
  expect_equal(prep$fp$num_knots, 7)
  expect_equal(dim(prep$fp$rvec_spldes), c(481, 7))
  expect_equal(prep$fp$rtpenord, 2L)
  expect_true(is.null(prep$fp$iota))
})

test_that("can prep leapfrog inputs for rlogistic model", {
  pjnz <- system_file("extdata/testpjnz", "Mozambique_Maputo_Cidade2018.PJNZ")
  inputs <- leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE)
  inputs$proj_years <- 49
  eppd <- prep_epp_data(pjnz)
  prep <- prep_fp_fitmod_lf(inputs, eppd[["Maputo Cidade"]], eppmod = "rlogistic")

  # TODO: What to test here?
  expect_equal(prep$fp$eppmod, "rlogistic")
  expect_equal(prep$fp$proj_steps[1] , 1970.5)
  expect_equal(length(prep$fp$proj_steps), 481)
  # TODO: We need all the below for leapfrog?
  expect_equal(prep$fp$ancrt, "both")
  expect_true(is.null(prep$fp$ancrtsite_beta))
  expect_true(prep$fp$ancsitedata)
  expect_equal(prep$fp$incidmod, "eppspectrum")

  # TODO: Test something about the returned items
  expect_equal(names(prep$likdat), c("hhs_dat", "ancsite_dat", "ancrtcens_dat"))

  # Test fit params
  expect_equal(prep$fp$ts_epidemic_start, 1975)
})

bw_theta <- c(-0.407503322169364, -2.76794181367538, -1.26018073624346, 1995.96447776502,
              -0.00307437171215574, 0.0114118307148102, 0.00760958379603691, 0.02,
              2.24103194827232, -0.0792123921862689, -5.01917961803606, 0.359444135205712,
              -6.10051517060137)

test_that("can create rvec and iota for rhybrid model", {
  pjnz <- system_file("extdata/testpjnz", "Mozambique_Maputo_Cidade2018.PJNZ")
  inputs <- leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE)
  inputs$proj_years <- 49
  eppd <- prep_epp_data(pjnz)
  prep <- prep_fp_fitmod_lf(inputs, eppd[["Maputo Cidade"]], eppmod = "rhybrid")

  params <- fnCreateParam_lf(bw_theta, prep$fp)

  expect_equal(names(params),
               c("rvec", "iota", "ancbias", "v_infl", "log_frr_adjust",
                 "frr_cd4", "frr_art", "ancrtcens_vinfl", "ancrtsite_beta"))
  expect_length(params$rvec, 472)
  expect_equal(params$frr_cd4, numeric(0))
  expect_equal(params$frr_art, numeric(0))
})
