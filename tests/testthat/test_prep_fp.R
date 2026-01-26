test_that("can prep leapfrog inputs for rhybrid model", {
  pjnz <- system_file("extdata/testpjnz", "Mozambique_Maputo_Cidade2018.PJNZ")
  inputs <- leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE)
  eppd <- prep_epp_data(pjnz)
  epp_t0 <- read_epp_t0(pjnz)
  prep <- prep_fp_fitmod_lf(inputs, eppd[["Maputo Cidade"]],
                            epp_t0["Maputo Cidade"],
                            eppmod = "rhybrid")

  # TODO: What to test here?
  expect_equal(prep$fp$proj_years, 52)
  expect_equal(prep$fp$sim_years, 48)
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
  eppd <- prep_epp_data(pjnz)
  epp_t0 <- read_epp_t0(pjnz)

  prep <- prep_fp_fitmod_lf(inputs, eppd[["Maputo Cidade"]],
                            epp_t0["Maputo Cidade"],
                            eppmod = "rspline")

  # TODO: What to test here?
  expect_equal(prep$fp$proj_years, 52)
  expect_equal(prep$fp$sim_years, 52)
  expect_equal(prep$fp$eppmod, "rspline")
  expect_equal(prep$fp$proj_steps[1] , 1970.5)
  expect_equal(length(prep$fp$proj_steps), 511)
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
  expect_equal(dim(prep$fp$rvec_spldes), c(511, 7))
  expect_equal(prep$fp$rtpenord, 2L)
  expect_true(is.null(prep$fp$iota))
})

test_that("can prep leapfrog inputs for rlogistic model", {
  pjnz <- system_file("extdata/testpjnz", "Mozambique_Maputo_Cidade2018.PJNZ")
  inputs <- leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE)
  eppd <- prep_epp_data(pjnz)
  epp_t0 <- read_epp_t0(pjnz)

  prep <- prep_fp_fitmod_lf(inputs, eppd[["Maputo Cidade"]],
                            epp_t0["Maputo Cidade"],
                            eppmod = "rlogistic")

  # TODO: What to test here?
  expect_equal(prep$fp$proj_years, 52)
  expect_equal(prep$fp$sim_years, 52)
  expect_equal(prep$fp$eppmod, "rlogistic")
  expect_equal(prep$fp$proj_steps[1] , 1970.5)
  expect_equal(length(prep$fp$proj_steps), 511)
  # TODO: We need all the below for leapfrog?
  expect_equal(prep$fp$ancrt, "both")
  expect_true(is.null(prep$fp$ancrtsite_beta))
  expect_true(prep$fp$ancsitedata)
  expect_equal(prep$fp$incidmod, "eppspectrum")

  # TODO: Test something about the returned items
  expect_equal(names(prep$likdat), c("hhs_dat", "ancsite_dat", "ancrtcens_dat"))

  # Test fit params
  expect_equal(prep$fp$ts_epidemic_start, 1975.5)
})

test_that("can create rvec and iota for rhybrid model", {
  pjnz <- system_file("extdata/testpjnz", "Mozambique_Maputo_Cidade2018.PJNZ")
  inputs <- leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE)
  eppd <- prep_epp_data(pjnz)
  epp_t0 <- read_epp_t0(pjnz)
  prep <- prep_fp_fitmod_lf(inputs, eppd[["Maputo Cidade"]],
                            epp_t0["Maputo Cidade"],
                            eppmod = "rhybrid")

  params <- fnCreateParam_lf(theta_rhybrid, prep$fp)

  expect_equal(names(params),
               c("rvec", "iota", "ancbias", "v_infl", "log_frr_adjust",
                 "frr_cd4", "frr_art", "ancrtcens_vinfl", "ancrtsite_beta"))
  expect_length(params$rvec, 472) # some test of values to within some tolerance?
  expect_length(params$frr_cd4, 2912)
  expect_length(params$frr_art, 8736)
})
