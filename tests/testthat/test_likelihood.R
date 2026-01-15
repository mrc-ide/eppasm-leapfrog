test_that("can run likelihood with leapfrog structured data", {
  pjnz <- system_file("extdata", "testpjnz", "Botswana2018.PJNZ")
  inputs <- leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE)
  eppd <- prep_epp_data(pjnz)
  epp_t0 <- read_epp_t0(pjnz)
  prep <- prep_fp_fitmod_lf(inputs, eppd$Urban, epp_t0["Urban"],
                            eppmod = "rhybrid")

  l <- likelihood_lf(theta_rhybrid, prep$fp, prep$likdat)

  # TODO: better tests here
  expect_true(!is.null(l))
  expect_length(l, 1)
})
