# Numerical equivalence between the native leapfrog engine (simmod_lf, driven
# from leapfrog's own data structure) and the original eppasm package, which is
# our independent reference. We compare on the direct-incidence workflow because
# there is no transmission-rate fitting to align first, so any difference is
# purely down to the model engines and their inputs.
#
# NOTE: the two engines are fed from *different* PJNZ readers
# (leapfrog::process_pjnz here vs eppasm's own reader in prepare_directincid), so
# perfect agreement also depends on those readers extracting identical inputs.
# See the ART-era characterisation test below.

test_that("native leapfrog reproduces eppasm incidence (direct incidence)", {
  skip_if_not_installed("eppasm")

  pjnz <- system_file("extdata/testpjnz", "Botswana2018.PJNZ")

  params <- process_directincid_lf(pjnz)
  sim <- simmod_lf(params)
  n <- params$sim_years

  fp <- eppasm::prepare_directincid(pjnz)
  mod <- eppasm::simmod(fp, VERSION = "C")

  # Direct incidence is an input, so 15-49 incidence must match to numerical
  # precision. This is the strongest invariant of the comparison.
  expect_equal(as.numeric(sim$incid_15to49)[seq_len(n)],
               as.numeric(eppasm::incid(mod))[seq_len(n)],
               tolerance = 1e-8)
})

test_that("native leapfrog reproduces eppasm pre-ART prevalence (direct incidence)", {
  skip_if_not_installed("eppasm")

  pjnz <- system_file("extdata/testpjnz", "Botswana2018.PJNZ")

  params <- process_directincid_lf(pjnz)
  sim <- simmod_lf(params)

  fp <- eppasm::prepare_directincid(pjnz)
  mod <- eppasm::simmod(fp, VERSION = "C")

  years <- params$projection_start_year + seq_len(params$sim_years) - 1L
  preart <- which(years <= 1996)

  # Before ART scale-up the two engines only differ in HIV natural history,
  # which is shared, so 15-49 prevalence should match closely.
  expect_equal(as.numeric(sim$prev_15to49)[preart],
               as.numeric(eppasm::prev(mod))[preart],
               tolerance = 1e-4)
})

test_that("native vs eppasm ART-era prevalence gap is within known bound", {
  # CHARACTERISATION TEST (not a correctness assertion). From ~1997 onward the
  # native leapfrog and eppasm prevalence trajectories diverge, peaking around
  # 1.5% absolute during ART scale-up while incidence stays identical. This
  # pins the current gap so a regression (or a fix) is noticed. The most likely
  # cause is a difference in ART inputs extracted by the two PJNZ readers rather
  # than the engines themselves; investigate before tightening this bound.
  skip_if_not_installed("eppasm")

  pjnz <- system_file("extdata/testpjnz", "Botswana2018.PJNZ")

  params <- process_directincid_lf(pjnz)
  sim <- simmod_lf(params)
  n <- params$sim_years

  fp <- eppasm::prepare_directincid(pjnz)
  mod <- eppasm::simmod(fp, VERSION = "C")

  max_prev_diff <- max(abs(sim$prev_15to49[seq_len(n)] -
                             as.numeric(eppasm::prev(mod))[seq_len(n)]))

  expect_lt(max_prev_diff, 0.02)
})
