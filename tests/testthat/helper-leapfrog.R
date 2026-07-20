theta_rhybrid <- c(-0.407503322169364, -2.76794181367538, -1.26018073624346, 1995.96447776502,
                   -0.00307437171215574, 0.0114118307148102, 0.00760958379603691, 0.02,
                   2.24103194827232, -0.0792123921862689, -5.01917961803606, 0.359444135205712,
                   -6.10051517060137)

# Build leapfrog params ready for simmod_lf() for a direct-incidence run. This is
# the native equivalent of eppasm::prepare_directincid(): process the PJNZ into
# leapfrog's own data structure and add the few book-keeping fields that
# simmod_lf() reads. No transmission fitting is involved.
process_directincid_lf <- function(pjnz) {
  params <- suppressMessages(
    leapfrog::process_pjnz(pjnz, use_coarse_age_groups = TRUE))
  ss <- leapfrog::get_leapfrog_ss(LEAPFROG_MODEL_CONFIG)
  params$proj_years <- as.integer(dim(params$Sx)[length(dim(params$Sx))])
  params$sim_years <- params$proj_years
  params$ss$hiv_steps_per_year <- ss$HIV_STEPS_PER_YEAR
  params
}

expect_mod_equivalent <- function(mod, mod_expected) {
  expect_equal(class(mod), class(mod_expected))
  expect_equal(dim(mod), dim(mod_expected))
  attribute_names <- setdiff(names(attributes(mod_expected)), c("dim", "class"))
  known_missing_attributes <- c("natdeaths_noart", "natdeaths_art", "popadjust",
                                "pregprevlag", "entrantprev")
  known_different_attributes <- c("rvec_ts")
  attribuites_to_check <- setdiff(
    attribute_names,
    c(known_missing_attributes, known_different_attributes))

  for (attr_name in attribuites_to_check) {
    expected <- attr(mod_expected, attr_name)
    received <- attr(mod, attr_name)
    if (is.array("expected")) {
      expect_equal(dim(received), dim(expected), info = attr_name)
    } else {
      expect_equal(length(received), length(expected), info = attr_name)
    }
  }
}
