# Setup the r-hybrid model

Setup the r-hybrid model

## Usage

``` r
prepare_rhybrid_lf(
  params,
  ts_epidemic_start = 1975.5,
  rw_start = params$rw_start,
  rw_trans = params$rw_trans,
  rw_dk = params$rw_dk
)
```

## Arguments

- rw_start:

  time when random walk starts

- rw_trans:

  number of years to transition from logistic differences to RW
  differences. If NULL, defaults to 5 steps

- fp:

  Leapfrog parameters object

- tsEpidemicStart:

  time step at which epidemic is seeded
