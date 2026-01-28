# Age-specific prevalence among pregnant women

Age-specific prevalence among pregnant women

## Usage

``` r
agepregprev_lf(
  mod,
  fp,
  aidx = 3:9 * 5 - fp$ss$AGE_START + 1L,
  yidx = 1:fp$proj_years,
  agspan = 5,
  expand = FALSE
)
```

## Arguments

- mod:

  todo

- fp:

  Leapfrog fixed parameters

- expand:

  whether to expand aidx, yidx, sidx, and agspan
