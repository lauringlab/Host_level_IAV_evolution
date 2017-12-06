mutational model
================
JT McCrone
4/18/2017

    ## Loading required package: knitr

    ## Loading required package: ggplot2

    ## Loading required package: magrittr

    ## Loading required package: tidyverse

    ## ── Attaching packages ────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ tibble  1.3.4     ✔ purrr   0.2.4
    ## ✔ tidyr   0.7.2     ✔ dplyr   0.7.4
    ## ✔ readr   1.1.1     ✔ stringr 1.2.0
    ## ✔ tibble  1.3.4     ✔ forcats 0.2.0

    ## ── Conflicts ───────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ tidyr::extract()   masks magrittr::extract()
    ## ✖ dplyr::filter()    masks stats::filter()
    ## ✖ dplyr::lag()       masks stats::lag()
    ## ✖ purrr::set_names() masks magrittr::set_names()

    ## Loading required package: HIVEr

    ## Loading required package: extrafont

    ## Registering fonts with R

    ## Loading required package: wesanderson

    ## Loading required package: grid

    ## Loading required package: cowplot

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Duplicated column names deduplicated: 'X1' => 'X1_1' [51]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Duplicated column names deduplicated: 'X1' => 'X1_1' [2]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_integer(),
    ##   X1_1 = col_integer(),
    ##   HOUSE_ID = col_integer(),
    ##   ENROLLID = col_character(),
    ##   SPECID = col_character(),
    ##   onset = col_date(format = ""),
    ##   collect = col_date(format = ""),
    ##   vaccination_status = col_integer(),
    ##   pcr_result = col_character(),
    ##   LAURING_ID = col_character(),
    ##   DPI = col_integer(),
    ##   season = col_character(),
    ##   log_copy_num = col_double(),
    ##   gc_ul = col_double(),
    ##   HIGHSD = col_character(),
    ##   sequenced = col_logical(),
    ##   home_collected = col_integer(),
    ##   snv_qualified = col_logical()
    ## )

<img src="mutational_model_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-3-1.png" style="display: block; margin: auto;" /> \# Estimating *N*<sub>*e*</sub> and *μ*<sub>*e*</sub> The goal here is to allow for the mutations to enter the system at a given rate so that we can go backwards in time when we need to. Here we are relying on the assumption that no polymorphisms are present at the start of infection

Here the probability function for 0 &lt; *f* &lt; 1 is given by

$$
g(f,t) = \\frac{2 \\mu N\_e}{f} e^{- \\frac{2N\_ef}{t}}
$$
 I am assumining anything present at less than 0.001 is present at 0.

6 hours
=======

|     mu|   Ne|         LL|
|------:|----:|----------:|
|  4e-06|   36|  -4686.507|

<img src="mutational_model_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

12 hours
========
