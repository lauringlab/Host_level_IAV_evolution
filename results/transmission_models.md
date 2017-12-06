Transmission models
================
JT McCrone
4/7/2017

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

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_logical(),
    ##   X1 = col_integer(),
    ##   chr = col_character(),
    ##   pos = col_integer(),
    ##   SPECID1 = col_character(),
    ##   SPECID2 = col_character(),
    ##   X1_1 = col_integer(),
    ##   HOUSE_ID = col_integer(),
    ##   ENROLLID1 = col_character(),
    ##   ENROLLID2 = col_character(),
    ##   onset1 = col_date(format = ""),
    ##   onset2 = col_date(format = ""),
    ##   transmission = col_date(format = ""),
    ##   pair_id = col_double(),
    ##   collect1 = col_date(format = ""),
    ##   collect2 = col_date(format = ""),
    ##   mutation = col_character(),
    ##   ref = col_character(),
    ##   var = col_character(),
    ##   season = col_character(),
    ##   pcr_result = col_character()
    ##   # ... with 2 more columns
    ## )

    ## See spec(...) for full column specifications.

Presence/Absence model
======================

Let \(A_1\) and \(A_2\) be a alleles in the donor. Then there are three possible outcomes of transmission. Either \(A_1\) is transmitted, \(A_2\) is transmitted, or both \(A_1\) and \(A_2\) are transmitted. The probability of only one allele being transmitted (let's call it \(A_i\)) given a bottleneck size of \(N_b\) is \[
P(A_i) = f_i^{N_b}
\]

Where \(p_i\) is the frequency of the allele in the donor. In otherwords this is simply the probability of only drawing \(A_i\) in \(N_b\) draws.

The probability of both alleles being transmitted is given by

\[
P(A_1,A_2) = 1- \big(f_1^{N_b}+f_2^{N_b}\big)
\]

where \(f_1\) and \(f_2\) are the frequencies of the alleles respectively. This is simply the probability of not picking only \(A_1\) or only \(A_2\) in \(N_b\) draws.

We can then define the probability of observing the data at each polymorphic site \(j\) as \(P_j\) where \(P_j=P(A_i)\) if only one allele is transmitted and \(P_j=P(A_1,A_2)\) if two alleles are transmitted.

The likelihood of a bottleneck size \(N_b\) is then given by

\[
L(N_b) = \prod^j P_j
\]

or the probability of observing the data at each polymorphic site. Thus the log likelihood is given by

\[
\text{log}(L(N_b)) = \sum^j\text{log}(P_j)
\]

### Fitting the Presence/Absence model

In this fit we take the minority frequency to the correct and set the major frequency to 1-minority

    ## # A tibble: 1 x 4
    ##   lambda lower_95 upper_95  mean_Nb
    ##    <dbl>    <dbl>    <dbl>    <dbl>
    ## 1   1.12     0.51     1.99 1.662411

### Log likelihood plot

<img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

    ## Loading required package: cowplot

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

### Fits by pair

For fun <img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

### Simulation

<img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

all transmitted \>10%
=====================

What if all variants above 10% were transmitted? What would the bottleneck size be?

<img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

    ## # A tibble: 1 x 4
    ##   lambda lower_95 upper_95  mean_Nb
    ##    <dbl>    <dbl>    <dbl>    <dbl>
    ## 1   5.85     3.96     8.44 5.866896

<img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

<img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-13-1.png" style="display: block; margin: auto;" /> This sort of makes sense. If there was a larger bottleneck then we would expect the low - frequency variants would have higher probabilities of transmission as well.

No frequency cut off
====================

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_logical(),
    ##   X1 = col_integer(),
    ##   chr = col_character(),
    ##   pos = col_integer(),
    ##   SPECID1 = col_character(),
    ##   SPECID2 = col_character(),
    ##   X1_1 = col_integer(),
    ##   HOUSE_ID = col_integer(),
    ##   ENROLLID1 = col_character(),
    ##   ENROLLID2 = col_character(),
    ##   onset1 = col_date(format = ""),
    ##   onset2 = col_date(format = ""),
    ##   transmission = col_date(format = ""),
    ##   pair_id = col_double(),
    ##   collect1 = col_date(format = ""),
    ##   collect2 = col_date(format = ""),
    ##   mutation = col_character(),
    ##   ref = col_character(),
    ##   var = col_character(),
    ##   season = col_character(),
    ##   pcr_result = col_character()
    ##   # ... with 2 more columns
    ## )

    ## See spec(...) for full column specifications.

    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson
    ## Nonfinite probabilities given by the zerotruncated Poisson - using a Poisson

    ## # A tibble: 1 x 4
    ##   lambda lower_95 upper_95 mean_Nb
    ##    <dbl>    <dbl>    <dbl>   <dbl>
    ## 1    118      112      124     118

The warning reduces to the 3 sites wher ethe total frequency is very close to 99%. Our assumption that major allele = 1-minor should hold just fine.

    ## # A tibble: 3 x 4
    ## # Groups:   pair_id, chr [3]
    ##   pair_id   chr   pos  sum.freq
    ##     <dbl> <chr> <int>     <dbl>
    ## 1   108.0   PB2   150 0.9759139
    ## 2   108.5   PB2   150 0.9776132
    ## 3   181.0    PA   513 0.9884250

    ## Joining, by = c("chr", "pos", "pair_id")

    ## # A tibble: 6 x 9
    ##   SPECID1 SPECID2 HOUSE_ID   chr   pos   ref   var       freq1       freq2
    ##     <chr>   <chr>    <int> <chr> <int> <chr> <chr>       <dbl>       <dbl>
    ## 1  HS1450  MH8309     5289    PA   513     C     C 0.983606557 1.000000000
    ## 2  HS1450  MH8309     5289    PA   513     C     T 0.004818486 0.000000000
    ## 3  HS1394  MH7843     5002   PB2   150     G     A 0.974780391 0.974650206
    ## 4  HS1394  MH7843     5002   PB2   150     G     G 0.001133466 0.002962963
    ## 5  MH7843  HS1394     5002   PB2   150     G     A 0.974650206 0.974780391
    ## 6  MH7843  HS1394     5002   PB2   150     G     G 0.002962963 0.001133466

### Log likelihood plot

<img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

No cut simulation
-----------------

<img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

No cut no infer
===============

Let's get rid of sites with infered minor alleles. <img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-20-1.png" style="display: block; margin: auto;" /><img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-20-2.png" style="display: block; margin: auto;" />

<img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-22-1.png" style="display: block; margin: auto;" />

    ## # A tibble: 1 x 4
    ##   lambda lower_95 upper_95  mean_Nb
    ##    <dbl>    <dbl>    <dbl>    <dbl>
    ## 1   1.67     0.91     2.71 2.057276

<img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-25-1.png" style="display: block; margin: auto;" />

Beta binomial model
===================

The Beta binomial model is explained in detail in Leonard . It is similar to the presence/absence model in that transmission is modeled as a simple sampling process; however, it loosens the assumption that the frequencies in the recipient are constant overtime. Instead frequencies of transmitted variants are allowed to change between transmission and sampling according the a beta distribution. The distribution is not dependent on the amount of time that passes between transmission and sampling, and the frequency in the donor is assumed to be the same between sampling and transmission.

The equations below are very similar to those present by Leonard with two exceptions. (1) We fit a distribution to the bottleneck sizes in our cohort, and (2) because we know our sensitivity to detect rare variants based on the expected frequency of the iSNV and the titer of the sample we can include the possiblity that iSNV are transmitted but are missed due to less than perfect sensitivity.

\[
L(N_b)_i=\sum_{k=0}^{N_b}\text{p_beta}\Big( _{R,i}|k,N_b-k\Big)\text{p_bin}\Big(k|N_b,v_{D,i}\Big)
\]

and

I will start with the most conservative assumption. We will always round the titer and frequency down to the nearest standard and apply that accuracy. Also I'm assuming the accuracy is perfect in the donor.

So now the likelihood function of lost variants is given by

\[
L(N_b)_i=\sum_{k=0}^{N_b}\Big[\text{p_beta_cdf}\Big( v_{R,i}<T|k,N_b-k\Big)\text{p_bin}\Big(k|N_b,v_{D,i}\Big)+\sum_{f_i}^{[0.02,0.05,0.1)}\text{p_beta}\big(f_i<v_{R,i}<f_{i+1}\big|k,N_b-k\big)\text{p_bin}\Big(k|N_b,v_{D,i}\Big)\big(1-\text{sensitivity}|\text{titer}_R,f_i)\Big]
\] In other words what is the probability the (variant was not transmitted or transmitted but remains \<2%) or ( the variant was transmitted and is present within a given frequency range and we don't find it given the lower end of that frequency range and the titer of the sample.)

Fitting
-------

    ## # A tibble: 1 x 4
    ##   lambda lower_95 upper_95  mean_Nb
    ##    <dbl>    <dbl>    <dbl>    <dbl>
    ## 1   1.22     0.57     2.17 1.731062

### Loglikelihodd plot

<img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-27-1.png" style="display: block; margin: auto;" />

### Fitting by person

<img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-28-1.png" style="display: block; margin: auto;" />

Simulations
-----------

<img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-29-1.png" style="display: block; margin: auto;" /> \#\# lambda = 10

<img src="transmission_models_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-30-1.png" style="display: block; margin: auto;" />

|Model|AIC|
|:----|--:|
|Presence/Absence|76.69269|
|BetaBinomial|83.02161|

|Model|lambda|lower\_95|upper\_95|mean\_Nb|
|:----|-----:|--------:|--------:|-------:|
|Presence Absence|1.12|0.51|1.99|1.662411|
|BetaBinomial|1.22|0.57|2.17|1.731062|

Figure 3
========
