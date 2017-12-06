Summary statistics
================
JT McCrone
10/30/2017

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

    ## Loading required package: doMC

    ## Loading required package: foreach

    ## 
    ## Attaching package: 'foreach'

    ## The following objects are masked from 'package:purrr':
    ## 
    ##     accumulate, when

    ## Loading required package: iterators

    ## Loading required package: parallel

First we'll read in the all\_meta file. Which now contains every sample we ever touched. We will filter it so that it contains 1 entry per person. This handels the cases in 2014-2015 where we have mulitple samples/ person. Although the multiple samples are important when we pick which sample to use in looking at transmission, here we just are looking at dates on onset so one sample/ person will do.

When there are 2 samples we pick the one we sequenced unless we sequenced both then we take the one that qualified for snv identification with the titer as a tie braker. If we didn't sequence either sample we pick the one to include randomly.

I should note we are taking one sample per person per infecting strain. There are cases where one person was sick twice in a season with H1N1 and H3N2. Those each counted in the meta\_one. The pcr\_result is used to identify transmission pairs later.

### Cohort stats

How many H3N2 and H1N1 did we sequence

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

### Table 1 stats

Take the meta and split by strain and season.

|pcr\_result|2010-2011|2011-2012|2012-2013|2013-2014|2014-2015|
|:----------|--------:|--------:|--------:|--------:|--------:|
|A/H1N1|26|1|3|47|0|
|A/H3N2|58|22|66|1|166|

### Transmission rules.

These apply to all cases where 2 individuals are sick within the same household within a week of eachtoher (difference in date of onset \<= 7 days).

In the event of multiple possible donors we assume the donor is the individual with symptom onset nearest to the recipeient.

The donor and recipeint are never allowed to have symptom onset on the same day unless they are the only cases in the house. In this case we will randomize the pair and estimate a bottleneck in both directions.

If there are two possible donors with the same date of onset then we throw out any possible pair to that recipient.

Also we require that pairs have an L1-norm below the 5% percentile of the non household pairs.

First we will just need one sample per person per pcr result. Date of onset is important for this analysis not date of sampling - that comes in later.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Duplicated column names deduplicated: 'X1' => 'X1_2' [2]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_integer(),
    ##   X1_2 = col_integer(),
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

Tranmission pairs are cases where 2 people are sick with the same strain from the same household with 7 days of eachother (inclusive)

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_logical(),
    ##   X1 = col_integer(),
    ##   HOUSE_ID = col_integer(),
    ##   pcr_result = col_character(),
    ##   season = col_character(),
    ##   ENROLLID1 = col_character(),
    ##   ENROLLID2 = col_character(),
    ##   onset1 = col_date(format = ""),
    ##   onset2 = col_date(format = ""),
    ##   transmission = col_date(format = ""),
    ##   gc_ul1 = col_double(),
    ##   gc_ul2 = col_double(),
    ##   pair_id = col_integer()
    ## )

    ## See spec(...) for full column specifications.

This yields 185 possible pairs. And 124 valid pairs

Summary table of transmission events
------------------------------------

However, we did not sequence every sample and every sample we sequenced did not yeild usable SNV data. We handeled 493 samples. We sequenced 367 samples or 74% of the samples. 255 samples (52%) had titers \> 1000 genomes/ul and were sequenced.

We had samples from 390 inidividuals. We sequenced samples from 298 or 76% of the individuals.

So the best we could hope for would be 87

When we further restrict our analysis and require both samples to be greater than 10<sup>3</sup> then we have 53 possible pairs. The 71 I quoted earlier included sequenced early samples and all 2014-2015 samples we handeled. I have remade the meta data csv for this analysis.

L1-norm
-------

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Duplicated column names deduplicated: 'X1' => 'X1_1' [51]

Now we will make all possible pair comparisions within a season and strain. Note - becuause we are using the meta\_one data frame we are not accounting for intrahost dynamics (in season 2014-2015). We have set up meta\_one to give us the highest quality iSNV data.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_integer(),
    ##   season = col_character(),
    ##   pcr_result = col_character(),
    ##   SPECID1 = col_character(),
    ##   SPECID2 = col_character(),
    ##   ENROLLID1 = col_character(),
    ##   ENROLLID2 = col_character(),
    ##   time_onset = col_integer(),
    ##   time_collect = col_integer(),
    ##   L1_norm = col_double(),
    ##   HOUSE_ID1 = col_integer(),
    ##   HOUSE_ID2 = col_integer(),
    ##   Household = col_logical(),
    ##   valid = col_logical(),
    ##   quality_distance = col_logical()
    ## )

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

|L1\_norm|threshold|valid\_pairs|
|-------:|--------:|-----------:|
|0.00000|0.00|0|
|15.98791|0.05|47|
|22.13775|0.10|47|
|32.09874|0.15|47|
|45.40643|0.20|48|
|66.70239|0.25|48|
|70.45920|0.30|48|
|73.13025|0.35|48|
|75.71780|0.40|48|
|77.82552|0.45|50|
|79.13835|0.50|51|
|80.76210|0.55|51|
|82.55170|0.60|51|
|84.48883|0.65|52|
|86.79044|0.70|52|
|89.03344|0.75|52|
|92.04885|0.80|52|
|96.47817|0.85|52|
|162.17490|0.90|52|
|293.94104|0.95|52|
|544.45895|1.00|52|

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-2.png" style="display: block; margin: auto;" />

    ## [1] 1592

    ## [1] 47

    ## [1] 15.98791

Overall Summary
---------------

|Class|All.Samples|Sequenced.samples|Titers1e3|SNV.sequenced|
|:----|----------:|----------------:|--------:|:------------|
|Households|240|191|135|133|
|Isolates|493|367|255|249|
|Individuals|390|298|205|200|
|Transmission pairs|124|87|53|52 (47)|
|Households with potential pairs|85|64|39|38|
|Longitudinal sampling|103|69|50|49|

|Individuals|2010-2011|2011-2012|2012-2013|2013-2014|2014-2015|
|----------:|--------:|--------:|--------:|--------:|--------:|
|2|13|2|9|7|23|
|3|5|2|3|3|11|
|4|NA|NA|1|2|4|
|\#\# Looking at|transmission|||||

Here are the functions used for plotting.

Valid pairs
===========

Transmission rules. These apply to all cases where 2 individuals are sick within the same household within a week of eachtoher (difference in date of onset \<= 7 days).

In the event of multiple possible donors we assume the donor is the individual with symptom onset neast to the recipeient

The donor and recipeint are never have symptoms on the same day unless they are the only cases in the house. In this case we will randomize the pair and estimate a bottleneck both ways.

We'll take a look at the pairs we are refering to and then get a list of all those that quailify.

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-15-1.png" style="display: block; margin: auto;" /><img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-15-2.png" style="display: block; margin: auto;" /> This is the in set figure for some talks.

Houses of interest

Here is a final summary table. These are the valid pairs

    ## Loading required package: tables

    ## Loading required package: Hmisc

    ## Loading required package: lattice

    ## Loading required package: survival

    ## Loading required package: Formula

    ## 
    ## Attaching package: 'Hmisc'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, src, summarize

    ## The following objects are masked from 'package:base':
    ## 
    ##     format.pval, round.POSIXt, trunc.POSIXt, units

    ##                                                                         
    ##            Group                                                        
    ##            All                 Sequenced              Titer$>$1e3       
    ##            A/H1N1 A/H3N2 total A/H1N1    A/H3N2 total A/H1N1      A/H3N2
    ##  2010-2011  5      19     24    2        11     13    2            3    
    ##  2011-2012  0       6      6    0         5      5    0            1    
    ##  2012-2013  0      18     18    0        11     11    0            2    
    ##  2013-2014 18       0     18   13         0     13    6            0    
    ##  2014-2015  0      58     58    0        45     45    0           39    
    ##  total     23     101    124   15        72     87    8           45    
    ##                                  
    ##                                  
    ##        SNV qualified             
    ##  total A/H1N1        A/H3N2 total
    ##   5    2              2      4   
    ##   1    0              1      1   
    ##   2    0              2      2   
    ##   6    6              0      6   
    ##  39    0             39     39   
    ##  53    8             44     52

SNV summary
===========

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:Hmisc':
    ## 
    ##     combine, src, summarize

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

### Titers

There are two samples here with NA DPI, I believe we don't have meta data on when the symptoms began. Also there are two samples with negative DPI. These are all removed in the DPI plots unless otherwise obvious.

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 13 rows containing non-finite values (stat_boxplot).

    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/all_titers-1.png" style="display: block; margin: auto;" /> There are warning messages for 13 samples we did not get usable numbers back from the qPCR. They have NA in the log\_copy\_num column but R gives them 0 in the gc\_ul. These are removed on the log scale.

### Frequency distribution

|class\_factor|mutations|
|:------------|--------:|
|Nonsynonymous|277|
|Synonymous|434|

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/frequency_distribution-1.png" style="display: block; margin: auto;" />

    ## # A tibble: 10 x 2
    ##      bin     ratio
    ##    <int>     <dbl>
    ##  1     1 0.6776860
    ##  2     2 1.0000000
    ##  3     3 0.7368421
    ##  4     4 0.9230769
    ##  5     5 0.3750000
    ##  6     6 0.5454545
    ##  7     7 1.0000000
    ##  8     8 0.2608696
    ##  9     9 0.2592593
    ## 10    10 0.5000000

0.6104079

### Distribution across the genome

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/isnv_genome-1.png" style="display: block; margin: auto;" />

### Distribution across Individuals - using every sample

These will need to be rare in mulitple individuals to show up here at \>1

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 26 rows containing missing values (geom_text).

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/mutations_per_sample-1.png" style="display: block; margin: auto;" /> These are the mutations found in multiple individuals <img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-19-1.png" style="display: block; margin: auto;" /><img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-19-2.png" style="display: block; margin: auto;" /> These counts don't sum to the total above, but that is because this is looking at individuals above is looking at all sequenced samples (multiple/person) sometimes. If I sum these counts over SPECID then the sums are equal.

### Diversity in samples

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/snv_sample-1.png" style="display: block; margin: auto;" />

### Distribution of snv in samples (Could be multiple samples/person)

    ##   0%  25%  50%  75% 100% 
    ##    0    1    2    3   57

    ## [1] 2.855422

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-22-1.png" style="display: block; margin: auto;" /><img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-22-2.png" style="display: block; margin: auto;" /><img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-22-3.png" style="display: block; margin: auto;" />

Figure 1
========

    ## Loading required package: cowplot

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 13 rows containing non-finite values (stat_boxplot).

    ## notch went outside hinges. Try setting notch=FALSE.

    ## notch went outside hinges. Try setting notch=FALSE.

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-23-1.png" style="display: block; margin: auto;" />

Separate Panels

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 13 rows containing non-finite values (stat_boxplot).

    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.

Figure 3
========

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-25-1.png" style="display: block; margin: auto;" />

Supplemental Figure 4
=====================

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-27-1.png" style="display: block; margin: auto;" />

Looking at the outliers
=======================

Let's take a look at the samples with many snv.

    ## Loading required package: ggjoy

    ## Loading required package: ggridges

    ## The ggjoy package has been deprecated. Please switch over to the
    ## ggridges package, which provides the same functionality. Porting
    ## guidelines can be found here:
    ## https://github.com/clauswilke/ggjoy/blob/master/README.md

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-29-1.png" style="display: block; margin: auto;" />

    ## Picking joint bandwidth of 0.00905

<img src="Summary_stats_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-29-2.png" style="display: block; margin: auto;" /> This looks like mixed infections. We find many mutations at similar frequencies.
