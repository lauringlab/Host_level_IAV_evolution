Benchmarking accuracy
================
JT McCrone
April 3,2017

Here is the accuracy from our benchmarking experiments.

|  gc\_ul|   freq|  sensitivity|  specificity|   TP|   FP|        PPV|
|-------:|------:|------------:|------------:|----:|----:|----------:|
|   1e+05|  0.050|         1.00|       0.9999|   20|    1|  0.9523810|
|   1e+05|  0.020|         0.85|       0.9999|   17|    6|  0.7391304|
|   1e+05|  0.010|         0.95|       0.9995|   19|   20|  0.4871795|
|   1e+05|  0.005|         0.35|       0.9999|    7|    5|  0.5833333|
|   1e+04|  0.050|         0.95|       0.9999|   19|    4|  0.8260870|
|   1e+04|  0.020|         0.90|       0.9999|   18|    6|  0.7500000|
|   1e+04|  0.010|         0.80|       0.9998|   16|   10|  0.6153846|
|   1e+04|  0.005|         0.40|       0.9999|    8|    2|  0.8000000|
|   1e+03|  0.050|         0.80|       0.9999|   16|    1|  0.9411765|
|   1e+03|  0.020|         0.45|       0.9999|    9|    3|  0.7500000|
|   1e+03|  0.010|         0.20|       0.9997|    4|   12|  0.2500000|
|   1e+03|  0.005|         0.10|       0.9999|    2|    4|  0.3333333|

While this is correct it does not account for the frequency of the False positives which is usually lower than the TP until we reach the rare variants. *Note : We did not set any frequency thresholds in the benchmarking paper. Here we will apply frequency thresholds. When samples are sequenced twice we will require both frequencies to be greater than 2%. That will disqualify some of the FP in the table above. This is consistent with how we have processed the patient sample. Well not exactly but I'll have to rerun the variant caller for that.And then this will be correct*

    ## [1] "We removed 723 variants of 774 and 51 remain."

    ## [1] "We removed 825 variants of 909 and 84 remain."

Here I have applied the stringency mentioned above and have calculated the FP as the number of FP in the give sample with a frequency above the expected frequency of the TP. This lowers the sensitivity by basically give perfect specificity.

|  exp.freq|  freq|   TP|  sensitivity|   FP|  gc\_ul|        PPV|
|---------:|-----:|----:|------------:|----:|-------:|----------:|
|      0.05|  0.05|   17|         0.85|    0|   1e+05|  1.0000000|
|      0.02|  0.02|    3|         0.15|    0|   1e+05|  1.0000000|
|      0.05|  0.05|   17|         0.85|    0|   1e+04|  1.0000000|
|      0.02|  0.02|    3|         0.15|    0|   1e+04|  1.0000000|
|      0.05|  0.05|   14|         0.70|    0|   1e+03|  1.0000000|
|      0.02|  0.02|    3|         0.15|    0|   1e+03|  1.0000000|
|      0.01|  0.01|    1|         0.05|    2|   1e+03|  0.3333333|

So in almost every case
$$
PPV=\\frac{TP}{TP+FP}=1
$$
 If we say something is there then you better believe it is. Unless it's present at 1% in a 10<sup>3</sup> sample, then it's probably a false positive. So it's probably not worth rerunning the anlaysis at this point yet. There are some transmitted variants in this category but they are mostly inferred reciprocal variants which need to be handled differently. I will create a similar accuracy table for inferred reciprocal variants, and then I'll include that accuracy and the above accuracy in the models.
