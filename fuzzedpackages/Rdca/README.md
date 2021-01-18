
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Rdca

`Rdca` is a Decline Curve Analysis (DCA) package for oil and gas
reservoirs. It generates a table of rate, cumulative, nominal decline
rate, and derivative of loss-ratio over time in a data frame format. It
also provides an optimization tool to fit a DCA model on production
data. The package currently supports Arps ‘exponential’, ‘harmonic’,
‘hyperbolic’, and ‘modified\_hyperbolic’ models.

## Installation

You can install the released version of Rdca from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("Rdca")
```

### `Arps Exponential Examples`

``` r
library(Rdca)
library(magrittr)
library(ggplot2)
library(ggpubr)

dcl_param_exp <- decline_param(input_unit = "Field", output_unit = "Field", fluid = "oil", 
                               model = "exponential", qi = 1000, Di = 0.0015, b = 0, q_abnd = NULL)

dcl_param_exp
#> $input_unit
#> [1] "Field"
#> 
#> $output_unit
#> [1] "Field"
#> 
#> $fluid
#> [1] "oil"
#> 
#> $qi
#> [1] 1000
#> 
#> $Di
#> [1] 0.0015
#> 
#> $b
#> [1] 0
#> 
#> attr(,"class")
#> [1] "exponential" "decline"

decline_time_exp <- decline_time(c(1:7300), unit = "day")   

str(decline_time_exp)
#> List of 3
#>  $ t             : int [1:7300] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ unit          : chr "day"
#>  $ reference_date: Date[1:1], format: "2020-05-13"
#>  - attr(*, "class")= chr [1:2] "day" "time"

decline_predict_exp <- decline_predict(dcl_param_exp, decline_time_exp)

head(decline_predict_exp, 10)
#>          Date Time_(day) q_(bbl/day)   Q_(bbl) D_(1/day) Beta
#> 1  2020-05-13          1    998.5011  999.2504    0.0015    0
#> 2  2020-05-14          2    997.0045 1997.0030    0.0015    0
#> 3  2020-05-15          3    995.5101 2993.2601    0.0015    0
#> 4  2020-05-16          4    994.0180 3988.0240    0.0015    0
#> 5  2020-05-17          5    992.5281 4981.2968    0.0015    0
#> 6  2020-05-18          6    991.0404 5973.0808    0.0015    0
#> 7  2020-05-19          7    989.5549 6963.3783    0.0015    0
#> 8  2020-05-20          8    988.0717 7952.1914    0.0015    0
#> 9  2020-05-21          9    986.5907 8939.5225    0.0015    0
#> 10 2020-05-22         10    985.1119 9925.3736    0.0015    0
```

### `Arps Harmonic Examples`

``` r
library(Rdca)
library(ggplot2)
library(ggpubr)

dcl_param_harm <- decline_param(input_unit = "SI", output_unit = "SI", fluid = "oil", 
                               model = "harmonic", qi = 1000, Di = 0.075, b = 1, q_abnd = 50)

dcl_param_harm
#> $input_unit
#> [1] "SI"
#> 
#> $output_unit
#> [1] "SI"
#> 
#> $fluid
#> [1] "oil"
#> 
#> $qi
#> [1] 1000
#> 
#> $Di
#> [1] 0.075
#> 
#> $b
#> [1] 1
#> 
#> $q_abnd
#> [1] 50
#> 
#> attr(,"class")
#> [1] "harmonic" "decline"

decline_time_harm <- decline_time(c(1:360), unit = "month")   

str(decline_time_harm)
#> List of 3
#>  $ t             : int [1:360] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ unit          : chr "month"
#>  $ reference_date: Date[1:1], format: "2020-05-13"
#>  - attr(*, "class")= chr [1:2] "month" "time"

decline_predict_harm <- decline_predict(dcl_param_harm, decline_time_harm)

head(decline_predict_harm, 10)
#>          Date Time_(month) q_(m3/month)    Q_(m3) D_(1/month) Beta
#> 1  2020-05-13            1     930.2326  964.2755  0.06976744    1
#> 2  2020-06-12            2     869.5652 1863.4926  0.06521739    1
#> 3  2020-07-12            3     816.3265 2705.8779  0.06122449    1
#> 4  2020-08-12            4     769.2308 3498.1902  0.05769231    1
#> 5  2020-09-11            5     727.2727 4246.0497  0.05454545    1
#> 6  2020-10-12            6     689.6552 4954.1808  0.05172414    1
#> 7  2020-11-11            7     655.7377 5626.5921  0.04918033    1
#> 8  2020-12-11            8     625.0000 6266.7151  0.04687500    1
#> 9  2021-01-11            9     597.0149 6877.5089  0.04477612    1
#> 10 2021-02-10           10     571.4286 7461.5438  0.04285714    1
#>    time_abnd_(months) EUR_(m3)
#> 1            253.3333  39943.1
#> 2            253.3333  39943.1
#> 3            253.3333  39943.1
#> 4            253.3333  39943.1
#> 5            253.3333  39943.1
#> 6            253.3333  39943.1
#> 7            253.3333  39943.1
#> 8            253.3333  39943.1
#> 9            253.3333  39943.1
#> 10           253.3333  39943.1
```

### `Arps Hyperbolic Examples`

``` r
library(Rdca)
library(ggplot2)
library(ggpubr)

dcl_param_hyp <- decline_param(input_unit = "Field", output_unit = "Field", fluid = "gas", 
                               model = "hyperbolic", qi = 100000, Di = 0.0055, b = 0.85, 
                               q_abnd = 2000)

dcl_param_hyp
#> $input_unit
#> [1] "Field"
#> 
#> $output_unit
#> [1] "Field"
#> 
#> $fluid
#> [1] "gas"
#> 
#> $qi
#> [1] 1e+05
#> 
#> $Di
#> [1] 0.0055
#> 
#> $b
#> [1] 0.85
#> 
#> $q_abnd
#> [1] 2000
#> 
#> attr(,"class")
#> [1] "hyperbolic" "decline"

decline_time_hyp <- decline_time(seq(as.Date("2000/1/1"), as.Date("2030/12/31"), "days"), 
                                 unit = "date")   

str(decline_time_hyp)
#> List of 3
#>  $ t             : num [1:11323] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ unit          : chr "date"
#>  $ reference_date: Date[1:1], format: "2000-01-01"
#>  - attr(*, "class")= chr [1:2] "day" "time"

decline_predict_hyp <- decline_predict(dcl_param_hyp, decline_time_hyp)

head(decline_predict_hyp, 10)
#>          Date Time_(day) q_(MSCF/day) Q_(MMSCF)   D_(1/day) Beta
#> 1  2000-01-01          1     99452.78  99.72593 0.005474407 0.85
#> 2  2000-01-02          2     98911.08 198.90741 0.005449051 0.85
#> 3  2000-01-03          3     98374.81 297.54991 0.005423929 0.85
#> 4  2000-01-04          4     97843.90 395.65882 0.005399038 0.85
#> 5  2000-01-05          5     97318.26 493.23947 0.005374374 0.85
#> 6  2000-01-06          6     96797.83 590.29708 0.005349934 0.85
#> 7  2000-01-07          7     96282.51 686.83683 0.005325716 0.85
#> 8  2000-01-08          8     95772.26 782.86379 0.005301716 0.85
#> 9  2000-01-09          9     95266.98 878.38300 0.005277931 0.85
#> 10 2000-01-10         10     94766.61 973.39938 0.005254359 0.85
#>    time_abnd_(days) EUR_(MMSCF)
#> 1          5733.712    53805.81
#> 2          5733.712    53805.81
#> 3          5733.712    53805.81
#> 4          5733.712    53805.81
#> 5          5733.712    53805.81
#> 6          5733.712    53805.81
#> 7          5733.712    53805.81
#> 8          5733.712    53805.81
#> 9          5733.712    53805.81
#> 10         5733.712    53805.81
```
