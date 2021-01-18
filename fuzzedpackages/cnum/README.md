
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Chinese Numerals in **R** 中文數字處理

<!-- badges: start -->

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/cnum)](https://cran.r-project.org/package=cnum)
[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

This R package provides useful functions to work with Chinese numerals
in R, such as conversion between Chinese numerals and Arabic numerals as
well as detection and extraction of Chinese numerals in character
objects and string.

*cnum* supports:

  - Conversion of any numbers of absolute value up to 10<sup>18</sup>
  - Conversion of decimal numbers
  - Traditional and Simplified Chinese
  - Financial numeral characters
  - Four scale naming systems including the casual and SI standard used
    in and outside mainland China

*cnum* 是協助處理中文數字的R套件，提供轉換、識別及抽取中文數字的函數。

本套件支援：

  - 任何數值不大於10<sup>18</sup>（一百京）之正負數轉換
  - 小數點轉換
  - 正體及簡體中文
  - 大小寫數字
  - 四種分別通用於中國大陸內外的日常及國際單位制進位命名法

*cnum* 是协助处理中文数字的R包，提供转换、识别及抽取中文数字的函数。

本包支援：

  - 任何数值不大於10<sup>18</sup>（一百京）之正负数转换
  - 小数点转换
  - 繁体及简体中文
  - 大小写数字
  - 四种分别通用於中国内地內外的日常及国际單位制进位命名法

## Installing

``` r
install.packages("cnum")

# To install the development version:
#install.packages("devtools")
devtools::install_github("elgarteo/cnum")
```

## Using

To convert from Arabic numerals to Chinese numerals:

``` r
num2c(721)
#> [1] "七百二十一"
num2c(-6)
#> [1] "負六"
num2c(3.14)
#> [1] "三點一四"
num2c(721, literal = TRUE) # instead of "七百二十一"
#> [1] "七二一"
num2c(1.45e12, financial = TRUE) # 大寫數字
#> [1] "壹兆肆仟伍佰億"
num2c(6.85e12, lang = "sc", mode = "casualPRC")
#> [1] "六万亿八千五百亿"
num2c(1.5e9, mode = "SIprefixPRC", single = TRUE) # instead of "一吉五百兆"
#> [1] "一點五吉"
```

To convert from Chinese numerals to Arabic numerals:

``` r
c2num("七百二十一")
#> [1] 721
c2num("負六")
#> [1] -6
c2num("三點一四")
#> [1] 3.14
c2num("七二一", literal = TRUE)
#> [1] 721
c2num("壹兆肆仟伍佰億", financial = TRUE)
#> [1] 1.45e+12
c2num("六万亿八千五百亿", lang = "sc", mode = "casualPRC")
#> [1] 6.85e+12
c2num("一點五吉", mode = "SIprefix")
#> [1] 1.5e+09
```

To detect Chinese numerals in character objects or string:

``` r
is_cnum("七百二十一")
#> [1] TRUE
is_cnum(c("非數字", "七百二十一"))
#> [1] FALSE  TRUE
is_cnum("七千三") # passes the casual test
#> [1] TRUE
is_cnum("七千三", strict = TRUE) # fails the strict test
#> [1] FALSE
is_cnum("壹兆肆仟伍佰億", financial = TRUE)
#> [1] TRUE
is_cnum("六万亿八千五百亿", lang = "sc", mode = "casualPRC")
#> [1] TRUE

has_cnum("加價一百八十元")
#> [1] TRUE
has_cnum(c("非數字", "加價一百八十元"))
#>         非數字 加價一百八十元 
#>          FALSE           TRUE
```

To extract Chinese numerals from string:

``` r
extract_cnum("分別加價一百八十元和五十六元")
#> [[1]]
#> [1] "一百八十" "五十六"
extract_cnum("十四亿人", lang = "sc")
#> [[1]]
#> [1] "十四亿"
extract_cnum("六万亿元", lang = "sc", mode = "casualPRC")
#> [[1]]
#> [1] "六万亿"
extract_cnum(c("加價一百八十元", "加價五十六元"))
#> [[1]]
#> [1] "一百八十"
#> 
#> [[2]]
#> [1] "五十六"
```

The default language is Traditional Chinese. This can be changed by
setting:

``` r
options(cnum.lang = "sc")
```

## To-Do

  - Add support for fractions and percentage
  - Add support for SI naming in decimal places (e.g. 1<sup>-6</sup> as
    一微)
