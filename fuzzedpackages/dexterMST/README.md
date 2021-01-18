<!-- README.md is generated from README.Rmd. Please edit that file -->

DexterMST
=========

DexterMST is an R package acting as a companion to dexter and adding
facilities to manage and analyze data from multistage tests (MST). It
includes functions for importing and managing test data, assessing and
improving the quality of data through basic test and item analysis, and
fitting an IRT model, all adapted to the peculiarities of MST designs.
DexterMST typically works with project database files saved on disk.

Installation
------------

``` r
install.packages('dexterMST')
```

If you encounter a bug, please post a minimal reproducible example on
[github](https://github.com/jessekps/dexter/issues). We post news and
examples on a [blog](http://dexterities.netlify.com), it’s also the
place for general questions.

Example
-------

Here is an example for a simple two-stage test.

``` r
library(dexterMST)
library(dplyr)
# start a project
db = create_mst_project(":memory:")

items = data.frame(item_id=sprintf("item%02i",1:70), item_score=1, delta=sort(runif(70,-1,1)))

design = data.frame(item_id=sprintf("item%02i",1:70),
                    module_id=rep(c('M4','M2','M5','M1','M6','M3', 'M7'),each=10))

routing_rules = routing_rules = mst_rules(
 `124` = M1[0:5] --+ M2[0:10] --+ M4, 
 `125` = M1[0:5] --+ M2[11:15] --+ M5,
 `136` = M1[6:10] --+ M3[6:15] --+ M6,
 `137` = M1[6:10] --+ M3[16:20] --+ M7)


scoring_rules = data.frame(
  item_id = rep(items$item_id,2), 
  item_score= rep(0:1,each=nrow(items)),
  response= rep(0:1,each=nrow(items))) # dummy respons
  

db = create_mst_project(":memory:")
add_scoring_rules_mst(db, scoring_rules)

create_mst_test(db,
                test_design = design,
                routing_rules = routing_rules,
                test_id = 'sim_test',
                routing = "all")
```

We can now plot the design

``` r
# plot test designs for all tests in the project
design_plot(db)
```

We now simulate data:

``` r
theta = rnorm(3000)

dat = sim_mst(items, theta, design, routing_rules,'all')
dat$test_id='sim_test'
dat$response=dat$item_score

add_response_data_mst(db, dat)
```

``` r
# IRT, extended nominal response model
f = fit_enorm_mst(db)

head(f)
```

| item\_id |  item\_score|        beta|   SE\_beta|
|:---------|------------:|-----------:|----------:|
| item01   |            1|  -1.0967010|  0.0629430|
| item02   |            1|  -0.9396378|  0.0626093|
| item03   |            1|  -0.9362441|  0.0626050|
| item04   |            1|  -0.9226755|  0.0625888|
| item05   |            1|  -0.7974781|  0.0625298|
| item06   |            1|  -0.8549515|  0.0625365|

``` r
# ability estimates per person
rsp_data = get_responses_mst(db)
abl = ability(rsp_data, parms = f)
head(abl)
```

| booklet\_id | person\_id |  booklet\_score|       theta|
|:------------|:-----------|---------------:|-----------:|
| 136         | 1          |              14|   0.1155725|
| 136         | 10         |              20|   0.9510407|
| 136         | 100        |              15|   0.2505641|
| 124         | 1000       |              11|  -1.0122974|
| 125         | 1001       |              19|   0.3546790|
| 124         | 1002       |              15|  -0.4465699|

``` r
# ability estimates without item Item01
abl2 = ability(rsp_data, parms = f, item_id != "item01")

# plausible values
pv = plausible_values(rsp_data, parms = f, nPV = 5)
head(pv)
```

| booklet\_id | person\_id |  booklet\_score|         PV1|         PV2|         PV3|         PV4|         PV5|
|:------------|:-----------|---------------:|-----------:|-----------:|-----------:|-----------:|-----------:|
| 136         | 1          |              14|  -0.3310002|   0.1372317|  -0.1772800|  -0.5183684|  -0.0290442|
| 136         | 10         |              20|   0.8863248|   0.8113941|   0.6248103|   0.8908880|   1.0080691|
| 136         | 100        |              15|   0.8551317|   0.3975271|   0.4908577|  -0.1864416|  -0.0545868|
| 136         | 1003       |              16|   0.6045560|   0.1657788|   0.5029436|   0.3651413|   0.5977085|
| 136         | 1006       |              24|   1.5130258|   1.4255989|   0.7971183|   1.1151426|   1.3334579|
| 136         | 1008       |              15|   0.1326590|  -0.2617173|   0.0628930|   0.6675144|   0.2771222|

Contributing
------------

Contributions are welcome but please check with us first about what you
would like to contribute.
