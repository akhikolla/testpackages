# coxrt - Cox Regression for Right Truncated Data

The  `coxrt` package fits Cox regression based on retrospectively ascertained time-to-event
data. The method uses Inverse-Probability-Weighting (IPW) estimating equations with stabilized weights as described
in Vakulenko-Lagun, Mandel, and Betensky, 
*Inverse-Probability-Weighting Methods for Cox Regression with Right Truncated Data* (under revision for *Biometrics*, 2019).

The `coxrt` package provides two functions: `coxph.RT` that assumes positivity and
`coxph.RT.a0` that allows for adjustment of estimation under violation of positivity. 

The  `coxrt` package can be installed by
```{r}
devtools::install_github("Bella2001/coxrt")
```
The examples of when and how to use `coxrt` package are [here](https://bella2001.github.io/coxrt/).
 
 
 
