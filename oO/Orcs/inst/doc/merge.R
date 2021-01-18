## ----pkg, echo=FALSE, message=FALSE-------------------------------------------
library(Orcs)

## ----mergeList, eval=TRUE-----------------------------------------------------
## sample data
set.seed(10)
ls_df <- list(data.frame(a = 1:10, b = 1:10),
              data.frame(a = 5:14, c = 11:20),
              data.frame(a = sample(20, 10), d = runif(10)))

## merge data frames in one go
merge(ls_df, by = "a", all = TRUE)

