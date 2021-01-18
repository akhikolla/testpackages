test_arithmetic <- function(){
wide_df = data.frame(
  subj = c("Tom", "Mary", "Jack"),
  t0 = c( 50, 42, 80),
  t1 = c( 48, 42, 75),
  t2 = c( 46, 44, 72),
  t3 = c( 42, 42, 73)
)

code = "
  subject = subj
  time = 0
  bw = t0
  push!()

  time = 1
  bw = t1 
  push!()

  time = 2
  bw = t2
  push!()

  time = 3
  bw = t3
"

long_df = datasailr::sail(wide_df , code = code, fullData=FALSE)
tom_long_df = long_df[long_df$subject=="Tom", ]

RUnit::checkEqualsNumeric( 4 * nrow(wide_df), nrow(long_df))
RUnit::checkEqualsNumeric( tom_long_df$bw , c(50, 48, 46, 42))

}
