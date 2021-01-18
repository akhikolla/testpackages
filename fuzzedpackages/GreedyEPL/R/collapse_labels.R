CollapseLabels <- function(decision)
{
  as.numeric(p__CollapseLabelsCpp(as.numeric(decision) - 1)$decision) + 1
}