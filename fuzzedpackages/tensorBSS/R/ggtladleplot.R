ggtladleplot <- function(x, crit = "gn", type="l", ylab = crit,  xlab = "component", main = deparse(substitute(x)), ...)
{
  comp <- NULL
  DF <- do.call(rbind, sapply(x$ResMode,List2dataframe, simplify=FALSE))
  crit <- match.arg(crit, c("gn", "fn", "phin", "lambda"))
  type <- ifelse(type == "l", 1, 0)
  ggplot(DF, aes(x = comp, y = eval(parse(text = crit)))) +
    geom_point() +
    geom_line(alpha = type) +
    facet_grid(. ~ mode) +
    labs(x = xlab, y = ylab, title = main) +
    ggtitle(main) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}