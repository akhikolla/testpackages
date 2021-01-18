
#' @importFrom  stats kmeans


kmeans.framed.clust <-
  function(X, K, frame.width, first.frame, last.frame)
  {
    X <- sort(X)

    result <- kmeans(X[1:frame.width], K)

    value <- result$tot.withinss

    ID <- 0

    size <- result$size[unique(result$cluster)]

    Border <- cumsum(size) - 1


    first.frame <- first.frame + 1

    last.frame <- last.frame + 1

    if ((first.frame + 1) <= last.frame)
    {
      for (i in (first.frame + 1):last.frame)
      {
        points = c(X[i:(i + frame.width - 1)])

        result_new <- kmeans(points, K)

        if (result_new$tot.withinss < value)
        {
          value <- result_new$tot.withinss

          ID <- i - 1

          size <- result$size[unique(result$cluster)]

          Border <- cumsum(size) + ID - 1

          result <- result_new
        }
      }
    }

    df <-
      list(
        "cluster" = result$cluster,
        "centers" = result$centers,
        "withinss" = result$withinss,
        "size" = result$size,
        "totss" = result$totss,
        "tot.withinss" = value,
        "betweenss" = result$betweenss,
        "ID" = ID,
        "Border" = Border
      )

    return(df)
  }
