
##########################################
#' Exploratory test for Differential Item Functioning
#' 
#' Compares two parameter objects and produces a test for DIF
#' based on equality of relative item difficulties category locations
#'
#' @param db an dexterMST db handle
#' @param person_property name of a person property defined in your dexterMST project
#' @param predicate logical predicate to select data to include in the analysis
#' 
#' @examples
#' \dontrun{
#' 
#' dif = DIF_mst(db, person_property = 'test_mode')
#' print(dif)
#' plot(dif)
#' 
#' }
#' @references 
#' Bechger, T. M. and Maris, G (2015); A Statistical Test for Differential Item Pair Functioning. 
#' Psychometrika. Vol. 80, no. 2, 317-340.
#' 
DIF_mst = function(db, person_property, predicate=NULL) 
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env() 
  
  items = get_rsp_data(db, qtpredicate = qtpredicate, env = env, columns = c('item_id', person_property)) %>% 
    distinct()
  
  pp_vals = sort(unique(items[[person_property]]))
  
  if(length(pp_vals) != 2)
    stop('person property must have two unique values in your data')
  
  spl = split(pull(items, 'item_id'), items[[person_property]])
  common_items = intersect(spl[[1]], spl[[2]])
  
  removed_items = items %>%
    distinct(.data$item_id) %>% 
    anti_join(tibble(item_id = common_items), by = 'item_id')
  
  # to~do: common item scores
  if(nrow(removed_items) > 0)
  {
    if(length(common_items) == 0)
      stop('There are no common items between your two populations')
    message(paste(nrow(removed_items), 
                  'items have been removed from the analysis because they do not occur in both groups'))

  }
  
  if(is.null(qtpredicate))
  {
    models = lapply(pp_vals, function(pp)
      {
        qtp = as.call(list(as.name('=='), as.name(person_property),pp))
        fit_enorm_mst_(db, qtpredicate=qtp,env=env)
      })
  } else
  {
    models = lapply(pp_vals, function(pp)
    {
      qtp = as.call(list(as.name('=='), as.name(person_property),pp))
      qtp = as.call(list(as.name('&'), qtp, qtpredicate))
      fit_enorm_mst_(db, qtpredicate = qtp, env = env)
    })
  }
  
  new_items = coef(models[[1]]) %>%
    inner_join(tibble(item_id = common_items), by='item_id') %>%
    select(.data$item_id, .data$item_score) %>%
    arrange(.data$item_id, .data$item_score) 
  
  
  if(nrow(removed_items) > 0)
  {
    models = lapply(models, function(m)
    {
      items = coef(m) %>%
        arrange(.data$item_id, .data$item_score) %>%
        mutate(rn = row_number()) %>%
        inner_join(tibble(item_id = common_items), by='item_id') %>%
        select(.data$item_id, .data$item_score, .data$rn) %>%
        arrange(.data$rn) 
      
      m$est$beta = m$est$beta[items$rn, , drop=FALSE]
      m$est$acov.beta = m$est$acov.beta[items$rn, items$rn]
      m
    })
  }
  
  names(models) = pp_vals

  ## 4. Call overallDIF_ and PairDIF_
  DIF_stats = dexter.OverallDIF_(models[[1]]$est$beta, models[[2]]$est$beta, 
                           models[[1]]$est$acov.beta, models[[2]]$est$acov.beta)
  
  D = dexter.PairDIF_(models[[1]]$est$beta, models[[2]]$est$beta, 
               models[[1]]$est$acov.beta, models[[2]]$est$acov.beta)
  
  
  
  ## 5. Report D and DIF_stats and inputs
  ou = list(DIF_overall = DIF_stats, DIF_pair = D, 
            items = new_items, group_labels = names(models))
  
  class(ou) = append('DIF_stats_mst', class(ou))
  return(ou)
}


print.DIF_stats_mst <- function(x, ...)
{
  specify_decimal <- function(x, k) format(round(x, k), nsmall=k)
  tmp = specify_decimal(x$DIF_overall$p,3)
  if (tmp=="0.000") tmp="< 0.0006"
  cat(paste0("Test for DIF:"," Chi-square = ", as.character(round(x$DIF_overall$stat, digits=3)),
               ", df = ", 
               as.character(x$DIF_overall$df),
               ", p = ", tmp))  
}

#' plot method for DIF_mst
#' 
#' @param x object produced by DIF_mst
#' @param items character vector of item id's for a subset of the plot. Useful if you have many items. 
#' If NULL all items are plotted.
#' @param itemsX character vector of item id's for the X axis
#' @param itemsY character vector of item id's for the Y axis
#' @param ... further arguments to plot
#'    
plot.DIF_stats_mst = function(x, items = NULL, itemsX = items, itemsY = items, ...)
{

  if(is.null(itemsX)) itemsX = unique(x$items$item_id)
  if(is.null(itemsY)) itemsY = unique(x$items$item_id)
  
  if(length(setdiff(c(itemsX, itemsY), x$items$item_id)) > 0)
  {
    cat('items not found in DIF object:\n')
    print(setdiff(c(itemsX, itemsY), x$items))
    stop('some of the item_ids you specified are not present in the DIF object')
  }
  
  x$items = x$items %>%
    mutate(rn = row_number())
  
  itemsX = x$items %>%
    semi_join(tibble(item_id=itemsX), by='item_id') %>%
    arrange(.data$rn)
  
  itemsY = x$items %>%
    semi_join(tibble(item_id=itemsY), by='item_id') %>%
    arrange(.data$rn)
  
  DIF_pair = x$DIF_pair[itemsX$rn, itemsY$rn]

  
  if(nrow(distinct(x$items,.data$item_id)) == nrow(x$items))
  {
    yLabels = pull(itemsY, 'item_id')
    xLabels = pull(itemsX, 'item_id')
  } else
  {
    yLabels = paste(itemsY$item_id, itemsY$item_score)
    xLabels = paste(itemsX$item_id, itemsX$item_score)

  }

  min_ = min(x$DIF_pair) # keep color range from complete object
  max_ = max(x$DIF_pair)
  default.args = list(main = paste(x$group_labels[1],'vs.',x$group_labels[2]),
                      axes=FALSE, zlim=c(min_,max_),xlab='',ylab='')
  
  graphics::layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  
  tmp = rainbow(256)[1:128]
  ColorRamp=c(tmp, tmp[128:1])
  ColorLevels <- seq(min(x$DIF_pair), max(x$DIF_pair), length=length(ColorRamp))
  
  # Reverse Y axis
  yLabels <- rev(yLabels)
  DIF_pair <- DIF_pair[nrow(DIF_pair) : 1,]
  
  # Data Map
  oldpar = par(mar = c(6,8,2.5,2))
  on.exit({par(oldpar)}, add=TRUE)
  do.call(image,
          merge_arglists(list(...),
                         override = list(x = 1:length(yLabels), y = 1:length(xLabels), z=t( DIF_pair),
                                         col=ColorRamp),
                         default = default.args))
  
  
  axis(1, at=1:length(yLabels), labels=yLabels, las= 3, cex.axis=0.6, hadj=1,padj=1)
  axis(2, at=1:length(xLabels), labels=xLabels, las=1, cex.axis=0.6, hadj=1,padj=1)
  
  #Color Scale
  par(mar=c(6,3,2,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  graphics::layout(1)

}