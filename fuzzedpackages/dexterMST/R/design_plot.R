


#' Plot the routing design of mst tests
#' 
#' 
#' @param db dexterMST project database connection
#' @param predicate logical predicate to select data (tests, booklets,responses) to include in the design plot
#' @param by_booklet plot and color the paths in a test per booklet
#' @param ... further arguments to \code{\link[igraph]{plot.igraph}}
#' 
#' @details
#' You can use this function to plot routing designs for tests before or after
#' they are administered. There are some slight differences.
#' 
#' If you have entered response data already, the thickness of the line will indicate the numbers
#' of respondents that took the respective paths through the test. Paths not taken will not be drawn. You can use the
#' predicate (see examples) to include or exclude items, tests and respondents.
#' 
#' If you have not entered response data, all lines will have equal thickness. Variables you can use in the predicate
#' are limited to test_id and booklet_id in this case.   
#' 
#' 
#' 
#' 
#' @examples
#' \dontrun{
#' # plot test designs for all tests in the project
#' design_plot(db)
#' 
#' # plot design for a test with id 'math'
#' design_plot(db, test_id == 'math')
#' 
#' # plot design for test math with item 'circumference' turned off
#' # (this plot will only work if you have response data)
#' design_plot(db, test_id == 'math' & item_id != 'circumference')
#' 
#' }
#' 
design_plot = function(db, predicate = NULL, by_booklet=FALSE,...)
{
  user.args = list(...)
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env() 
  
  # to~do make work for a fit object
  
  # linear ignored for now
  
  
  if(is.null(qtpredicate))
  {
    routing = dbReadTable(db, 'booklet_design') %>%
      mutate(bid = paste(dense_rank(.data$test_id), dense_rank(.data$booklet_id)))
    
    admin =  dbGetQuery(db, 'SELECT test_id, booklet_id, COUNT(*) AS n FROM Administrations GROUP BY test_id, booklet_id;')
    
    if(NROW(admin)==0)
    {
      n_rsp = distinct(routing, .data$bid)
      n_rsp$n=0L
    } else
    {
      n_rsp = admin %>%
        mutate(bid = paste(dense_rank(.data$test_id), dense_rank(.data$booklet_id))) %>%
        select(-.data$test_id,-.data$booklet_id)
    }
    
  } else
  {
    mst_inputs = get_mst_data(db, qtpredicate, env)
    # to~do: niet helemaal correct - > module_id is geen identifier meer omdat voor sommige leerlingen met en somiige
    # zonder bepaalde items kan izjn
    #temp fix, bettr done in data selection
    if(!'module_id' %in% names(mst_inputs$routing))
    {
      md = dbGetQuery(db,'select test_id,booklet_id, module_id, module_nbr from booklet_design')
      mst_inputs$routing = inner_join(mst_inputs$routing,md,by=c('test_id','booklet_id','module_nbr'))
      
    }
    routing = mst_inputs$routing
    
    n_rsp = mst_inputs$suf_stats$plt %>%
      distinct(.data$bid, .data$booklet_score, .keep_all=TRUE) %>%
      group_by(.data$bid) %>%
      summarise(n=sum(.data$N)) %>%
      ungroup()
  }
  
  
  am = routing %>%
    group_by(.data$test_id, .data$bid) %>%
    arrange(.data$module_nbr) %>%
    mutate(mfrom = .data$module_id, mto = lead(.data$module_id),
          label = paste(.data$module_exit_score_min, .data$module_exit_score_max, sep=':')) %>%
    slice(-n()) %>%
    ungroup()
  
  
  # to~do: can show paths not taken, should show desiogn if no data yet
  if(by_booklet)
  {
    am_summed = am %>%
      inner_join(n_rsp,by='bid') %>%
      mutate(width = .data$n/max(.data$n,1) * 10) %>%
      mutate(width = pmax(1L,.data$width))
  } else
  {
    am_summed = am %>%
      inner_join(n_rsp,by='bid') %>%
      group_by(.data$test_id, .data$mfrom, .data$mto, .data$label) %>%
      summarise(width = sum(.data$n)) %>%
      ungroup() %>%
      mutate(width = .data$width/max(.data$width,1) * 10) %>%
      mutate(width = pmax(1L,.data$width))
  }


  
  op = par(mar=c(0.1,0.1,0.1,0.1))
  on.exit({par(op)}, add=TRUE)
  
  

  lapply(split(am_summed, am_summed$test_id), function(amsb)
  {
    g = graph.data.frame(amsb[,c('mfrom','mto','label','width')]) %>% 
      simplify(remove.multiple=FALSE, remove.loops=TRUE)
    
    ## Decompose the graph, individual layouts
    comp = decompose.graph(g)
    roots = sapply(lapply(comp, topological.sort), head, n=1)
    coords = mapply(FUN = layout.reingold.tilford, comp,
                    root = roots, SIMPLIFY=FALSE)
    
    ## Put the nodes side by side, roots on the top
    width = sapply(coords, function(x) { r <- range(x[, 1]); r[2] - r[1] })
    gap = 0.5
    shift = c(0, cumsum(width[-length(width)] + gap))
    ncoords = mapply(FUN=function(mat, shift) {
      mat[,1] = mat[,1] - min(mat[,1]) + shift
      mat[,2] = mat[,2] - max(mat[,2])
      mat
    }, coords, shift, SIMPLIFY=FALSE)
    
    ## Put together the coordinates for the original graph,
    ## based on the names of the vertices
    lay = matrix(0, ncol=2, nrow = vcount(g))
    for (i in seq_along(comp)) {
      lay[match(V(comp[[i]])$name, V(g)$name),] <- ncoords[[i]]
    }
    
    ## Plot everything
    
    #if(!is.null(mst_inputs$module_design_history))
    if(FALSE)#to~do: design history not yet implemented in new one, for v 1.0
    {
      mdl = mst_inputs$module_design_history %>%
        group_by(.data$miid, .data$module_id) %>%
        summarise(vlabel = paste0(.data$item_id[.data$rsp_incl == 0], collapse = ',')) %>%
        ungroup() %>%
        mutate(vlabel = paste0(.data$module_id, ifelse(.data$vlabel=='','','\n-'), .data$vlabel), miid = as.character(.data$miid))
      
      vlab = as_tibble(vertex_attr(g)) %>%
        mutate(indx = row_number()) %>%
        inner_join(mdl, by = c(name='miid')) %>%
        arrange(.data$indx) %>%
        pull(.data$vlabel)
      
      do.call(plot, modifyList(list(x=g, layout=lay, vertex.shape='rectangle',
                                    vertex.size=45,vertex.size2=20,vertex.label=vlab),
                               user.args))
    } else
    {
      #decide curves
      curves = as.data.frame(igraph::as_edgelist(g), stringsAsFactors=FALSE) %>% 
        mutate(rn=row_number()) %>%
        mutate(direction=if_else(.data$V1>.data$V2,1,-1), m1=pmin(.data$V1,.data$V2), m2=pmax(.data$V1,.data$V2)) %>% 
        group_by(.data$m1,.data$m2) %>% 
        mutate(curvature = .data$direction * if.else(n()==1,0,seq(-.5,.5,length.out=n()))) %>%
        arrange(.data$rn) %>%
        pull(.data$curvature)

      clr = NULL
      if(by_booklet)
      {
        clr = c("#FB8072", "#80B1D3", "#FDB462", "#8DD3C7", "#BEBADA", "#B3DE69", 
                "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", 
                "#F781BF", "#999999", "#BC80BD", "#CCEBC5", "#A50F15", "#08306B", 
                "#00441B", "#54278F", "#fc4E2A", "#525252", "#66C2A4")
        clr = rep_len(clr,n_distinct(amsb$bid))[dense_rank(amsb$bid)]
      }
      

      do.call(plot, modifyList(list(x=g, layout=lay, vertex.shape='rectangle',
                                    vertex.size=45,vertex.size2=20, edge.curved=curves, edge.color=clr),
                               user.args))
    }
  })

  invisible(NULL)
}

