
ffactor = function (x, levels=NULL, as_int=FALSE) 
{
  if(is.null(levels))
  {
    fast_factor(x, as_int)
  } else
  {
    fast_factor_lev(x, levels, as_int)
  }
}

# can/likely to return unused levels
bid = function(respData, db)
{
  tlev = dbGetQuery(db,'SELECT test_id FROM Tests ORDER BY test_id;')$test_id
  if(length(tlev)==0)
    stop('database contains no data')
  
  tblev = dbGetQuery(db,'SELECT test_id, booklet_id FROM Booklets ORDER BY test_id, booklet_id;')
  
  tblev$bid = 1:NROW(tblev)
  lev = paste(tblev$test_id, tblev$booklet_id, sep='-')
  if(anyDuplicated(lev))
  {
    tchar = max(nchar(tblev$test_id))
    lev = sprintf(paste0('% -',tchar,'s-%s'),tblev$test_id, tblev$booklet_id)
  }
  
  class(tblev$bid) = 'factor'
  levels(tblev$bid) = lev
  
  if(length(tlev)==1)
  {
    out = ffactor(respData$booklet_id,levels=tblev$booklet_id)
  } else
  {
    out = bid_c(respData$test_id, respData$booklet_id, tblev$test_id, tblev$booklet_id)
  }
  class(out) = 'factor'
  levels(out) = lev
  list(bid=out, bid_tb = tblev)
}


#' Extract response data from a dexterMST database
#' 
#'
#' @param db a dexterMST project database connection
#' @param predicate an expression to select data on
#' @param columns the columns you wish to select, can include any column in the project
#' @return a data.frame of responses
#' 
get_responses_mst = function(db, predicate = NULL, 
                             columns=c('person_id', 'test_id', 'booklet_id', 'item_id', 'item_score'))
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env() 
  
  get_rsp_data(db, qtpredicate, env=env, columns = columns)
}


get_cte = function(db, columns)
{
  cte = "Responses"
  columns = setdiff(columns, dbListFields(db, 'Responses'))
  
  tbls = list(
    Scoring_rules = 'INNER JOIN Scoring_rules USING(item_id, response)',
    Persons = 'INNER JOIN Persons USING(person_id)',
    Items = 'INNER JOIN Items USING(item_id)',
    Administrations = 'INNER JOIN Administrations USING(person_id, test_id, booklet_id)',
    Tests = 'INNER JOIN Tests USING(test_id)',
    Booklets = 'INNER JOIN Booklets USING(test_id, booklet_id)',
    Modules = 'INNER JOIN Modules USING(test_id, module_id)',
    Module_design = 'INNER JOIN Module_design USING(test_id, module_id, item_id)',
    Booklet_design = 'INNER JOIN Booklet_design USING(test_id, booklet_id, module_id)')
  
  for(tb in names(tbls))
  {
    if(length(columns)==0)
      break
    flds = dbListFields(db, tb)
    added = intersect(columns, flds)
    if(length(added)>0)
    {
      columns = setdiff(columns, added)
      cte = c(cte, tbls[[tb]])
    }
  }

  paste0(cte,collapse=" ")
}


get_rsp_data = function(db, qtpredicate=NULL, env=NULL,
                        columns=c('person_id', 'test_id', 'booklet_id', 'item_id', 'item_score'),
                        filter_col = NULL)
{

  pred = qtpredicate_to_sql(qtpredicate, db, env)
  
  used_columns = union(columns,pred$db_vars)

  
  cte = get_cte(db, used_columns)
    

  respData = NULL
  if(pred$succes)
  {
    if(is.null(filter_col))
    {
      respData = try(dbGetQuery(db, 
                            paste("SELECT", 
                                  paste0(columns, collapse=','),
                                  "FROM",
                                  cte,
                                  pred$where,';')),silent=TRUE)
    } else
    {
      # filtercol necessitates a predicate
      respData = try(dbGetQuery(db, 
                            paste("SELECT", 
                                  paste0(columns, collapse=','),',',
                                  'CASE WHEN (',pred$sql,') THEN 1 ELSE 0 END AS ', filter_col,
                                  " FROM ",
                                  cte)), silent = TRUE)
    }
  }

  
  if(!pred$success || inherits(respData,'try-error'))
  {
    #print('sql failed')
    respData = dbGetQuery(db,
                            paste("SELECT",
                                  paste0(used_columns, collapse=','),
                                  " FROM ",
                                  cte),';') 
    if(is.null(filter_col))
    {
      respData =respData[eval_tidy(qtpredicate, data = respData, env=env), columns]
    } else 
    {
      respData[[filter_col]] = as.integer(eval_tidy(qtpredicate, data = respData, env=env))
    }
  }

  respData
}



mst_safe = function(db, qtpredicate,env){
  blacklist = 
    setdiff(
      Reduce(union, list(
        dbListFields(db, 'Scoring_rules'),
        dbListFields(db, 'Booklet_design'),
        dbListFields(db, 'Module_design'),
        dbListFields(db, 'Items'),
        dbListFields(db, 'Scoring_rules'))),
      c('test_id', 'booklet_id')) 
  

  pred = qtpredicate_to_sql(qtpredicate, db, env)
  
  length(intersect(tolower(pred$db_vars),tolower(blacklist))) == 0
}


get_mst_variables = function (db) 
{
  lapply(dbListTables(db), 
         function(tbl) 
         {
           res = dbSendQuery(db, paste("SELECT * FROM", tbl, 
                                            "WHERE 0=1;"))
           r = dbColumnInfo(res)
           dbClearResult(res)
           return(r)
         }) %>% 
    bind_rows() %>% 
    distinct() %>% 
    arrange(.data$name)
}

get_mst_data = function(db, qtpredicate=NULL, env=NULL)
{
  if(is.null(env)) env = caller_env() 
  
  if(is.null(qtpredicate) || mst_safe(db, qtpredicate,env))
  {
    safe_mst_data(db, qtpredicate, env)
  } else
  {
    unsafe_mst_data(db, qtpredicate, env)
  }
}



safe_mst_data = function(db, qtpredicate, env)
{
  # avoids all time consuming computations for mutilating routing rules, etc.
  
  columns = c('person_id', 'test_id', 'booklet_id', 'item_id', 'item_score')
  
  respData = get_rsp_data(db, qtpredicate=qtpredicate, env=env,
                          columns=columns,
                          filter_col = NULL)
  
  respData$person_id = ffactor(respData$person_id, 
                               levels=dbGetQuery(db, 'SELECT person_id FROM Persons ORDER BY person_id;')$person_id,
                               as_int=TRUE)
  respData$item_id = ffactor(respData$item_id, 
                               levels=dbGetQuery(db, 'SELECT item_id FROM Items ORDER BY item_id;')$item_id)
  
  bid_info = bid(respData, db)
  respData$bid = bid_info$bid
  respData$test_id = NULL
  respData$booklet_id = NULL
  
  if(!is_person_booklet_sorted(respData$bid, respData$person_id))
    respData = arrange(respData, .data$person_id, .data$bid)
  
  respData$booklet_score = mutate_booklet_score(respData$person_id, respData$bid, respData$item_score)
  
  max_score=dbGetQuery(db,"SELECT MAX(item_score) FROM scoring_rules;")[1,1]
  
  suf_stats = suf_stats_nrm(respData, max_score)
  
  bid_tb = semi_join(bid_info$bid_tb, suf_stats$plt, by='bid')
  
  routing = dbGetQuery(db, 'SELECT * FROM Tests INNER JOIN booklet_design USING(test_id) 
                              ORDER BY test_id, booklet_id, module_nbr;') %>%
    inner_join(bid_tb, by=c('test_id','booklet_id'))
  

  booklet_design = dbGetQuery(db, 'SELECT test_id, booklet_id, item_id, module_nbr
                              FROM booklet_design INNER JOIN module_design USING(test_id, module_id)') %>%
    inner_join(bid_tb, by=c('test_id','booklet_id')) %>%
    arrange(.data$bid,.data$module_nbr, .data$item_id)
  
  booklet_design$item_id = ffactor(booklet_design$item_id, levels(respData$item_id))
  

  list(x = respData, 
       booklet_design= booklet_design,
       routing = routing,
       suf_stats = suf_stats)
}



unsafe_mst_data = function(db, qtpredicate,  env)
{
 
  columns=c('person_id', 'test_id', 'booklet_id', 'module_nbr','item_id', 'item_score')
  
  respData = get_rsp_data(db, qtpredicate=qtpredicate, env=env,
                          columns=columns,
                          filter_col = 'rsp_incl')
  
  respData$person_id = ffactor(respData$person_id, 
                               levels=dbGetQuery(db, 'SELECT person_id FROM Persons ORDER BY person_id;')$person_id,
                               as_int=TRUE)
  
  respData$item_id = ffactor(respData$item_id, 
                             levels=dbGetQuery(db, 'SELECT item_id FROM Items ORDER BY item_id;')$item_id)
  
  bid_info = bid(respData, db)
  respData$bid = bid_info$bid
  respData$test_id = NULL
  respData$booklet_id = NULL
  
  if(!is_person_booklet_sorted(respData$bid, respData$person_id))
    respData = arrange(respData, .data$person_id, .data$bid)
  
  respData$booklet_score = 0L
  
  
  nmod = dbGetQuery(db, 'SELECT test_id, booklet_id, COUNT(*) AS nmod FROM Booklet_design GROUP BY test_id, booklet_id;') %>%
    inner_join(bid_info$bid_tb, by=c('test_id','booklet_id')) %>%
    arrange(as.integer(.data$bid)) %>%
    pull(.data$nmod)
  
  nmod = as.integer(c(-1L,nmod)) # 1 indexing
  
  dsg = make_booklets_unsafe(respData$person_id, respData$bid,
                             respData$module_nbr, respData$item_id, respData$item_score,
                             respData$booklet_score, respData$rsp_incl,
                             nmod)
  
  
  names(respData)[names(respData) == 'booklet_id'] = 'bid'

  respData = respData[respData$rsp_incl==1L, names(respData) !='rsp_incl']
  
  max_score = dbGetQuery(db,"SELECT MAX(item_score) FROM scoring_rules;")[1,1]
  
  suf_stats = suf_stats_nrm(respData, max_score)
  
  item_ms = suf_stats$ssIS %>%
    group_by(.data$item_id) %>%
    summarise(item_maxscore = max(.data$item_score))
  
  bdes = dbGetQuery(db,"SELECT test_id, booklet_id, module_nbr, item_id
                            FROM booklet_design 
                            INNER JOIN Module_design USING(test_id, module_id);")
  
  bdes$item_id = ffactor(bdes$item_id,levels=levels(suf_stats$ssIS$item_id))
  
  class(dsg$booklets_items$item_id) = 'factor'
  levels(dsg$booklets_items$item_id)= levels(suf_stats$ssIS$item_id)
  
  bdes = bdes  %>% 
    inner_join(item_ms, by='item_id') %>%
    inner_join(bid_info$bid_tb, by=c('test_id','booklet_id')) %>%
    select(old_bid='bid', .data$module_nbr, .data$item_id, .data$item_maxscore) %>%
    mutate(old_bid=as.integer(.data$old_bid)) %>%
    inner_join(dsg$booklets_items, by=c('old_bid','item_id'))
  
  mod_max = bdes %>%
    group_by(.data$bid, .data$module_nbr) %>%
    summarise(mod_maxscore = sum(.data$item_maxscore))
  
  routing = dbGetQuery(db,'SELECT test_id, booklet_id, routing, module_exit_score_min, module_exit_score_max, module_nbr
                      FROM Tests INNER JOIN booklet_design USING(test_id);') %>%
    inner_join(bid_info$bid_tb, by=c('test_id','booklet_id')) %>%
    rename(old_bid='bid') %>%
    mutate(old_bid=as.integer(.data$old_bid))
  
  
  routing = dsg$booklets_modules %>%
    inner_join(routing, by=c('old_bid','module_nbr')) %>%
    inner_join(mod_max, by=c('bid','module_nbr'))  %>%
    arrange(.data$bid, .data$module_nbr) %>%
    group_by(.data$bid) %>%
    do({
      tmp=.
      if(tmp$routing[1] == 'last')
      {
        tmp %>%
          mutate(module_exit_score_min = pmax(0L, .data$module_exit_score_min - .data$mod_excluded),
                 module_exit_score_max = pmin(.data$mod_maxscore, .data$module_exit_score_max - .data$mod_excluded))
      } else
      {
        # pmax(0, prblbly superfluous, you cannot get smaller than 0 imo
        tmp = tmp %>%
          mutate(module_exit_score_min = pmax(0L, .data$module_exit_score_min - cumsum(.data$mod_excluded)),
                 module_exit_score_max =  .data$module_exit_score_max - cumsum(.data$mod_excluded)) 
        
        tmp$module_exit_score_max[1] = min(tmp$mod_maxscore[1],tmp$module_exit_score_max[1])
        if(nrow(tmp)>1)
        {
          for(i in 2:nrow(tmp))
          {
            tmp$module_exit_score_max[i] = pmin(tmp$module_exit_score_max[i], 
                                                tmp$module_exit_score_max[i-1] + tmp$mod_maxscore[i])
          }
        }
        tmp
      }
    }) %>%
    ungroup() %>%
    select(-.data$old_bid,-.data$mod_maxscore)

  bid_lev = routing %>%
    group_by(.data$bid, .data$test_id, .data$booklet_id) %>%
    summarise(pf = paste(.data$mod_excluded, collapse='-')) %>%
    ungroup()
  
  bid_lev$str_id = paste0(bid_lev$test_id,'-', bid_lev$booklet_id, ' (',bid_lev$pf,')')
  if(anyDuplicated(bid_lev$str_id))
  {
    tchar=max(nchar(bid_lev$test_id))
    bid_lev$str_id = sprintf(paste0('% -',tchar,'s-%s (%s)'),bid_lev$test_id, bid_lev$booklet_id, bid_lev$pf)
  }
  
  
  
  class(bdes$bid) = 'factor'
  class(routing$bid) = 'factor'
  class(suf_stats$plt$bid) = 'factor'
  class(respData$bid) = 'factor'
  
  levels(bdes$bid) = bid_lev$str_id 
  levels(routing$bid) = bid_lev$str_id 
  levels(suf_stats$plt$bid) = bid_lev$str_id 
  levels(respData$bid) = bid_lev$str_id 
  
  list(x = respData, 
       booklet_design= bdes,
       routing = routing,
       suf_stats = suf_stats)
  
 }  
