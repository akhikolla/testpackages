globalVariables(c("."))


#to~do: is het mogelijk dat check op onmogelijke toetspaden een warning wordt?


#' Define routing rules
#'
#' Define routing rules for use in \code{\link{create_mst_test}}
#'
#' @param ... routing rules defined using a a dot-like syntax, read --+ as an arrow and [:] as a range of score to move to the next stage
#' @return data.frame with columns...
#' @details 
#' Each scoring rule in `...` defines one or more routing rules together making up a booklet.
#' For example, `route1 = a[0:5] --+ d[9:15] --+ f` means a start at module `a`, continue to module `d` when the score on 
#' `a` is between 0 and 5 (inclusive) and continue to `g` when the score on modules `a + b` is between 0 and 8 (for `All` routing)
#' or the score on just module 'b' is between 0 and 8 (for `Last` routing).
#' `route1` becomes the id of the specific path or booklet, which must be supplied with the data later.
#' 
#' A routing design for a linear (non-multistage) booklet can simply be entered as \code{mst_rules(my_booklet = my_single_module)}.
#' 
#' @seealso
#' \code{\link{create_mst_test}} for a description of all and last routing and \code{\link{add_response_data_mst}} to see how to enter data
#' 
#' @examples
#' # a (complicated) three stage (1-3-3) routing design with 9 booklets and 7 modules
#' 
#' routing_rules = mst_rules(bk1 = M1[0:61] --+ M2[0:136]   --+ M5,
#'                           bk2 = M1[0:61] --+ M2[137:183] --+ M6,
#'                           bk3 = M1[0:61] --+ M2[184:Inf] --+ M7,
#'
#'                           bk4 = M1[62:86] --+ M3[0:98]    --+ M5,
#'                           bk5 = M1[62:86] --+ M3[99:149]  --+ M6,
#'                           bk6 = M1[62:86] --+ M3[150:Inf] --+ M7,
#'
#'                           bk7 = M1[87:Inf] --+ M4[0:98]    --+ M5,
#'                           bk8 = M1[87:Inf] --+ M4[99:130]  --+ M6,
#'                           bk9 = M1[87:Inf] --+ M4[131:Inf] --+ M7)
#'       
mst_rules = function(...)
{
  mf = as.list(match.call())[-1]
  
  gop = function(x) {
    if (is.call(x)) {
      return(list(as.character(x[[1]]), lapply(x[-1], gop)))
    }
    else {
      return(NULL)
    }
  }
  ops = unlist(lapply(mf, gop))
  
  if (!all(ops %in% c("-", ":","[","+"))) {
    stop("Invalid operator in formula")
  }
  
  f = function(x) {
    if (is.call(x)) {
      if (length(x) == 3) {
        if(as.character(x[[1]]) == '[')
        {
          return(list(f(x[[2]]), op = as.character(x[[1]]), 
                      f(x[[3]]), op = ']'))
        }
        return(list(f(x[[2]]), op = as.character(x[[1]]), 
                    f(x[[3]])))
      }
      else {
        return(list(op = as.character(x[[1]]), f(x[[2]])))
      }
    }
    else {
      return(c(sym = as.character(x)))
    }
  }
  gj = function(x)
  {
    x = t(x[c(TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE)])

    do.call(bind_rows,
            lapply(split(x, ceiling(seq_along(x)/3)),
                   function(y){tibble(module_id = y[1], 
                                      exit_min = suppressWarnings(as.integer(y[2])), 
                                      exit_max = suppressWarnings(as.integer(y[3])))})) %>%
      mutate(module_nbr = row_number())
  }
  
  booklets =  mf %>%             
    lapply(function(x) unlist(f(x))) %>%
    lapply(check_rule) %>%  
    lapply(gj)

  if(any(duplicated(names(booklets))))
  {
    stop("booklet id's are not unique")
  }
  
  for(bk_name in names(booklets))
  {
    booklets[[bk_name]]$booklet_id = bk_name
  }

  res = do.call(rbind, booklets)
  check_routing(res)
  res
} 


#' open an existing mst project
#' 
#' @param pth path to project file
#' 
open_mst_project = function(pth)
{
  if (!file.exists(pth)) stop("There is no such file")
  db = dbConnect(SQLite(), pth)
  dbExecute(db, 'pragma foreign_keys=1;')
  db
}

#' create a new (empty) mst project
#' 
#' @param pth path and filename to save project file
#' 
#' @return handle to project database
#' 
create_mst_project = function( pth)
{
  if (file.exists(pth)) file.remove(pth)
  db = dbConnect(SQLite(), pth)

  dbRunScript(db, 'mst_sqlite.sql')
  db
}

#' Close an mst project
#'
#' @param db dextermst project db connection
#' 
close_mst_project = function(db) dbDisconnect(db) 


#' add scoring rules to an mst project
#' 
#' @param db a dextermst db connection
#' @param rules dataframe (item_id, response, item_score), 
#' listing all permissible responses to an item and their scores
#' 
add_scoring_rules_mst = function(db, rules)
{
  rules$item_id = as.character(rules$item_id)
  rules$response = as.character(rules$response)
  rules$item_score = as.integer(rules$item_score)
  
  problems = rules %>%
    union_all(dbGetQuery(db, 'SELECT item_id, response, item_score FROM Scoring_rules;')) %>%
    group_by(.data$item_id) %>%
    summarise(min_score = min(.data$item_score), 
              distinct_scores = n_distinct(.data$item_score),
              duplicated_responses = any(duplicated(.data$response))) %>%
    ungroup() %>%
    filter(.data$min_score > 0 | .data$distinct_scores < 2 | .data$duplicated_responses)
  
  if(nrow(problems) > 0)
  {
    message('problematic scoring rules, check output')
    print(as.data.frame(problems), row.names=FALSE)
    stop('invalid rules')
  }
  new_items = rules %>% 
    distinct(.data$item_id) %>%
    anti_join(dbGetQuery(db, 'SELECT item_id FROM Items;'), by='item_id')
    
  dbTransaction(db, {
    dbExecute(db, 'INSERT INTO Items(item_id) VALUES(:item_id);', new_items)
    n = dbExecute(db, 'INSERT INTO Scoring_rules(item_id, response, item_score) VALUES(:item_id, :response, :item_score);', rules)
  })
  cat(paste('Added', n, 'scoring rules\n'))
  invisible(NULL)
}



#' alter scoring rules in an mst project
#' 
#' @param db a dextermst db connection
#' @param rules data.frame (item_id, response, item_score), see dexter
#' 
#' @description 
#' It is only possible to change item_scores for existing items and responses 
#' through this function. Scoring rules can only be changed for items that are 
#' in the last module of a (mst) test.
#' 
alter_scoring_rules_mst = function(db, rules)
{
  rules$response = as.character(rules$response)
  rules$item_id = as.character(rules$item_id)
  rules$item_score = as.integer(rules$item_score)
  rules = rules[,c('item_id', 'response','item_score')] 
  
  updated_items = dbGetQuery(db, 
                             'SELECT item_id, response, item_score FROM Scoring_rules WHERE item_id=:item_id;',
                             distinct(rules, .data$item_id))
  
  problems = updated_items %>%
    anti_join(rules, by=c('item_id','response')) %>%
    union(rules) %>%
    group_by(.data$item_id) %>%
    summarise(min_score = min(.data$item_score), 
              distinct_scores = n_distinct(.data$item_score),
              duplicated_responses = any(duplicated(.data$response))) %>%
    ungroup() %>%
    filter(.data$min_score > 0 | .data$distinct_scores < 2 | .data$duplicated_responses)
  
  if(nrow(problems) > 0)
  {
    message('problematic scoring rules, check output')
    print(as.data.frame(problems))
    stop('invalid rules')
  }
  
  mdl = dbGetQuery(db, 
    'WITH mdl_cnt AS (
      SELECT test_id, booklet_id, COUNT(*) AS n_modules
        FROM Booklet_design
        GROUP BY test_id, booklet_id)
    SELECT test_id, booklet_id, module_id, module_nbr, item_id
        FROM booklet_design
          INNER JOIN module_design USING(test_id, module_id)
            INNER JOIN mdl_cnt USING(test_id, booklet_id)
        WHERE n_modules > module_nbr AND item_id=:item_id;', 
    distinct(rules, .data$item_id))
  
  if(nrow(mdl) > 0)
  {
    message(paste('scoring rules cannot be updated for this mst design because they',
                  'influence items that are not in the last module of a booklet, check output'))
    print(mdl)
    stop()
  }
  dbTransaction(db, {
  
    n = dbExecute(db, 
                  'UPDATE Scoring_rules SET item_score=:item_score WHERE response=:response AND item_id=:item_id;',
                  rules)
    
    new_max = rules %>%
      group_by(.data$item_id) %>%
      summarise(item_newmax = sum(.data$item_score)) %>%
      ungroup()
    
    bd_update = updated_items %>% 
      group_by(.data$item_id) %>%
      summarise(item_oldmax = max(.data$item_score)) %>%
      ungroup() %>%
      inner_join(new_max, by='item_id') %>%
      mutate(dif = .data$item_newmax - .data$item_oldmax) %>%
      inner_join(dbGetQuery(db, 'SELECT test_id, module_id, item_id FROM Module_design;'), 
                 by='item_id') %>%
      group_by(.data$test_id, .data$module_id) %>%
      summarise(module_dif = sum(.data$dif)) %>%
      ungroup() %>%
      inner_join(dbGetQuery(db, 'SELECT test_id, booklet_id, module_id, module_exit_score_max
                                  FROM Booklet_design;'), 
                 by=c('test_id','module_id')) %>%
      mutate(new_max = .data$module_exit_score_max + .data$module_dif) %>%
      filter(.data$module_exit_score_max != .data$new_max)
    
    dbExecute(db, 
              'UPDATE Booklet_design SET module_exit_score_max = :new_max
                   WHERE test_id = :test_id AND booklet_id = :booklet_id AND module_id = :module_id',
              bd_update)
  })
  cat(paste(n, 'rules updated\n'))
  invisible(NULL)
}


#' Define a new multi stage test
#' 
#' Before you can enter data, dexterMST needs to know the design of your test. 
#'  
#' @param db output of \code{\link{open_mst_project}} or \code{\link{create_mst_project}}
#' @param test_design data.frame with columns item_id, module_id, item_position
#' @param routing_rules output of \code{\link{mst_rules}} 
#' @param test_id id of the mst test
#' @param routing all or last routing (see details)
#' 
#' @details
#' In dexterMST we use the following terminology:
#' \describe{
#' \item{test}{collection of modules and rules to go from one module to the other.
#' A test must have one starting module}
#' \item{booklet}{a specific path through a mst test.}
#' \item{module}{ a block of items that is always administered together. Each item has a specific position in a module.}
#' \item{routing rules}{rules to go from one module to another based on score on the current and possibly previous modules}
#'  }
#'  
#'  Additionally, there are two possible types of routing:
#' \describe{
#' \item{all}{the routing rules are based on the sum of the current and previous modules}
#' \item{last}{the routing rules are based only on the current module}
#' }
#' The type of routing must be defined for a test as a whole so it is not possible to mix routing types. 
#' In CML (as opposed to MML) the routing rules are actually used in the calibration so it
#' is important they are correctly specified. DexterMST includes multiple checks, 
#' both when defining the test and when entering data, 
#' to make sure your routing rules are valid and your data conform to them. 
#' 
#' 
#' @examples
#' # extended example
#' # we: 
#' # 1) define an mst design
#' # 2) simulate mst data
#' # 3) create a project, enter scoring rules and define the MST test
#' # 4) do an analysis
#' 
#' library(dplyr)
#' 
#' items = data.frame(item_id=sprintf("item%02i",1:70), item_score=1, delta=sort(runif(70,-1,1)))
#'
#' design = data.frame(item_id=sprintf("item%02i",1:70),
#'                     module_id=rep(c('M4','M2','M5','M1','M6','M3', 'M7'),each=10))
#'
#' routing_rules = routing_rules = mst_rules(
#'  `124` = M1[0:5] --+ M2[0:10] --+ M4, 
#'  `125` = M1[0:5] --+ M2[11:15] --+ M5,
#'  `136` = M1[6:10] --+ M3[6:15] --+ M6,
#'  `137` = M1[6:10] --+ M3[16:20] --+ M7)
#'
#' theta = rnorm(3000)
#'
#' dat = sim_mst(items, theta, design, routing_rules,'all')
#' dat$test_id='sim_test'
#' dat$response=dat$item_score
#' 
#' 
#' scoring_rules = data.frame(
#'   item_id = rep(items$item_id,2), 
#'   item_score= rep(0:1,each=nrow(items)),
#'   response= rep(0:1,each=nrow(items))) # dummy respons
#'   
#' 
#' db = create_mst_project(":memory:")
#' add_scoring_rules_mst(db, scoring_rules)
#'
#' create_mst_test(db,
#'                 test_design = design,
#'                 routing_rules = routing_rules,
#'                 test_id = 'sim_test',
#'                 routing = "all")
#' 
#' add_response_data_mst(db, dat)
#' 
#' 
#' design_plot(db)
#' 
#' f = fit_enorm_mst(db)
#' 
#' head(coef(f))
#' 
#' abl = ability(get_responses_mst(db), f) %>%
#'    inner_join(tibble(person_id=as.character(1:3000), theta.sim=theta), by='person_id')
#' 
#' plot(abl$theta, abl$theta.sim)   
#' 
#' abl = filter(abl, is.finite(theta))
#' 
#' cor(abl$theta, abl$theta.sim)
#' 
create_mst_test = function(db, test_design, routing_rules, test_id, routing = c('all', 'last'))
{
  routing = match.arg(routing)
  
  if(!'item_position' %in% test_design)
  {
    test_design = test_design %>%
      group_by(.data$module_id) %>%
      mutate(item_position=row_number()) %>%
      ungroup()
  }
  
  test_design = test_design %>% 
    select(.data$item_id, .data$module_id, .data$item_position) %>%
    mutate_if(is.factor, as.character)
  
  routing_rules = routing_rules %>%
    select(.data$module_id, .data$exit_min, .data$exit_max, .data$module_nbr, .data$booklet_id) %>%
    mutate_if(is.factor, as.character)
  
  check_routing(routing_rules)
  
  if(!setequal(routing_rules[['module_id']], test_design[['module_id']]))
    stop("module id's in routing_rules differ from module_id's in test_design")
  
  unknown_items = test_design %>%
    distinct(.data$item_id) %>%
    anti_join(dbGetQuery(db, 'SELECT item_id FROM Items;'), by='item_id')
  
  if(nrow(unknown_items) > 0)
  {
    message('Unknown items, showing first 30')
    print(pull(unknown_items, 'item_id')[1:30])
    stop(paste(nrow(unknown_items), 'items in your test design have not been defined in the scoring rules.',
               'Add scoring rules first using the function `add_scoring_rules_mst`'))
  }
  
  itm_max = dbGetQuery(db, 'SELECT item_id, MAX(item_score) AS item_maxscore 
                              FROM Scoring_rules GROUP BY item_id;')
  
  module_max = test_design %>%
    inner_join(itm_max, by='item_id') %>%
    group_by(.data$module_id) %>%
    summarise(mod_max = sum(.data$item_maxscore)) %>%
    ungroup()

  if(routing=='all')
  {
    correct_exit_max = function(exit_max, mod_max)
    {
      for(i in seq_along(exit_max))
      {
        prevmax = ifelse(i==1, 0L, exit_max[i-1])
        if(is.na(exit_max[i]) || exit_max[i] > mod_max[i] + prevmax)
        {
          exit_max[i] = mod_max[i] + prevmax
        } 
      }
      exit_max
    }
    
    bd = routing_rules %>%
      inner_join(module_max, by='module_id') %>%
      mutate(old_min = coalesce(.data$exit_min,0L), old_max = .data$exit_max) %>%
      arrange(.data$booklet_id, .data$module_nbr) %>%
      group_by(.data$booklet_id) %>%
      mutate(exit_min = cummax(coalesce(.data$exit_min,0L)),
             exit_max = correct_exit_max(.data$exit_max, .data$mod_max)) %>%
      ungroup()
  } else
  {
    bd = routing_rules %>%
      inner_join(module_max, by='module_id') %>%
      mutate(old_min = coalesce(.data$exit_min,0L), old_max = .data$exit_max,
              exit_min = coalesce(.data$exit_min,0L), 
              exit_max = pmin(coalesce(.data$exit_max, .data$mod_max), .data$mod_max))
  }
  
  adjustments = bd %>%
    arrange(.data$booklet_id, .data$module_nbr) %>%
    group_by(.data$booklet_id) %>%
    mutate(nmod=n()) %>%
    ungroup() %>%
    mutate(adj_min = (.data$nmod != .data$module_nbr & .data$old_min != .data$exit_min),
           adj_max = (.data$nmod != .data$module_nbr & coalesce(.data$old_max,.data$exit_max) != .data$exit_max)) %>%
    group_by(.data$booklet_id) %>%
    filter(any(.data$adj_min | .data$adj_max)) %>%
    ungroup()
  
  if(nrow(adjustments) > 0)
  {
    message('One or more of your routing rules have been adjusted based on the maximum ',
            'possible scores in the modules given your item scoring rules.')
    
    adjustments %>%
      mutate(exit_min = as.character(.data$exit_min), exit_max = as.character(.data$exit_max)) %>%
      mutate(exit_min = if_else(.data$adj_min, red(.data$exit_min), .data$exit_min),
             exit_max = if_else(.data$adj_max, red(.data$exit_max), .data$exit_max)) %>%
      arrange(.data$module_nbr) %>%
      group_by(.data$booklet_id) %>%
      summarise(msg = paste0('`',.data$booklet_id[1], '` = ', 
                             paste0(.data$module_id, 
                              if_else(.data$module_nbr < .data$nmod,
                                      paste0('[',.data$exit_min,':',.data$exit_max,']'),'\n'),
                              collapse = ' --+ '))) %>%
      ungroup() %>%
      pull(.data$msg) %>%
      paste(collapse='') %>%
      cat()
  }
    
 
  if(any(bd$exit_min > bd$exit_max))
  {
    message('The routing designs for the folowing booklets are impossible')
    
    problems = bd %>%
      group_by(.data$booklet_id) %>%
      filter(any(.data$exit_min > .data$exit_max)) %>%
      summarise(rtng = paste0(.data$module_id, '[',.data$exit_min,':',
                                .data$exit_max,']', collapse = ' --+ ')) %>%
      ungroup()
      
    cat(paste0(problems$booklet_id, ': ', problems$rtng, collapse='\n'))
    stop('impossible routing')
  }

  
  dbTransaction(db,{
    dbExecute(db,'INSERT INTO Tests(test_id, routing) VALUES(:test_id,:routing);', 
                  tibble(test_id=test_id, routing=routing))
    
    dbExecute(db,'INSERT INTO Modules(test_id, module_id) VALUES(:test_id,:module_id);',
                  tibble(test_id=test_id, module_id = unique(bd$module_id)))
    
    dbExecute(db, 'INSERT INTO Booklets(test_id, booklet_id) VALUES(:test_id, :booklet_id);',
              tibble(test_id=test_id, booklet_id = unique(bd$booklet_id)))
    
    dbExecute(db, 'INSERT INTO Module_design(test_id,	module_id, item_id,	item_position)
                  VALUES(:test_id, :module_id, :item_id,	:item_position);',
              test_design %>%
                mutate(test_id = test_id) %>% 
                semi_join(bd, by='module_id') %>%
                select(.data$test_id, .data$module_id, .data$item_id,	.data$item_position))
    bd$test_id = test_id
    dbExecute(db, 
        'INSERT INTO Booklet_design(test_id, booklet_id,
    	                              module_id, module_nbr, module_exit_score_min,	module_exit_score_max) 
          VALUES(:test_id, :booklet_id, :module_id, :module_nbr, :exit_min, :exit_max);',
        select(bd,.data$test_id, .data$booklet_id, .data$module_id, .data$module_nbr, .data$exit_min, .data$exit_max))
  })  
 invisible(NULL)
}


#' Add multistage response data
#' 
#' Multistage response data can be entered in long format for one or multiple booklets simultaneously
#' or in wide format one booklet at a time. 
#' 
#' @param db a dextermst db handle
#' @param rsp_data data.frame with columns (person_id, test_id, booklet_id, item_id, response)
#' @param booklet_data data.frame with a column person_id and other columns which names correspond to item_id's
#' @param test_id id of a test known in the database
#' @param booklet_id id of a booklet known in the database
#' @param auto_add_unknown_rules if FALSE, unknown responses (i.e. not defined in the scoring rules) will 
#' generate an error and the function will abort. If TRUE unknown responses will be automatically added to 
#' the scoring rules with a score of 0
#' 
#' @details
#' Users familiar with dexter might expect to be able to enter new booklets here. Because mst tests 
#' have a more complicated design that cannot be (easily) derived from the data, 
#' in dexterMST the test designs have to be entered beforehand. 
#' 
#' @seealso
#' \code{\link{create_mst_test}}
#' 
add_response_data_mst = function(db, rsp_data, auto_add_unknown_rules = FALSE)
{
  if(length(setdiff(c('person_id', 'test_id', 'booklet_id', 'item_id', 'response'), colnames(rsp_data)))> 0) 
    stop('columns (person_id, test_id, booklet_id, item_id, response) are required')
  
  rsp_data = rsp_data %>% 
    select(.data$person_id, .data$test_id, .data$booklet_id, .data$item_id, .data$response)
  
  md = dbGetQuery(db,'SELECT test_id, booklet_id, module_id, item_id 
                        FROM Booklet_design INNER JOIN Module_design USING(test_id, module_id);')
  
  rsp_data = rsp_data %>% 
    mutate_if(is.factor, as.character) %>%
    mutate(response = as.character(.data$response), person_id=as.character(.data$person_id)) %>%
    left_join(md, by=c('test_id','booklet_id', 'item_id'))
  
  rsp_data[is.na(rsp_data$response),'response'] = 'NA'
  
  # design errors
  if(anyNA(rsp_data$module_id))
  {
    message('Some booklet-item combinations are unknown (showing first 30)')
    rsp_data %>%
      filter(is.na(.data$module_id)) %>%
      slice(1:30) %>%
      as.data.frame() %>%
      print(row.names=FALSE)
    stop('rsp_data does not match your test design')
  }
    
  unknown_responses = rsp_data %>%
    distinct(.data$item_id, .data$response) %>%
    anti_join(dbGetQuery(db,'SELECT item_id, response FROM scoring_rules;'), by=c('item_id','response'))
  
  dbTransaction(db,
  {
    if(nrow(unknown_responses)>0)
    {
      if(auto_add_unknown_rules)
      {
        dbExecute(db, 
                  'INSERT INTO Scoring_rules(item_id, response, item_score) VALUES(:item_id, :response, 0);',
                  unknown_responses)
      } else
      {
        message(paste('Some responses do not occur in your scoring rules (showing first 30).',
                      'Set auto_add_unknown_rules to TRUE to code them as zero automatically'))
        unknown_responses %>%
          arrange(.data$item_id, .data$response) %>%
          slice(1:30) %>%
          as.data.frame(row.names=FALSE) %>%
          print()
        stop('unknown responses, see output')
      }
    }
    
    booklet_design = dbGetQuery(db, 
                                'SELECT *
                                  FROM Tests INNER JOIN Booklet_design USING(test_id)')
    
    routing_errors = rsp_data %>%
      inner_join(dbGetQuery(db,'SELECT item_id, response, item_score FROM scoring_rules;'), by=c('item_id','response')) %>%
      group_by(.data$person_id, .data$test_id, .data$booklet_id, .data$module_id) %>%
      summarise(sum_score = sum(.data$item_score)) %>%
      ungroup() %>%
      inner_join(booklet_design, by = c('test_id', 'booklet_id', 'module_id')) %>%
      arrange(.data$person_id, .data$test_id, .data$booklet_id, .data$module_nbr) %>%
      group_by(.data$person_id, .data$test_id, .data$booklet_id) %>%
      mutate(exit_score = if_else(.data$routing=='all',cumsum(.data$sum_score), .data$sum_score)) %>%
      ungroup() %>%
      filter(.data$exit_score < .data$module_exit_score_min | .data$exit_score > .data$module_exit_score_max)
  
    if(nrow(routing_errors) > 0)
    {
      message('encountered uneligible scores for mst routing design (showing first 30)')
      print(as.data.frame(slice(routing_errors, 1:30)), row.names=FALSE)
      stop('invalid scores')
    }
    
    bkn = dbGetQuery(db, 'SELECT test_id, booklet_id, COUNT(*) AS n_responses_required
                        FROM Booklet_design 
                          INNER JOIN Module_design USING(test_id, module_id)
                      GROUP BY test_id, booklet_id;')
    
    missing_data_errors = rsp_data %>%
      group_by(.data$test_id, .data$booklet_id, .data$person_id) %>%
      summarise(n_responses = n()) %>%
      ungroup() %>%
      inner_join(bkn, by=c('test_id','booklet_id')) %>%
      filter(.data$n_responses != .data$n_responses_required)
    
    if(nrow(missing_data_errors) > 0)
    {
      message('too few responses for some students (showing first 30)')
      print(as.data.frame(slice(missing_data_errors, 1:30)),row.names=FALSE)
      stop('missing data is not allowed')
    }
    
    prs = rsp_data %>% 
      distinct(.data$person_id) %>%
      anti_join(dbGetQuery(db, 'SELECT person_id FROM persons;'), by='person_id')
    
    
    if(nrow(prs) > 0)
    {
      dbExecute(db, 'INSERT INTO Persons(person_id) VALUES(:person_id);', prs)
    }
    
    n1 = dbExecute(db, 'INSERT INTO Administrations(person_id, test_id, booklet_id) VALUES(:person_id, :test_id, :booklet_id);',
                    distinct(rsp_data, .data$person_id, .data$test_id, .data$booklet_id))
    
    n2 = dbExecute(db, 'INSERT INTO Responses(person_id, test_id, booklet_id, module_id, item_id, response)
                        VALUES(:person_id, :test_id, :booklet_id, :module_id, :item_id, :response);', 
                  select(rsp_data, .data$person_id, .data$test_id, .data$booklet_id, .data$module_id, .data$item_id, .data$response))
  })
  cat(paste('Added', n2, 'responses for', n1, 'administrations\n'))
}

#' @rdname add_response_data_mst
add_booklet_mst = function(db, booklet_data, test_id, booklet_id, auto_add_unknown_rules = FALSE)
{
  if(!'person_id' %in% colnames(booklet_data))
    stop('booklet data must contain a column person id')
  
  if(any(duplicated(booklet_data$person_id)))
    stop('column person_id must contain unique values')
  
  rsp_data = gather(booklet_data, key='item_id', value='response', -.data$person_id)
  
  rsp_data$test_id = test_id
  rsp_data$booklet_id = booklet_id
  
  add_response_data_mst(db, rsp_data, auto_add_unknown_rules)
}



#'retrieve information from a mst database
#'
#'@param db dexterMST project database connection
#'
get_booklets_mst = function(db)
{
  dbGetQuery(db,
    "WITH bkn AS (SELECT test_id, booklet_id, COUNT(*) AS n_respondents FROM Administrations GROUP BY test_id, booklet_id),
          bki AS (SELECT test_id, booklet_id, COUNT(*) AS n_items, MAX(module_nbr) AS n_modules 
                    FROM Booklet_design INNER JOIN Module_design USING(test_id, module_id)
                      GROUP BY test_id, booklet_id)
     SELECT test_id, booklet_id, n_modules, n_items, COALESCE(n_respondents,0) AS n_respondents
      FROM bki 
        LEFT OUTER JOIN bkn USING(test_id, booklet_id);")
}

#' @rdname get_booklets_mst
get_design_mst = function(db)
{
  dbGetQuery(db, 
    'SELECT test_id, booklet_id, module_id, item_id 
      FROM Booklet_design 
        INNER JOIN Module_design USING(test_id, module_id)
    ORDER BY test_id, booklet_id, module_nbr, item_position;') %>%
    group_by(.data$test_id, .data$booklet_id)  %>%
    mutate(item_position = row_number()) %>%
    ungroup()
}

#' @rdname get_booklets_mst
get_routing_rules_mst = function(db)
{
  dbGetQuery(db, 
             'SELECT * FROM booklet_design INNER JOIN Tests USING(test_id)
             ORDER BY test_id, booklet_id, module_nbr;')
}

#' @rdname get_booklets_mst       
get_scoring_rules_mst = function(db)
{
  dbReadTable(db, 'Scoring_rules')
}

#' @rdname get_booklets_mst  
get_items_mst = function(db)
{
  dbReadTable(db, 'Items')
}

#' @rdname get_booklets_mst  
get_persons_mst = function(db)
{
  dbGetQuery(db,'SELECT * FROM Persons;')
}

#' Add item properties to an dextermst project
#'
#' @param db dexterMST project database
#' @param item_properties data.frame with a column item_id and other columns containing the item properties
#'   
add_item_properties_mst = function(db, item_properties)
{
  item_properties = item_properties %>% mutate_if(is.factor, as.character)
  colnames(item_properties) = tolower(colnames(item_properties))
  
  if(!'item_id' %in% colnames(item_properties))
    stop('item_properties must contain a column named item_id')
  
  for(nprop in setdiff(colnames(item_properties), dbListFields(db, 'items')))
  {
    dbExecute(db, paste('ALTER TABLE Items ADD COLUMN', nprop, 
                        sql_data_type(item_properties[[nprop]]),
                        'DEFAULT NULL;'))
  }
  iprop = setdiff(colnames(item_properties), 'item_id')
  dbTransaction(db,{
    n = dbExecute(db, paste0('UPDATE Items SET ', paste0(iprop,'=:',iprop, collapse=','),
                         ' WHERE item_id=:item_id;'),
                item_properties)
  })
  cat(paste(ncol(item_properties)-1, 'item properties updated for', n, 'items'))
  
}

#' Add person properties to a mst project
#'
#' @param db dextermst project database
#' @param person_properties data.frame with a column person_id and 
#' other columns containing the person properties
#'   
add_person_properties_mst = function(db, person_properties)
{
  person_properties = person_properties %>% mutate_if(is.factor, as.character)
  colnames(person_properties) = tolower(colnames(person_properties))
  
  if(!'person_id' %in% colnames(person_properties))
    stop('person_properties must contain a column named person_id')
  
  blacklist = unique(unlist(lapply(c('Scoring_rules','Items', 'Administrations', 
                                     'Responses', 'module_design', 'booklet_design'),
                                   dbListFields, conn=db)))
  blacklist = blacklist[blacklist != 'person_id']
  
  if(length(intersect(colnames(person_properties), blacklist)) > 0)
    stop(paste('columns',
               paste(intersect(colnames(person_properties), blacklist), collapse=','),
               'cannot be used as person properties'))
  
  
  for(nprop in setdiff(colnames(person_properties), dbListFields(db, 'persons')))
  {
    dbExecute(db, paste('ALTER TABLE persons ADD COLUMN', nprop, 
                        sql_data_type(person_properties[[nprop]]),
                        'DEFAULT NULL;'))
  }
  pprop = colnames(person_properties)[colnames(person_properties) != 'person_id']

  dbTransaction(db,{
    n = dbExecute(db, paste0('UPDATE persons SET ', paste0(pprop,'=:',pprop, collapse=','),
                          ' WHERE person_id=:person_id;'),
                person_properties)
  })
  cat(paste(ncol(person_properties)-1, 'person properties updated for', n, 'persons'))
}
