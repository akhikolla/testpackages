dbRunScript <- function(db, fn)
{
  # run sql script included in the package
  # The R dbi api does not provide for execution of scripts.
  # A usual method to get around this is to split the input script on ;
  # however, there are numerous exceptions where this would not work (e.g. strings, triggers)
  # the kludge is to use a custom split string in our sql scripts, which is:
  # --#split#--
  
  fn = system.file("extdata", fn, package = "dexterMST", mustWork = TRUE)
  
  script = strsplit(paste0(readLines(fn, warn = FALSE), collapse='\n'),'--#split#--')[[1]]
  
  for (statement in script)
  {
    dbExecute(db,statement)
  }
}


dbTransaction = function(db, expr, on_error = stop, on_error_rollback=TRUE)
{
  if(is(db, 'SQLiteConnection')) dbExecute(db,'pragma foreign_keys=1;')
  dbBegin(db)
  tryCatch(expr, error=function(e){if(on_error_rollback) dbRollback(db); on_error(e);}, finally=NULL)
  tryCatch(dbCommit(db), error=function(e){if(on_error_rollback) dbRollback(db); on_error(e);}, finally=NULL)
}




sql_coldef_from_pragma = function(pr_info, tbname)
{
  a = pr_info %>%
    mutate(stm=
      paste('ALTER TABLE', tbname, 'ADD COLUMN',
            .data$name,
            .data$type,
            if_else(.data$notnull==1 | !is.na(.data$dflt_value),'NOT NULL',''),
            case_when(is.na(.data$dflt_value) ~ '',
                      grepl('null',.data$dflt_value, ignore.case=TRUE) ~ 'DEFAULT NULL',
                      grepl('(int)|(double)|(real)',.data$type, ignore.case=TRUE) ~ 
                        paste('DEFAULT',.data$dflt_value),
                      TRUE ~ paste('DEFAULT',sql_quote(.data$dflt_value, "'"))),
            ';') ) %>%
    pull('stm')
  
  b = pr_info %>%
    filter(grepl('(int)|(double)|(real)',.data$type, ignore.case=TRUE) & !is.na(.data$dflt_value)) %>%
    mutate(stm = paste('UPDATE',tbname,'SET', .data$name,'=',.data$name, '+ 0;')) %>%
    pull('stm')
  
  c(a,b)
  
}

#' import data from a dexter project
#' 
#' This function will import items, scoring rules, persons, test designs and responses from
#' a dexter database into the dexterMST database. 
#' 
#' @param db dextermst project db connection
#' @param dexter_db path to a dexter database file or open dexter db connection
#' @param dx_response_prefix string to prefix responses from dexter with (usually not necessary, see details)
#' 
#' @details
#' DexterMST has no problem calibrating data from linear tests. However, dexter and dexterMST have 
#' differently structured project databases. If you already have response data from linear tests in
#' a dexter database, you can easily import it into your dexterMST database from there.
#' 
#' The dexterMST variables test_id, module_id and booklet_id will all be set to the dexter variable
#' booklet_id (i.e. a linear test becomes a multistage test with one booklet and one module only).
#' 
#' It is assumed that items with equal id's in your dexter and dexterMST project refer to the same items. 
#' If an item in dexter has different score categories compared to an existing item with the same item_id in dexterMST
#' an error will be generated. If the same response to the same item has a different score, this will also generate an error.
#' However, it is possible for an item in dexter to have scoring rules for responses not defined in dexterMST and vice versa.
#' 
#' In the unusual and unfortunate situation that the same response to the same item should have a different score
#' in dexter than in dexterMST, you can use the parameter dx_response_prefix to prefix the responses in dexter with
#' some unique combination of characters, e.g. "dexter". In practice this sometimes happens when old archived data
#' is only available in scored form (i.e. response 0 has score 0, response 1 has score 1) and new data is available in
#' raw form but the actual response can also be 0 or 1, etc. causing a conflict.
#' 
#' @examples
#' \dontrun{
#' library(dexter)
#' dbDex = start_new_project(verbAggrRules, "verbAggression.db", 
#'   person_properties=list(gender="unknown"))
#' add_booklet(dbDex, verbAggrData, "agg")
#' add_item_properties(dbDex, verbAggrProperties)
#' db = create_mst_project(':memory:')
#' import_from_dexter(db, dbDex)
#' f_mst = fit_enorm_mst(db)
#' f_dexter = fit_enorm(dbDex)
#' close_mst_project(db)
#' close_project(dbDex)
#' }
#' 
#' 
import_from_dexter = function(db, dexter_db, dx_response_prefix = '' )
{
  if(inherits(dexter_db, 'DBIConnection'))
  {
    dxdb = dexter_db
  } else
  {
    dxdb = open_project(dexter_db)
  }
  
  
  dbTransaction(db, 
  {
    # new iprop
    new_iprop = setdiff(dbListFields(dxdb,'dxItems'), dbListFields(db,'Items'))
    
    # items and properties
    dbGetQuery(dxdb, 'pragma table_info(dxItems);') %>%     
      filter(!.data$name %in% dbListFields(db,'items')) %>%
      sql_coldef_from_pragma(tbname = 'Items') %>%
      lapply(dbExecute, conn = db)
  
    # add new items with properties  
    dxitems = get_items(dxdb) %>%
      anti_join(dbGetQuery(db,'SELECT item_id FROM Items;'), by='item_id')
    
    dbExecute(db, 
              paste0('INSERT INTO Items(', 
                     paste(colnames(dxitems), collapse=','), 
                     ') VALUES(',
                     paste0(':',colnames(dxitems), collapse=','),
                    ');'),
              dxitems)
    
    if(length(new_iprop) > 0)
    {
      existing_items =  get_items(dxdb) %>%
        inner_join(dbGetQuery(db,'SELECT item_id FROM Items;'), by='item_id')

      dbExecute(db, 
                paste0('UPDATE Items SET ', 
                       paste0(new_iprop, '=:',new_iprop, collapse=','), 
                       ' WHERE item_id=:item_id;'),
                       existing_items)
    }
    
    
    
    dxrules = get_rules(dxdb) %>% 
      mutate(response = paste0(dx_response_prefix, .data$response))
    
    mbscores = dbGetQuery(db, 'SELECT DISTINCT item_id, item_score FROM Scoring_rules;') %>%
      semi_join(dxrules, by='item_id') %>%
      add_column(in_mst=1L)

    mismatch = dxrules %>%
      distinct(.data$item_id, .data$item_score) %>%
      semi_join(mbscores, by='item_id') %>%
      add_column(in_dexter = 1L) %>%
      full_join(mbscores, by=c('item_id','item_score')) %>%
      filter(is.na(.data$in_mst) | is.na(.data$in_dexter))
    
    if(nrow(mismatch) > 0)
    {
      cat("score categories for items in dexter that don't have a match in acetcpp:\n\n")
      mismatch %>% 
        group_by(.data$item_id) %>%
        arrange(.data$item_score) %>%
        summarise(scores_in_dexter = paste(.data$item_score[!is.na(.data$in_dexter)], collapse=', '),
                  scores_in_mst = paste(.data$item_score[!is.na(.data$in_mst)], collapse=', ')) %>%
        ungroup() %>%  
        as.data.frame() %>% 
        print(row.names=FALSE)
      stop('score categories in dexter differ from score categories for the same items in acetcpp')
    }
  
    dxrules = dxrules %>%
      anti_join(dbGetQuery(db,'SELECT item_id, response, item_score FROM Scoring_rules;'), 
                by=c('item_id','response','item_score'))
    
    double_responses = dxrules %>%
      inner_join(dbGetQuery(db,'SELECT item_id, response, item_score FROM Scoring_rules;'),
                 by=c('item_id','response'),
                 suffix=c('.dexter','.acetcpp'))
    
    if(nrow(double_responses) > 0)
    {
      cat("\nresponses that are score differently in your dexter and acetcpp project:\n\n")
      double_responses %>%
        arrange( .data$item_id, .data$response) %>%
        as.data.frame() %>% 
        print(row.names=FALSE)
        
      stop('different scores for the same responses')
    }
    
    
    dbExecute(db, 
              'INSERT INTO Scoring_rules(item_id, response, item_score) VALUES(:item_id, :response, :item_score);',
              dxrules)
      
    dxbooklets = dbReadTable(dxdb,'dxBooklets')
    
    dbGetQuery(dxdb, 'pragma table_info(dxBooklets);') %>%
      filter(!.data$name %in% dbListFields(db,'Booklets')) %>%
      sql_coldef_from_pragma(tbname = 'Tests') %>%
      lapply(dbExecute, conn = db)
    
    
    testcol = setdiff(colnames(dxbooklets),'booklet_id')
    if(length(testcol>0))
    {
      dbExecute(db,
                paste0(
                  'INSERT INTO Tests(test_id, ',
                  paste(testcol, collapse=','),
                  ') VALUES(:booklet_id,',
                  paste0(':', testcol, collapse=','),
                  ');'),
                dxbooklets)
    } else
    {
      dbExecute(db,
        'INSERT INTO Tests(test_id) VALUES(:booklet_id);',
        select(dxbooklets, .data$booklet_id))
    }
    
    dbExecute(db, 
              'INSERT INTO Booklets(test_id, booklet_id) VALUES(:booklet_id,:booklet_id);',
              select(dxbooklets,.data$booklet_id))
      
    dbExecute(db, 
              'INSERT INTO Modules(test_id, module_id) VALUES(:booklet_id,:booklet_id);',
              select(dxbooklets,.data$booklet_id))
      
    dxdesign = dbReadTable(dxdb, 'dxBooklet_design')
    
    dbGetQuery(dxdb, 'pragma table_info(dxBooklet_design);') %>%
      filter(!.data$name %in% dbListFields(db,'Module_design')) %>%
      filter(!.data$name == 'booklet_id') %>%
      sql_coldef_from_pragma(tbname = 'Module_design') %>%
      lapply(dbExecute, conn = db)
    
    new_design_columns = setdiff(colnames(dxdesign),c('booklet_id','item_id','item_position'))
    
    dbExecute(db,
              paste0(
              'INSERT INTO Module_design(test_id,module_id,item_id,item_position',
              ifelse(length(new_design_columns)>0,
                     paste0(',',new_design_columns,collapse=''),
                     ''),
              ') VALUES(:booklet_id, :booklet_id, :item_id, :item_position',
              ifelse(length(new_design_columns)>0,
                    paste0(',:',new_design_columns,collapse=''),
                    '')
              ,');'),
              dxdesign)
  
    dxbk = dbGetQuery(dxdb,
               'WITH imax AS (
                 SELECT item_id, MAX(item_score) AS mxi
                  FROM dxScoring_rules
                  GROUP BY item_id)
               SELECT booklet_id, SUM(mxi) AS bkmax
                FROM dxBooklet_design
                  INNER JOIN imax USING(item_id)
                GROUP BY booklet_id;')
    
    dbExecute(db,
              'INSERT INTO Booklet_design
                (test_id, booklet_id, module_id, module_nbr, module_exit_score_min, module_exit_score_max)
                VALUES(:booklet_id, :booklet_id, :booklet_id, 1, 0, :bkmax);',
              dxbk)
        
    # persons and properties
    dxpersons = get_persons(dxdb)

    
    dbGetQuery(dxdb, 'pragma table_info(dxPersons);') %>%
      filter(!.data$name %in% dbListFields(db,'Persons')) %>%
      sql_coldef_from_pragma(tbname = 'Persons') %>%
      lapply(dbExecute, conn = db)
    
    
    
    dbExecute(db, 
              paste0('INSERT INTO Persons(', 
                     paste(colnames(dxpersons), collapse=','), 
                     ') VALUES(',
                     paste0(':',colnames(dxpersons), collapse=','),
                     ');'),
              dxpersons)
    
    dxadmin = dbGetQuery(dxdb, 'SELECT person_id,  booklet_id FROM dxAdministrations;') 
    
    n1 = dbExecute(db,
              'INSERT INTO Administrations(person_id, test_id, booklet_id) 
                VALUES(:person_id, :booklet_id, :booklet_id);',
              dxadmin)
    
    dxrsp  = dbGetQuery(dxdb, 'SELECT person_id, booklet_id, item_id, response FROM dxResponses;') %>%
      mutate(response = paste0(dx_response_prefix, .data$response))
    
    n2 = dbExecute(db, 
              'INSERT INTO Responses(person_id, test_id, booklet_id, module_id, item_id, response)
              VALUES(:person_id, :booklet_id, :booklet_id, :booklet_id, :item_id, :response);',
              dxrsp)
  })
  if(!inherits(dexter_db, 'DBIConnection'))  close_project(dxdb)
    
  dbExecute(db,'VACUUM;')
  cat(paste('imported', n2, 'responses from', n1, 'persons\n'))
}


sql_data_type = function(value)
{
  if(inherits(value,'Date')) return(' DATE ')
  else if(inherits(value,'factor')) return(' TEXT ')
  else if(inherits(value,'POSIXlt') || inherits(value,'POSIXt')) return(' DATETIME ')
  else if(typeof(value) == 'integer') return(' INTEGER ')
  else if(is.numeric(value)) return(' DOUBLE PRECISION ')
  else return(" TEXT ")
}
