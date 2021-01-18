setSNPTable <-function(snpInfo,table='allchrpos',db='snplistdb') {
    drv   <- dbDriver("SQLite")
    dbFile<- sprintf("%s.sqlite",db)
    conn  <- dbConnect(drv, dbname = dbFile)

    tables <- dbGetQuery(conn,"SELECT name FROM sqlite_master WHERE type='table';")
    if(nrow(tables)!=0 && any(table==tables)) {
        cat("remove the existing table :",table, "\n\n")
        dbExecute(conn,sprintf("DROP TABLE %s;",table))
    }  
    
    field.types <- list(chr="TEXT",pos="INTEGER",rsid="TEXT")
    if(suppressWarnings(isFile(snpInfo[1]))){ # Use of '[1]' prevents warning from multiple comparisons
        dbWriteTable(conn,table,snpInfo,row.names=FALSE, 
                     header=FALSE,field.types=field.types,sep ="\t")
    } 
    else {
        if ( !all( names(field.types) %in% names(snpInfo)) ) {
            stop("The column names should include 'chr','pos','rsid'.")
        }
        snpInfo$chr <- as.character(snpInfo$chr)
        snpInfo$pos <- as.numeric(snpInfo$pos)
        dbWriteTable(conn,table,snpInfo,row.names=FALSE)
    }

    cat("Create Table :",table,"\n")
    print( dbGetQuery(conn,sprintf("SELECT * FROM %s LIMIT 10;",table)) )
    cat(".....\n\n")

    r <- dbGetQuery(conn,sprintf("SELECT COUNT(*) FROM %s;",table))
    dbDisconnect(conn)
    
    return(r)
}

