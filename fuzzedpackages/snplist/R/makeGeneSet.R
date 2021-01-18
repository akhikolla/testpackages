makeGeneSet <-function(annoInfo=NULL, margin=0,
                      annoTable='anno',geneTable='gene',allTable='allchrpos',
                      db='snplistdb',dbCleanUp=FALSE) {
    drv <- dbDriver("SQLite")
    dbFile <- sprintf("%s.sqlite",db)
    conn <- dbConnect(drv, dbname = dbFile)
    
    chrPosView <- sprintf("%sChrPos",annoTable)
    RSIDtoGeneTable <- sprintf("%sToGene",annoTable)

    #-----------------------------------------------------------------------
    tables <- dbGetQuery(conn,"SELECT name FROM sqlite_master WHERE type='table';")
    views <- dbGetQuery(conn,"SELECT name FROM sqlite_master WHERE type='view';")
    if(nrow(tables)!=0 && any(RSIDtoGeneTable==tables)) {
        cat("Remove the existing table",RSIDtoGeneTable, "\n")
        dbExecute(conn,sprintf("DROP TABLE %s;",RSIDtoGeneTable))
    }    
    if(nrow(views)!=0 && any(chrPosView==views)) {
        cat("Remove the existing view",chrPosView, "\n")
        dbExecute(conn,sprintf("DROP VIEW %s;",chrPosView))
    }    
    if(nrow(tables)!=0 && any(annoTable==tables)) {
        cat("Remove the existing table",annoTable, "\n")
        dbExecute(conn,sprintf("DROP TABLE %s;",annoTable))
    }
    #if(nrow(views)!=0 && any(annoTable==views)) {
    #    cat("Remove the existing table",annoTable, "\n")
    #    dbExecute(conn,sprintf("DROP VIEW %s;",annoTable))
    #}
    
    #-----------------------------------------------------------------------
    if(is.null(annoInfo)) {
        dbExecute(conn, 
            sprintf("CREATE TABLE %s AS SELECT rsid from %s",
            annoTable,allTable))
    }
    else if(suppressWarnings(isFile(annoInfo[1]))){ # Use of '[1]' prevents warning from multiple comparisons
        dbExecute(conn,sprintf("CREATE TABLE %s(rsid,primary key(rsid))",annoTable))
        field.types <- list(rsid="TEXT")
        dbWriteTable(conn,annoTable,annoInfo, 
                     row.names=FALSE,append=TRUE,header=FALSE,field.types=field.types)
    } 
    else {
        dbExecute(conn,sprintf("CREATE TABLE %s(rsid,primary key(rsid))",annoTable))
        if(class(annoInfo)!="data.frame") {
            annoInfo <- data.frame("rsid"=annoInfo)
        }
        else if ( !all('rsid' %in% names(annoInfo)) ) {
            stop("The column names should have 'rsid'.")
        }
        dbWriteTable(conn,annoTable,annoInfo,row.names=FALSE,append=TRUE)
    }
    
    cat("create table",annoTable,"\n")
    print( dbGetQuery(conn,sprintf("SELECT COUNT(*) FROM %s;",annoTable)) )
    print( dbGetQuery(conn,sprintf("SELECT * FROM %s LIMIT 5;",annoTable)) )
    cat(".....\n\n")

    #-----------------------------------------------------------------------
    dbExecute(conn, 
        sprintf("CREATE VIEW %s AS SELECT * from %s INNER JOIN %s USING(rsid)",
        chrPosView,allTable,annoTable))

    cat("Create view",chrPosView,"\n")
    #print( dbGetQuery(conn,sprintf("SELECT COUNT(*) FROM %s;",chrPosView)) )
    print( dbGetQuery(conn,sprintf("SELECT * FROM %s LIMIT 5;",chrPosView)) )
    cat(".....\n\n")

    #-----------------------------------------------------------------------    
    dbExecute(conn, 
        paste("CREATE TABLE ",RSIDtoGeneTable," AS 
                SELECT ", chrPosView,".rsid AS rsid, ", geneTable,".gene AS gene
                FROM   ", chrPosView,",", geneTable, "
                WHERE  ", chrPosView,".chr = ", geneTable,".chr 
                   AND ", chrPosView,".pos > ", geneTable,".start - ", margin," 
                   AND ", chrPosView,".pos < ", geneTable,".end   + ", margin,";
                ",sep = ""))

    cat("Create table",RSIDtoGeneTable,"\n")
    print( dbGetQuery(conn,sprintf("SELECT COUNT(*) FROM %s;",RSIDtoGeneTable)) )
    print( dbGetQuery(conn,sprintf("SELECT * FROM %s LIMIT 10;",RSIDtoGeneTable)) )
    cat(".....\n\n")

    #----------------------------------------------------------------------- 
    rsidToGene<-dbReadTable(conn,RSIDtoGeneTable)
    if(dbCleanUp) {
        cat("Remove the table",RSIDtoGeneTable,"\n")
        dbExecute(conn,sprintf("DROP TABLE %s;",RSIDtoGeneTable))
        cat("Remove the table",chrPosView,"\n")
        dbExecute(conn,sprintf("DROP VIEW %s;",chrPosView))
        cat("Remove the table",annoTable,"\n")
        dbExecute(conn,sprintf("DROP TABLE %s;",annoTable))
    }
    dbDisconnect(conn)

    n<-nrow(rsidToGene)
    geneset<-vector("list")
    for(i in 1:n) {
        gene<-rsidToGene[i,"gene"]
        rsid<-rsidToGene[i,"rsid"]
        geneset[[gene]]<-c(geneset[[gene]],rsid)
    }

    cat("Result :\n")
    cat("Gene sets :",length(geneset),"\n")
    cat("rsIDs :",length(unlist(geneset)),"\n")
    
    return(geneset)
}

