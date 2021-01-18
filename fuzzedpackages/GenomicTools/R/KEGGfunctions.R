# Functions based on the REST API from KEGG, description can be found here:
# http://www.kegg.jp/kegg/rest/keggapi.html

getKEGGOrganisms <- function(url="http://rest.kegg.jp/list/organism"){
  tmp <- readLines(url)
  tmp <- strsplit(tmp,"\t")
  out <- data.frame(KEGG=sapply(tmp,"[",1),
                    code=sapply(tmp,"[",2),
                    name=sapply(tmp,"[",3),
                    ontology=sapply(tmp,"[",4))
  out
}

getKEGGPathwayOverview <- function(code="hsa"){
  tmp <- readLines(paste("http://rest.kegg.jp/list/pathway/",code,sep=""))
  tmp <- strsplit(tmp,"\t")
  out <- data.frame(pathway=sapply(tmp,"[",1),
                    description=sapply(tmp,"[",2))
  out
}

getKEGGPathway <- function(pathway){
  tmp <- readLines(paste("http://rest.kegg.jp/get/",pathway,sep=""))
  ENTRY <- tmp[which(grepl("ENTRY",tmp)==TRUE)]
  NAME <- tmp[which(grepl("NAME",tmp)==TRUE)]
  DESCRIPTION <- tmp[which(grepl("DESCRIPTION",tmp)==TRUE)]
  CLASS <- tmp[which(grepl("CLASS",tmp)==TRUE)]
  PATHWAY_MAP <- tmp[which(grepl("PATHWAY_MAP",tmp)==TRUE)]
  ORGANISM <- tmp[which(grepl("ORGANISM",tmp)==TRUE)]
  GENE <- which(grepl("GENE",tmp)==TRUE)
  COMPOUND <- which(grepl("COMPOUND",tmp)==TRUE)
  leadingChr <- substr(tmp,0,1)
  # Now get the gene rows
   if(length(GENE)>0){
      geneStart <- GENE
      checkThis <- GENE + 1
      while(leadingChr[checkThis]==" "){
        checkThis <- checkThis + 1
      }
      geneEnd <- checkThis - 1
      GENE <- tmp[geneStart:geneEnd]
   # Now bring the data into some better form
      GENE[1] <- gsub("GENE", "    ",GENE[1])
      GENE <- strsplit(GENE,";")
      GENEtmp1 <- strsplit(trim.leading(sapply(GENE,"[",1)),"  ")
      GENEtmp2 <- trim.leading(sapply(GENE,"[",2))
      GENE <- data.frame(GeneID = sapply(GENEtmp1,"[",1),
                         GeneName = sapply(GENEtmp1,"[",2),
                         GeneDesc = GENEtmp2)
   } else {
      GENE <- data.frame(GeneID = "Not available",
                         GeneName = "Not available",
                         GeneDesc = "Not available")
   }     
  # Now get the compound rows
    if(length(COMPOUND)>0){
      compoundStart <- COMPOUND
      checkThis <- COMPOUND + 1
      while(leadingChr[checkThis]==" "){
        checkThis <- checkThis + 1
      }
      compoundEnd <- checkThis - 1
      COMPOUND <- tmp[compoundStart:compoundEnd]
    } else {
      COMPOUND <- "Not available"
    }
  out <- list(Entry = ENTRY,
              Name = NAME,
              Description = DESCRIPTION,
              Class = CLASS,
              Pathway_map= PATHWAY_MAP,
              Organism = ORGANISM,
              Gene = GENE,
              Compound = COMPOUND)
  out
}

getKEGGPathwayImage <- function(pathway, folder=NULL){
  pathway <- gsub("path:","", pathway)
  if(is.null(folder)) folder <- getwd()
  filename <- paste(pathway,".png",sep="")
  dlURL <- paste("http://rest.kegg.jp/get/",pathway,"/image",sep="")
  download.file(url=dlURL, destfile= file.path(folder,filename), mode="wb")
}