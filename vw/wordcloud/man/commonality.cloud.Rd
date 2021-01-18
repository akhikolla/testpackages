\name{commonality.cloud}
\alias{commonality.cloud}
\title{
  Plot a commonality cloud 
}
\description{
  Plot a cloud of words shared across documents
}
\usage{
commonality.cloud(term.matrix,comonality.measure=min,max.words=300,...)
}
\arguments{
  \item{term.matrix}{A term frequency matrix whose rows represent words and whose columns represent documents.}
  \item{comonality.measure}{A function taking a vector of frequencies for a single term, and returning a common frequency}
  \item{max.words}{Maximum number of words to be plotted. least frequent terms dropped}
  \item{...}{Additional parameters to be passed to wordcloud.}
}
\value{
nothing
 }
\examples{

if(require(tm)){
	data(SOTU)
	corp <- SOTU
	corp <- tm_map(corp, removePunctuation)
	corp <- tm_map(corp, content_transformer(tolower))
	corp <- tm_map(corp, removeNumbers)
	corp <- tm_map(corp, function(x)removeWords(x,stopwords()))

	term.matrix <- TermDocumentMatrix(corp)
	term.matrix <- as.matrix(term.matrix)
	colnames(term.matrix) <- c("SOTU 2010","SOTU 2011")
	comparison.cloud(term.matrix,max.words=40,random.order=FALSE)
	commonality.cloud(term.matrix,max.words=40,random.order=FALSE)
}

}