% 
% TODO: Write Documentation 
% 



\name{wordcloud}
\alias{wordcloud}
\title{
  Plot a word cloud
}
\description{
  Plot a word cloud
}
\usage{
wordcloud(words,freq,scale=c(4,.5),min.freq=3,max.words=Inf,
	random.order=TRUE, random.color=FALSE, rot.per=.1,
	colors="black",ordered.colors=FALSE,use.r.layout=FALSE,
	fixed.asp=TRUE, ...)
}
\arguments{
  \item{words}{the words}
  \item{freq}{their frequencies}
  \item{scale}{A vector of length 2 indicating the range of the size of the words. }
  \item{min.freq}{words with frequency below min.freq will not be plotted}
  \item{max.words}{Maximum number of words to be plotted. least frequent terms dropped}
  \item{random.order}{plot words in random order. If false, they will be plotted in decreasing frequency}
  \item{random.color}{choose colors randomly from the colors. If false, the color is chosen based on the frequency}
  \item{rot.per}{proportion words with 90 degree rotation}
  \item{colors}{color words from least to most frequent}
  \item{ordered.colors}{if true, then colors are assigned to words in order}
  \item{use.r.layout}{if false, then c++ code is used for collision detection, otherwise R is used}
  \item{fixed.asp}{if TRUE, the aspect ratio is fixed. Variable aspect ratio only supported if rot.per==0}
  \item{...}{Additional parameters to be passed to text (and strheight,strwidth).}
}
\details{
If freq is missing, then words can either be a character vector, or Corpus. If it is a
vector and freq is missing, standard stop words will be removed prior to plotting.
}
\value{
nothing
}
\seealso{
  \code{\link{text}}
}
\examples{

wordcloud(c(letters, LETTERS, 0:9), seq(1, 1000, len = 62))

if(require(tm)){

	##### 			from character 		#####
	wordcloud(
"Many years ago the great British explorer George Mallory, who 
was to die on Mount Everest, was asked why did he want to climb 
it. He said, \"Because it is there.\"

Well, space is there, and we're going to climb it, and the 
moon and the planets are there, and new hopes for knowledge 
and peace are there. And, therefore, as we set sail we ask 
God's blessing on the most hazardous and dangerous and greatest 
adventure on which man has ever embarked.",
	,random.order=FALSE)

\dontrun{
	data(crude)
	crude <- tm_map(crude, removePunctuation)
	crude <- tm_map(crude, function(x)removeWords(x,stopwords()))
	
	##### 			from corpus 		#####
	wordcloud(crude)
	
	
	##### 		from frequency counts 	#####
	tdm <- TermDocumentMatrix(crude)
	m <- as.matrix(tdm)
	v <- sort(rowSums(m),decreasing=TRUE)
	d <- data.frame(word = names(v),freq=v)
	
	wordcloud(d$word,d$freq)

	#A bigger cloud with a minimum frequency of 2
	wordcloud(d$word,d$freq,c(8,.3),2)

	#Now lets try it with frequent words plotted first
	wordcloud(d$word,d$freq,c(8,.5),2,,FALSE,.1)

	##### 			with colors 		#####
	if(require(RColorBrewer)){

		pal <- brewer.pal(9,"BuGn")
		pal <- pal[-(1:4)]
		wordcloud(d$word,d$freq,c(8,.3),2,,FALSE,,.15,pal)


		pal <- brewer.pal(6,"Dark2")
		pal <- pal[-(1)]
		wordcloud(d$word,d$freq,c(8,.3),2,,TRUE,,.15,pal)
		
		#random colors
		wordcloud(d$word,d$freq,c(8,.3),2,,TRUE,TRUE,.15,pal)
	}
	##### 			with font 			#####

	wordcloud(d$word,d$freq,c(8,.3),2,,TRUE,,.15,pal,
		vfont=c("gothic english","plain"))

	wordcloud(d$word,d$freq,c(8,.3),2,100,TRUE,,.15,pal,vfont=c("script","plain"))
	
	wordcloud(d$word,d$freq,c(8,.3),2,100,TRUE,,.15,pal,vfont=c("serif","plain"))
}
}


} 
