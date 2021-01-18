## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(warning=FALSE,
                      fig.height = 6,
                      fig.width = 8,
                      fig.retina=1,
                      fig.keep='high',
                      fig.align='center')

## -----------------------------------------------------------------------------
library(Radviz)

## -----------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

## -----------------------------------------------------------------------------
library(bodenmiller)
data(refPhenoMat)
data(refFuncMat)
data(refAnnots)
ref.df <- data.frame(refAnnots,
                     refPhenoMat,
                     refFuncMat)

## -----------------------------------------------------------------------------
trans <- function(coln) do.L(coln,fun=function(x) quantile(x,c(0.005,0.995)))

## -----------------------------------------------------------------------------
hist(ref.df$CD3)
abline(v=quantile(ref.df$CD3,c(0.005,0.995)),
       col=2,lty=2)

## -----------------------------------------------------------------------------
ct.S <- make.S(dimnames(refPhenoMat)[[2]])

## -----------------------------------------------------------------------------
## compute the similarity matrix
ct.sim <- cosine(as.matrix(ref.df[,row.names(ct.S)]))
## the current Radviz-independent measure of projection efficiency
in.da(ct.S,ct.sim)
## the current Radviz-independent measure of projection efficiency
rv.da(ct.S,ct.sim)

## -----------------------------------------------------------------------------
optim.ct <- do.optimRadviz(ct.S,ct.sim,iter=100,n=1000)
ct.S <- make.S(get.optim(optim.ct))

## ----echo=FALSE,results='asis'------------------------------------------------
ksink <- lapply(dimnames(refPhenoMat)[[2]],function(x) cat(' *',x,'\n'))

## ----echo=FALSE,results='asis'------------------------------------------------
ksink <- lapply(row.names(ct.S),function(x) cat(' *',x,'\n'))

## -----------------------------------------------------------------------------
ct.S <- recenter(ct.S,'CD3')

## ----echo=FALSE,results='asis'------------------------------------------------
ksink <- lapply(row.names(ct.S),function(x) cat(' *',x,'\n'))

## -----------------------------------------------------------------------------
ct.rv <- do.radviz(ref.df,ct.S,trans=trans)

## -----------------------------------------------------------------------------
summary(ct.rv)

## -----------------------------------------------------------------------------
head(ct.rv)

## -----------------------------------------------------------------------------
dim(ct.rv)

## -----------------------------------------------------------------------------
ct.rv

## -----------------------------------------------------------------------------
plot(ct.rv,anchors.only=FALSE)

## -----------------------------------------------------------------------------
plot(ct.rv)+
  geom_point()

## -----------------------------------------------------------------------------
plot(ct.rv)+
  geom_point(data=. %>% 
               arrange(CD4),
             aes(color=CD4))+
  scale_color_gradient(low='grey80',high="dodgerblue4")

## -----------------------------------------------------------------------------
smoothRadviz(ct.rv)

## -----------------------------------------------------------------------------
smoothRadviz(ct.rv)+
  geom_point(shape='.',alpha=1/5)

## -----------------------------------------------------------------------------
contour(ct.rv)

## -----------------------------------------------------------------------------
cur.pop <- 'igm+'
sub.rv <- subset(ct.rv,refAnnots$Cells==cur.pop)
smoothRadviz(ct.rv)+
  geom_density2d(data=sub.rv$proj$data,
                 aes(x=rx,y=ry),
                 color='black')

## -----------------------------------------------------------------------------
hexplot(ct.rv)

## -----------------------------------------------------------------------------
hexplot(ct.rv,color='CD4')

## -----------------------------------------------------------------------------
hexplot(ct.rv,color='pS6')

## -----------------------------------------------------------------------------
hexplot(ct.rv,color='pAkt')

## -----------------------------------------------------------------------------
hexplot(ct.rv,color='pErk')

## ----results='asis'-----------------------------------------------------------
ksink <- lapply(levels(refAnnots$Cells),function(x) cat(' *',x,'\n'))

## -----------------------------------------------------------------------------
bubbleRadviz(ct.rv,group = 'Cells')

## -----------------------------------------------------------------------------
bubbleRadviz(ct.rv,group = 'Cells',color='pS6')

## -----------------------------------------------------------------------------
data(untreatedPhenoMat)
data(untreatedFuncMat)
data(untreatedAnnots)
untreated.df <- bind_rows(ref.df %>% 
                            mutate(Treatment='unstimulated',
                                   Source=as.character(Source),
                                   Cells=as.character(Cells)),
                          data.frame(untreatedAnnots,
                                     untreatedPhenoMat,
                                     untreatedFuncMat) %>% 
                            mutate(Treatment=as.character(Treatment),
                                   Source=as.character(Source),
                                   Cells=as.character(Cells))) %>% 
  mutate(Treatment=factor(Treatment),
         Treatment=relevel(Treatment,'unstimulated'),
         Cells=factor(Cells))

## -----------------------------------------------------------------------------
tcells.df <- untreated.df %>% 
  filter(Cells %in% c('cd4+','cd8+'))
tcells.df %>% 
  count(Cells,Treatment)

## -----------------------------------------------------------------------------
func.S <- make.S(dimnames(refFuncMat)[[2]])
func.sim <- cosine(as.matrix(tcells.df[,row.names(func.S)]))
optim.func <- do.optimRadviz(func.S,func.sim,iter=100,n=1000)
func.S <- make.S(tail(optim.func$best,1)[[1]])
func.rv <- do.radviz(tcells.df,func.S,trans=trans)

## ----fig.width=12-------------------------------------------------------------
smoothRadviz(subset(func.rv,tcells.df$Treatment=='unstimulated'))+
  facet_grid(~Cells)

## ----fig.width=12-------------------------------------------------------------
plot(func.rv)+
  geom_density2d(aes(color=Treatment))+
  facet_grid(~Cells)

## -----------------------------------------------------------------------------
tcells.df %>% 
  select(Cells, Treatment, pS6, pSlp76) %>% 
  gather('Channel','value',one_of('pS6','pSlp76')) %>% 
  ggplot(aes(x=Treatment,y=value))+
  geom_boxplot(aes(fill=Treatment))+
  facet_grid(Channel~Cells)+
  theme_light()+
  theme(axis.text.x=element_text(angle=45,hjust=1))

