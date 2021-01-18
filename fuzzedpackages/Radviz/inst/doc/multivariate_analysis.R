## ----setup,echo=FALSE,include=FALSE-------------------------------------------
library(knitr)
library(bodenmiller)
library(ggplot2)
library(cytofan)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(scales)
library(igraph)
library(Radviz)

knitr::opts_chunk$set(warning=TRUE,
                      fig.keep='high',
                      fig.align='center',
                      fig.height = 5,
                      fig.width = 6)

## -----------------------------------------------------------------------------
data(refPhenoMat)
data(refFuncMat)
data(refAnnots)
ref.df <- data.frame(refAnnots,
                     refPhenoMat,
                     refFuncMat)

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

btcells.df <- untreated.df %>% 
  filter(Cells %in% c('cd8+','igm+')) %>% 
  mutate(Cells=droplevels(Cells)) %>% 
  group_by(Cells,Treatment) %>% 
  mutate(cellID=seq(length(Cells))) %>% 
  unite('cellID',one_of(c('Treatment','Cells','cellID')),sep = '_',remove = FALSE)

## -----------------------------------------------------------------------------
left_join(btcells.df %>% 
  count(Cells,Treatment),
  untreated.df %>% 
    count(Treatment,name = 'Total'),
  by=c('Treatment')) %>% 
  mutate(Fraction=round(100*n/Total,1))

## -----------------------------------------------------------------------------
btcells.df <- btcells.df %>% 
  group_by(Cells,Treatment) %>% 
  sample_n(500,replace = TRUE)

## ----fig.width=8,fig.height=8-------------------------------------------------
btcells.df %>% 
  gather('Channel','value',
         one_of(colnames(refPhenoMat),colnames(refFuncMat))) %>% 
  filter(Channel %in% colnames(refPhenoMat)) %>% 
  ggplot(aes(x=Channel,y=value))+
  geom_fan()+
  facet_grid(Treatment~Cells)+
  theme_light(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ----fig.width=8,fig.height=8-------------------------------------------------
btcells.df %>% 
  gather('Channel','value',
         one_of(colnames(refPhenoMat),colnames(refFuncMat))) %>% 
  filter(Channel %in% colnames(refFuncMat)) %>% 
  ggplot(aes(x=Channel,y=value))+
  geom_fan()+
  facet_grid(Treatment~Cells)+
  theme_light(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## -----------------------------------------------------------------------------
btcells.channels <- btcells.df %>% 
  gather('Channel','value',
         one_of(colnames(refPhenoMat),colnames(refFuncMat))) %>% 
  group_by(Channel,Treatment,Cells) %>% 
  summarize(value=quantile(value,0.8)) %>% 
  group_by(Channel) %>% 
  summarise(value=max(value)) %>% 
  filter(value>1) %>% 
  .$Channel

## ----results='asis'-----------------------------------------------------------
ksink <- lapply(btcells.channels,function(x) cat(' -',x,'\n'))

## -----------------------------------------------------------------------------
sim.mat <- cosine(as.matrix(btcells.df[,btcells.channels]))
classic.S <- make.S(btcells.channels)
classic.optim <- do.optimRadviz(classic.S,sim.mat)
classic.S <- make.S(get.optim(classic.optim))
btcells.rv <- do.radviz(btcells.df,classic.S)

## -----------------------------------------------------------------------------
plot(btcells.rv)+
  geom_point(aes(color=Treatment))

## -----------------------------------------------------------------------------
btcells.rv <- rescalePlot(btcells.rv)
plot(btcells.rv)+
  geom_point(aes(color=Treatment))

## -----------------------------------------------------------------------------
plot(btcells.rv)+
  geom_point(aes(color=Cells))

## ----fig.width=9,fig.height=4-------------------------------------------------
plot(btcells.rv)+
  geom_point(aes(color=Treatment))+
  facet_wrap(~Cells)+
  theme_radviz(base_size = 16)

## -----------------------------------------------------------------------------
btcells.df <- btcells.df %>% 
  unite('Condition',c('Treatment','Cells'),sep='_',remove = FALSE)
treat.S <- do.optimFreeviz(btcells.df[,btcells.channels],
                           classes = btcells.df$Condition)
btcells.fv <- do.radviz(btcells.df, treat.S)

## -----------------------------------------------------------------------------
plot(btcells.fv)+
  geom_point(aes(color=Treatment))

## ----fig.width=9,fig.height=4-------------------------------------------------
plot(btcells.fv)+
  geom_point(aes(color=Treatment))+
  facet_wrap(~Cells)+
  theme_radviz(base_size = 16)

## -----------------------------------------------------------------------------
btcells.fv <- rescalePlot(btcells.fv,fraction=0.5)

## -----------------------------------------------------------------------------
plot(btcells.fv)+
  geom_point(aes(color=Treatment))

## -----------------------------------------------------------------------------
contour(btcells.fv,
        color='Treatment')

## -----------------------------------------------------------------------------
plot(btcells.fv,
     anchors.filter = 0.5)+
  geom_point(aes(color=Treatment))

## -----------------------------------------------------------------------------
btcells.dist <- as.matrix(btcells.df[,btcells.channels])
rownames(btcells.dist) <- btcells.df$cellID
btcells.dist <- btcells.dist%*% t(btcells.dist)
btcells.dist <- btcells.dist/(sqrt(diag(btcells.dist)) %*% t(sqrt(diag(btcells.dist))))
btcells.dist[btcells.dist>1] <- 1
btcells.dist[btcells.dist<0] <- 0
btcells.dist <- 2*acos(btcells.dist)/pi

diag(btcells.dist) <- NA # avoid self loops

## -----------------------------------------------------------------------------
K <- floor(nrow(btcells.df)^0.5)
btcells.adj <- apply(btcells.dist,1,rank,na.last = TRUE)
btcells.adj[btcells.adj<=K] <- btcells.dist[btcells.adj<=K]
btcells.adj[btcells.adj>K] <- 0

## -----------------------------------------------------------------------------
bind_rows(data.frame(value=btcells.dist[sample(1000)],
                     Type='Overall',
                     stringsAsFactors = FALSE),
          data.frame(value=btcells.adj[btcells.adj!=0][sample(1000)],
                     Type='Nearest Neighbors',
                     stringsAsFactors = FALSE)) %>% 
  mutate(Type=factor(Type),
         Type=relevel(Type,'Overall')) %>% 
  filter(!is.na(value)) %>% 
  ggplot(aes(x=value))+
  geom_histogram(aes(fill=Type),
                 position = 'identity',
                 bins=50,
                 alpha=0.5)+
  theme_light(base_size=16)

## -----------------------------------------------------------------------------
btcells.weights <- btcells.adj
btcells.weights <- exp(-btcells.weights^2/(2*median(btcells.weights[btcells.weights!=0])^2))
btcells.weights[btcells.adj==0] <- 0

## -----------------------------------------------------------------------------
btcells.graph <- graph_from_adjacency_matrix(btcells.weights,
                                            mode='undirected',
                                            weighted = TRUE,
                                            diag = FALSE)

## -----------------------------------------------------------------------------
btcells.groups <- cluster_louvain(btcells.graph)
btcells.df <- btcells.df %>% 
  ungroup() %>% 
  mutate(Group=membership(btcells.groups),
         Group=as.numeric(Group),
         Group=factor(Group))

## ----fig.height=12, fig.width=9-----------------------------------------------
btcells.df %>% 
  gather('Channel','value',
         one_of(colnames(refPhenoMat),colnames(refFuncMat))) %>% 
  filter(Channel %in% btcells.channels) %>% 
  ggplot(aes(x=Channel,y=value))+
  geom_fan()+
  facet_grid(Group~.)+
  theme_light(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## -----------------------------------------------------------------------------
btcells.df %>% 
  count(Group,Cells) %>% 
  group_by(Group) %>% 
  mutate(n=n/sum(n)) %>% 
  spread(Cells, n)

## -----------------------------------------------------------------------------
btcells.df %>% 
  count(Group,Treatment) %>% 
  group_by(Group) %>% 
  mutate(n=n/sum(n)) %>% 
  spread(Treatment, n)

## -----------------------------------------------------------------------------
btcells.S <- do.optimGraphviz(btcells.df[,btcells.channels],btcells.graph)
btcells.gv <- do.radviz(btcells.df, btcells.S,
                       graph = btcells.graph)

## -----------------------------------------------------------------------------
plot(btcells.gv)+
  geom_point(aes(color=Group))

## -----------------------------------------------------------------------------
btcells.graph$layout <- layout_with_drl(btcells.graph)
community.cols <- hue_pal()(length(btcells.groups))
plot(btcells.graph,
     vertex.label=NA,
     vertex.color=community.cols[membership(btcells.groups)])

## -----------------------------------------------------------------------------
plot(btcells.gv)+
  geom_point(aes(color=Cells),alpha=0.5)

## -----------------------------------------------------------------------------
plot(btcells.gv)+
  geom_point(aes(color=Treatment))

