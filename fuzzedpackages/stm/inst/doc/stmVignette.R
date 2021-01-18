### R code from vignette source 'stmVignette.Rnw'

###################################################
### code chunk number 1: stmVignette.Rnw:73-74
###################################################
options(prompt = "R> ", continue = "+", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: stmVignette.Rnw:145-147
###################################################
library("stm")
set.seed(2138)


###################################################
### code chunk number 3: stmVignette.Rnw:156-157
###################################################
list(poliblog5k.docs[[1]][,1:5], poliblog5k.docs[[2]][,1:3])


###################################################
### code chunk number 4: stmVignette.Rnw:169-175 (eval = FALSE)
###################################################
## data <- read.csv("poliblogs2008.csv")
## processed <- textProcessor(data$documents, metadata = data)
## out <- prepDocuments(processed$documents, processed$vocab, processed$meta)
## docs <- out$documents
## vocab <- out$vocab
## meta  <-out$meta


###################################################
### code chunk number 5: stmVignette.Rnw:200-203 (eval = FALSE)
###################################################
## plotRemoved(processed$documents, lower.thresh = seq(1, 200, by = 100))
## out <- prepDocuments(processed$documents, processed$vocab,
##   processed$meta, lower.thresh = 15)


###################################################
### code chunk number 6: stmVignette.Rnw:228-232 (eval = FALSE)
###################################################
## poliblogPrevFit <- stm(documents = out$documents, vocab = out$vocab, 
##                        K = 20, prevalence =~ rating + s(day), 
##                        max.em.its = 75, data = out$meta, 
##                        init.type = "Spectral")


###################################################
### code chunk number 7: stmVignette.Rnw:237-238
###################################################
load(url("https://github.com/bstewart/stm/blob/gh-pages/files/VignetteObjects.RData?raw=true"))


###################################################
### code chunk number 8: stmVignette.Rnw:254-257 (eval = FALSE)
###################################################
## poliblogSelect <- selectModel(out$documents, out$vocab, K = 20,
##   prevalence =~ rating + s(day), max.em.its = 75,
##   data = out$meta, runs = 20, seed = 8458159)


###################################################
### code chunk number 9: stmVignette.Rnw:275-276
###################################################
plotModels(poliblogSelect, pch=c(1,2,3,4), legend.position="bottomright")


###################################################
### code chunk number 10: stmVignette.Rnw:285-286 (eval = FALSE)
###################################################
## selectedmodel <- poliblogSelect$runout[[3]] 


###################################################
### code chunk number 11: stmVignette.Rnw:296-298 (eval = FALSE)
###################################################
## storage <- searchK(out$documents, out$vocab, K = c(7, 10),
##   prevalence =~ rating + s(day), data = meta)


###################################################
### code chunk number 12: stmVignette.Rnw:321-322
###################################################
labelTopics(poliblogPrevFit, c(3, 7, 20))


###################################################
### code chunk number 13: stmVignette.Rnw:331-335
###################################################
thoughts3 <- findThoughts(poliblogPrevFit, texts = shortdoc,
  n = 2, topics = 3)$docs[[1]]
thoughts20 <- findThoughts(poliblogPrevFit, texts = shortdoc,
  n = 2, topics = 20)$docs[[1]]


###################################################
### code chunk number 14: stmVignette.Rnw:339-345
###################################################
thoughts3 <-
  c("Here's video of the ad we reported on below that the Obama campaign is running in Ohio responding to the earlier Swift-Boating spot tying Obama to former Weatherman Bill Ayers... With all our pr",
      "As noted here and elsewhere, the words 'William Ayers' appeared nowhere in yesterday's debate, despite the fact that the McCain campaign hinted for days that McCain would go hard at Obama's association")
thoughts20 <-
  c("Waxman calls for release of FBI interviews with Bush and Cheney.  In a letter to Attorney General Michael Mukasey today, Rep. Henry Waxman (D-CA), the Chairman of the Hous e Committee on Oversight",
      "Report: Bush 'Personally Directed' Gonzales To Strong-Arm Ashcroft At His BedsideIn his May 2007 testimony, describing the infamous strong-arming of John Ashcroft done by Andy Card and Alberto")


###################################################
### code chunk number 15: stmVignette.Rnw:351-354
###################################################
par(mfrow = c(1, 2),mar = c(.5, .5, 1, .5))
plotQuote(thoughts3, width = 30, main = "Topic 3")
plotQuote(thoughts20, width = 30, main = "Topic 20")


###################################################
### code chunk number 16: stmVignette.Rnw:372-376
###################################################
out$meta$rating <- as.factor(out$meta$rating)
prep <- estimateEffect(1:20 ~ rating + s(day), poliblogPrevFit,
  meta = out$meta, uncertainty = "Global")
summary(prep, topics=1)


###################################################
### code chunk number 17: stmVignette.Rnw:390-391
###################################################
plot(poliblogPrevFit, type = "summary", xlim = c(0, .3))


###################################################
### code chunk number 18: stmVignette.Rnw:413-420
###################################################
plot(prep, covariate = "rating", topics = c(3, 7, 20),
  model = poliblogPrevFit, method = "difference",
  cov.value1 = "Liberal", cov.value2 = "Conservative",
  xlab = "More Conservative ... More Liberal",
  main = "Effect of Liberal vs. Conservative",
  xlim = c(-.1, .1), labeltype = "custom",
  custom.labels = c('Obama', 'Sarah Palin','Bush Presidency'))


###################################################
### code chunk number 19: stmVignette.Rnw:435-442
###################################################
plot(prep, "day", method = "continuous", topics = 7,
  model = z, printlegend = FALSE, xaxt = "n", xlab = "Time (2008)")
monthseq <- seq(from = as.Date("2008-01-01"), 
  to = as.Date("2008-12-01"), by = "month")
monthnames <- months(monthseq)
axis(1,at = as.numeric(monthseq) - min(as.numeric(monthseq)), 
  labels = monthnames)


###################################################
### code chunk number 20: stmVignette.Rnw:454-457 (eval = FALSE)
###################################################
## poliblogContent <- stm(out$documents, out$vocab, K = 20,
##   prevalence =~ rating + s(day), content =~ rating,
##   max.em.its = 75, data = out$meta, init.type = "Spectral")


###################################################
### code chunk number 21: stmVignette.Rnw:465-466
###################################################
plot(poliblogContent, type = "perspectives", topics = 11)


###################################################
### code chunk number 22: stmVignette.Rnw:477-478
###################################################
plot(poliblogPrevFit, type = "perspectives", topics = c(12, 20))


###################################################
### code chunk number 23: stmVignette.Rnw:491-494 (eval = FALSE)
###################################################
## poliblogInteraction <- stm(out$documents, out$vocab, K = 20,
##   prevalence =~ rating * day, max.em.its = 75,
##   data = out$meta, init.type = "Spectral")


###################################################
### code chunk number 24: stmVignette.Rnw:501-515
###################################################
prep <- estimateEffect(c(20) ~ rating * day, poliblogInteraction,
  metadata = out$meta, uncertainty = "None")

plot(prep, covariate = "day", model = poliblogInteraction,
  method = "continuous", xlab = "Days", moderator = "rating",
  moderator.value = "Liberal", linecol = "blue", ylim = c(0, .12),
  printlegend = F)

plot(prep, covariate = "day", model = poliblogInteraction,
  method = "continuous", xlab = "Days", moderator = "rating",
  moderator.value = "Conservative", linecol = "red", add = T,
  printlegend = F)
legend(0, .08, c("Liberal", "Conservative"),
  lwd = 2, col = c("blue", "red"))


###################################################
### code chunk number 25: stmVignette.Rnw:532-533
###################################################
cloud(poliblogPrevFit, topic = 7, scale = c(2,.25))


###################################################
### code chunk number 26: stmVignette.Rnw:544-545
###################################################
mod.out.corr <- topicCorr(poliblogPrevFit)


###################################################
### code chunk number 27: stmVignette.Rnw:551-552
###################################################
plot(mod.out.corr)


###################################################
### code chunk number 28: stmVignette.Rnw:596-599
###################################################
plot(poliblogPrevFit$convergence$bound, type = "l",
  ylab = "Approximate Objective",
  main = "Convergence")


###################################################
### code chunk number 29: stmVignette.Rnw:651-653 (eval = FALSE)
###################################################
## set.seed(2138)
## heldout <- make.heldout(poliblog5k.docs, poliblog5k.voc)


###################################################
### code chunk number 30: stmVignette.Rnw:658-659 (eval = FALSE)
###################################################
## slam <- convertCorpus(heldout$documents, heldout$vocab, type="slam")


###################################################
### code chunk number 31: stmVignette.Rnw:663-671 (eval = FALSE)
###################################################
## mod1 <- CTM(slam, k = 100)
## control_CTM_VEM <- list(estimate.beta = TRUE, verbose = 1,
##   seed = as.integer(2138), nstart = 1L, best = TRUE,
##   var = list(iter.max = 20, tol = 1e-6),
##   em = list(iter.max = 1000, tol = 1e-3),
##   initialize = "random",
##   cg = list(iter.max = -1, tol = 1e-6))
## mod2 <- CTM(slam, k = 100, control = control_CTM_VEM)


###################################################
### code chunk number 32: stmVignette.Rnw:675-679 (eval = FALSE)
###################################################
## stm.mod1 <- stm(heldout$documents, heldout$vocab, K = 100,
##   init.type = "Spectral")
## stm.mod2 <- stm(heldout$documents, heldout$vocab, K = 100,
##   init.type = "Spectral", emtol = 1e-3)


