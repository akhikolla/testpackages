## ---- message=FALSE, warning=FALSE, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# knitr options
knitr::opts_chunk$set(fig.width=6, fig.height=4.5)
options(width=800)

## ---- message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# libs
library(bayes4psy)
library(cowplot)
library(dplyr)
library(ggplot2)

# load data
data_all <- after_images

# load stimuli
stimuli <- after_images_stimuli

## ---- message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# prepare data
data_red <- data_all %>% filter(stimuli == "red")
data_red <- data.frame(r=data_red$r,
                       g=data_red$g,
                       b=data_red$b)

# fit
fit_red <- b_color(colors=data_red, chains=1, iter=200, warmup=100)

# inspect
plot_trace(fit_red)
plot_hsv(fit_red)
# the command below is commented out for the sake of brevity
#print(fit_red)

## ---- message=FALSE, warning=FALSE, results = 'hide'--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# green
data_green <- data_all %>% filter(stimuli == "green")
data_green <- data.frame(r=data_green$r,
                         g=data_green$g,
                         b=data_green$b)
fit_green <- b_color(colors=data_green, chains=1, iter=200, warmup=100)

# blue
data_blue <- data_all %>% filter(stimuli == "blue")
data_blue <- data.frame(r=data_blue$r,
                        g=data_blue$g,
                        b=data_blue$b)
fit_blue <- b_color(colors=data_blue, chains=1, iter=200, warmup=100)

# yellow
data_yellow <- data_all %>% filter(stimuli == "yellow")
data_yellow <- data.frame(r=data_yellow$r,
                          g=data_yellow$g,
                          b=data_yellow$b)
fit_yellow <- b_color(colors=data_yellow, chains=1, iter=200, warmup=100)

# cyan
data_cyan <- data_all %>% filter(stimuli == "cyan")
data_cyan <- data.frame(r=data_cyan$r,
                        g=data_cyan$g,
                        b=data_cyan$b)
fit_cyan <- b_color(colors=data_cyan, chains=1, iter=200, warmup=100)

# magenta
data_magenta <- data_all %>% filter(stimuli == "magenta")
data_magenta <- data.frame(r=data_magenta$r,
                           g=data_magenta$g,
                           b=data_magenta$b)
fit_magenta <- b_color(colors=data_magenta, chains=1, iter=200, warmup=100)

## ---- message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# load theory predictions
trichromatic <- after_images_trichromatic
opponent_process <- after_images_opponent_process

# red
stimulus <- "red"
lines <- list()
lines[[1]] <- c(trichromatic[trichromatic$stimuli == stimulus, ]$h,
                trichromatic[trichromatic$stimuli == stimulus, ]$s,
                trichromatic[trichromatic$stimuli == stimulus, ]$v)
lines[[2]] <- c(opponent_process[opponent_process$stimuli == stimulus, ]$h,
                opponent_process[opponent_process$stimuli == stimulus, ]$s,
                opponent_process[opponent_process$stimuli == stimulus, ]$v)
    
points <- list()
points[[1]] <- c(stimuli[stimuli$stimuli == stimulus, ]$h_s,
                 stimuli[stimuli$stimuli == stimulus, ]$s_s,
                 stimuli[stimuli$stimuli == stimulus, ]$v_s)
    
plot_red <- plot_distributions_hsv(fit_red, points=points,
                                   lines=lines, hsv=TRUE)

plot_red <- plot_red + ggtitle("Red") +
  theme(plot.title = element_text(hjust = 0.5))

# green
stimulus <- "green"
lines <- list()
lines[[1]] <- c(trichromatic[trichromatic$stimuli == stimulus, ]$h,
                trichromatic[trichromatic$stimuli == stimulus, ]$s,
                trichromatic[trichromatic$stimuli == stimulus, ]$v)
lines[[2]] <- c(opponent_process[opponent_process$stimuli == stimulus, ]$h,
                opponent_process[opponent_process$stimuli == stimulus, ]$s,
                opponent_process[opponent_process$stimuli == stimulus, ]$v)

points <- list()
points[[1]] <- c(stimuli[stimuli$stimuli == stimulus, ]$h_s,
                 stimuli[stimuli$stimuli == stimulus, ]$s_s,
                 stimuli[stimuli$stimuli == stimulus, ]$v_s)

plot_green <- plot_distributions_hsv(fit_green, points=points,
                                     lines=lines, hsv=TRUE)
plot_green <- plot_green + ggtitle("Green") +
  theme(plot.title = element_text(hjust = 0.5))

# blue
stimulus <- "blue"
lines <- list()
lines[[1]] <- c(trichromatic[trichromatic$stimuli == stimulus, ]$h,
                trichromatic[trichromatic$stimuli == stimulus, ]$s,
                trichromatic[trichromatic$stimuli == stimulus, ]$v)
lines[[2]] <- c(opponent_process[opponent_process$stimuli == stimulus, ]$h,
                opponent_process[opponent_process$stimuli == stimulus, ]$s,
                opponent_process[opponent_process$stimuli == stimulus, ]$v)

points <- list()
points[[1]] <- c(stimuli[stimuli$stimuli == stimulus, ]$h_s,
                 stimuli[stimuli$stimuli == stimulus, ]$s_s,
                 stimuli[stimuli$stimuli == stimulus, ]$v_s)

plot_blue <- plot_distributions_hsv(fit_blue, points=points,
                                    lines=lines, hsv=TRUE)
plot_blue <- plot_blue + ggtitle("Blue") +
  theme(plot.title = element_text(hjust = 0.5))

# yellow
stimulus <- "yellow"
lines <- list()
lines[[1]] <- c(trichromatic[trichromatic$stimuli == stimulus, ]$h,
                trichromatic[trichromatic$stimuli == stimulus, ]$s,
                trichromatic[trichromatic$stimuli == stimulus, ]$v)
lines[[2]] <- c(opponent_process[opponent_process$stimuli == stimulus, ]$h,
                opponent_process[opponent_process$stimuli == stimulus, ]$s,
                opponent_process[opponent_process$stimuli == stimulus, ]$v)

points <- list()
points[[1]] <- c(stimuli[stimuli$stimuli == stimulus, ]$h_s,
                 stimuli[stimuli$stimuli == stimulus, ]$s_s,
                 stimuli[stimuli$stimuli == stimulus, ]$v_s)

plot_yellow <- plot_distributions_hsv(fit_yellow, points=points,
                                      lines=lines, hsv=TRUE)
plot_yellow <- plot_yellow + ggtitle("Yellow") +
  theme(plot.title = element_text(hjust = 0.5))


# cyan
stimulus <- "cyan"
lines <- list()
lines[[1]] <- c(trichromatic[trichromatic$stimuli == stimulus, ]$h,
                trichromatic[trichromatic$stimuli == stimulus, ]$s,
                trichromatic[trichromatic$stimuli == stimulus, ]$v)
lines[[2]] <- c(opponent_process[opponent_process$stimuli == stimulus, ]$h,
                opponent_process[opponent_process$stimuli == stimulus, ]$s,
                opponent_process[opponent_process$stimuli == stimulus, ]$v)

points <- list()
points[[1]] <- c(stimuli[stimuli$stimuli == stimulus, ]$h_s,
                 stimuli[stimuli$stimuli == stimulus, ]$s_s,
                 stimuli[stimuli$stimuli == stimulus, ]$v_s)

plot_cyan <- plot_distributions_hsv(fit_cyan, points=points,
                                    lines=lines, hsv=TRUE)
plot_cyan <- plot_cyan + ggtitle("Cyan") +
  theme(plot.title = element_text(hjust = 0.5))


# magenta
stimulus <- "magenta"
lines <- list()
lines[[1]] <- c(trichromatic[trichromatic$stimuli == stimulus, ]$h,
                trichromatic[trichromatic$stimuli == stimulus, ]$s,
                trichromatic[trichromatic$stimuli == stimulus, ]$v)
lines[[2]] <- c(opponent_process[opponent_process$stimuli == stimulus, ]$h,
                opponent_process[opponent_process$stimuli == stimulus, ]$s,
                opponent_process[opponent_process$stimuli == stimulus, ]$v)

points <- list()
points[[1]] <- c(stimuli[stimuli$stimuli == stimulus, ]$h_s,
                 stimuli[stimuli$stimuli == stimulus, ]$s_s,
                 stimuli[stimuli$stimuli == stimulus, ]$v_s)

plot_magenta <- plot_distributions_hsv(fit_magenta, points=points,
                                       lines=lines, hsv=TRUE)
plot_magenta <- plot_magenta + ggtitle("Magenta") +
  theme(plot.title = element_text(hjust = 0.5))

## ---- message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_grid(plot_red, plot_green, plot_blue,
          plot_yellow, plot_cyan, plot_magenta,
          ncol=3, nrow=2, scale=0.9)

