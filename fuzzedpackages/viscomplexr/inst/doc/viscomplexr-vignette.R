## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

# This option anti-aliases the plots made below under Windows
if(Sys.info()[["sysname"]] == "Windows") {
  knitr::opts_chunk$set(dev = "CairoPNG")
}

options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup, echo = FALSE------------------------------------------------------
library(viscomplexr)

## ---- figure_1, fig.width = 5, fig.height = 5, results = 'hide', fig.align='center', cache = FALSE, fig.show = 'hold', fig.cap = 'Phase portrait of the function $f(z)=z$ in the window $\\left|\\Re(z)\\right| < 8.5$ and $\\left|\\Im(z)\\right| < 8.5$.'----
phasePortrait("z", xlim = c(-8.5, 8.5), ylim = c(-8.5, 8.5),
              xlab = "real", ylab = "imaginary", main = "f(z) = z",
              nCores = 2) # Probably not required on your machine (see below)

# Note the argument 'nCores' which determines the number of parallel processes to
# be used. Setting nCores = 2 has been done here and in all subsequent 
# examples as CRAN checks do not allow more parallel processes. 
# For normal work, we recommend not to define nCores at all which will make 
# phasePortrait use all available cores on your machine.

# The progress messages phasePortrait is writing to the console can be 
# suppressed by including 'verbose = FALSE' in the call (see documentation).

## ----figure_2, fig.width=5, fig.height=5, results="hide", fig.align='center', fig.show='hold', cache=TRUE, fig.cap= "Different options for including reference lines with the argument `pType`."----
# divide graphics device into four regions and adjust plot margins 
op <- par(mfrow = c(2, 2), 
          mar   = c(0.25, 0.55, 1.10, 0.25)) 
# plot four phase portraits with different choices of pType
phasePortrait("z", xlim = c(-8.5, 8.5), ylim = c(-8.5, 8.5), pType = "p",
              main = "pType = 'p'",   axes = FALSE, nCores = 2)
phasePortrait("z", xlim = c(-8.5, 8.5), ylim = c(-8.5, 8.5), pType = "pa",
              main = "pType = 'pa'",  axes = FALSE, nCores = 2)
phasePortrait("z", xlim = c(-8.5, 8.5), ylim = c(-8.5, 8.5), pType = "pm",
              main = "pType = 'pm'",  axes = FALSE, nCores = 2)
phasePortrait("z", xlim = c(-8.5, 8.5), ylim = c(-8.5, 8.5), pType = "pma",
              main = "pType = 'pma'", axes = FALSE, nCores = 2)
par(op) # reset the graphics parameters to their previous values

## ----eval=FALSE, figure_3, fig.width=5, fig.height=5, results='hide', fig.align='center', cache=TRUE, fig.show='hold', fig.cap='Phase portrait of the function $f(z)=\\frac{(3+2\\mathrm{i}+z)(-5+5\\mathrm{i}+z)}{(-2-2\\mathrm{i}+z)^2}$ in the window $\\left|\\Re(z)\\right| < 8.5$ and $\\left|\\Im(z)\\right| < 8.5$.'----
#  op <- par(mar = c(5.1, 4.1, 2.1, 2.1), cex = 0.8) # adjust plot margins
#                                                    # and general text size
#  phasePortrait("(3+2i+z)*(-5+5i+z)/(-2-2i+z)^2",
#                xlim = c(-8.5, 8.5), ylim = c(-8.5, 8.5),
#                xlab = "real", ylab = "imaginary",
#                nCores = 2) # Increase or leave out for higher performance
#  par(op) # reset the graphics parameters to their previous values

## ----eval = FALSE, figure_4, fig.width=7, fig.height=2.8, results='hide', fig.align='center', , fig.show='hold', cache=TRUE, fig.cap='The function $f(z)=\\frac{(3+2\\mathrm{i}+z)(-5+5\\mathrm{i}+z)}{(-2-2\\mathrm{i}+z)^2}$ portrayed with three different settings of `pi2Div` and `pType = "pa"`.'----
#  # divide graphics device into three regions and adjust plot margins
#  op <- par(mfrow = c(1, 3), mar = c(0.2, 0.2, 0.4, 0.2))
#  for(n in c(6, 9, 18)) {
#    phasePortrait("(3+2i+z)*(-5+5i+z)/(-2-2i+z)^2", xlim = c(-8.5, 8.5), ylim = c(-8.5, 8.5),
#                  pi2Div = n, pType = "pa", axes = FALSE, nCores = 2)
#    # separate title call (R base graphics) for nicer line adjustment, just cosmetics
#    title(paste("pi2Div =", n), line = -1.2)
#  }
#  par(op) # reset graphics parameters to previous values

## ----figure_5, fig.width=7, fig.height=2.8, results='hide', fig.align='center', , fig.show='hold', cache=TRUE, fig.cap='The function $f(z)=\\frac{(3+2\\mathrm{i}+z)(-5+5\\mathrm{i}+z)}{(-2-2\\mathrm{i}+z)^2}$ portrayed with three different settings of `pi2Div` and `pType = "pma"`.'----
# divide graphics device into three regions and adjust plot margins 
op <- par(mfrow = c(1, 3), mar = c(0.2, 0.2, 0.4, 0.2))
for(n in c(6, 9, 18)) {
  phasePortrait("(3+2i+z)*(-5+5i+z)/(-2-2i+z)^2", xlim = c(-8.5, 8.5), ylim = c(-8.5, 8.5),
                pi2Div = n, pType = "pma", axes = FALSE, nCores = 2)
  # separate title call (R base graphics) for nicer line adjustment, just cosmetics
  title(paste("pi2Div =", n), line = -1.2) 
}
par(op) # reset graphics parameters to previous values

## ----eval=FALSE, figure_6, fig.width=7, fig.height=2.8, results='hide', fig.align='center', , fig.show='hold', cache=TRUE, fig.cap='The function $f(z)=\\frac{(3+2\\mathrm{i}+z)(-5+5\\mathrm{i}+z)}{(-2-2\\mathrm{i}+z)^2}$ portrayed with decoupled settings of `pi2Div` and `logBase`.'----
#  # divide graphics device into three regions and adjust plot margins
#  op <- par(mfrow = c(1, 3), mar = c(0.2, 0.2, 0.4, 0.2))
#  for(n in c(6, 9, 18)) {
#    phasePortrait("(3+2i+z)*(-5+5i+z)/(-2-2i+z)^2", xlim = c(-8.5, 8.5), ylim = c(-8.5, 8.5),
#                  pi2Div = n, logBase = sqrt(3), pType = "pma", axes = FALSE, nCores = 2)
#    # separate title call (R base graphics) for nicer line adjustment, just cosmetics
#    title(paste("pi2Div = ", n, ", logBase = 3^(1/3)", sep = ""), line = -1.2)
#  }
#  par(op) # reset graphics parameters to previous values

## ----eval=FALSE, figure_7, fig.width = 5, fig.height = 5, results = 'hide', fig.align='center', fig.show='hold', cache = TRUE, fig.cap = 'Phase portrait of the function $f(z)=\\mathrm{e}^z$ in the window $\\left|\\Re(z)\\right| < 8.5$ and $\\left|\\Im(z)\\right| < 8.5$ with iso-modulus lines.'----
#  op <- par(mar = c(5.1, 4.1, 2.1, 2.1), cex = 0.8) # adjust plot margins
#                                                    # and general text size
#  phasePortrait(exp, xlim = c(-8.5, 8.5), ylim = c(-8.5, 8.5), pType = "pm",
#                xlab = "real", ylab = "imaginary", nCores = 2)
#  par(op) # reset graphics parameters to previous values

## ----eval=FALSE, figure_8, fig.width=7, fig.height=2.8, results='hide', fig.align='center', fig.show='hold', cache=TRUE, fig.cap='The function $f(z)=\\mathrm{e}^z$ portrayed with the default coupling of `pi2Div` and `logBase` as implemented in *phasePortrait*.'----
#  # divide graphics device into three regions and adjust plot margins
#  op <- par(mfrow = c(1, 3), mar = c(0.2, 0.2, 0.4, 0.2))
#  for(n in c(6, 9, 18)) {
#    phasePortrait("exp(z)", xlim = c(-8.5, 8.5), ylim = c(-8.5, 8.5),
#                  pi2Div = n, pType = "pma", axes = FALSE, nCores = 2)
#    # separate title call (R base graphics) for nicer line adjustment, just cosmetics
#    title(paste("pi2Div = ", n, ", logBase = exp(2*pi/pi2Div)", sep = ""),
#          line = -1.2, cex.main = 0.9)
#  }
#  par(op) # reset graphics parameters to previous values

## ----eval=FALSE, figure_9, fig.width=5, fig.height=5, fig.show='hold', results='hide', cache=TRUE, fig.cap='Tuning reference zone contrast with the parameters `darkestShade` (column-wise, 0, 0.2, 0.4), and `lambda` (row-wise, 0.1, 1, 10).'----
#  op <- par(mfrow = c(3, 3), mar = c(0.2, 0.2, 0.2, 0.2))
#  for(lb in c(0.1, 1, 10)) {
#    for(dS in c(0, 0.2, 0.4)) {
#      phasePortrait("tan(z^2)", xlim = c(-1.7, 1.7), ylim = c(-1.7, 1.7),
#                    pType = "pm", darkestShade = dS, lambda = lb,
#                    axes = FALSE, xaxs = "i", yaxs = "i", nCores = 2)
#    }
#  }
#  par(op)

## ----figure_10, fig.width=7, fig.height=2.7, results='hide', fig.align='center', fig.show='hold', cache=TRUE, fig.cap='Three phase portraits with branch cuts (dashed line), illustrating the three values of $f(z)=z^{1/3}$, $z \\in \\mathbb{C} \\setminus \\lbrace 0 \\rbrace$. The transitions between the phase portraits are indicated by same-coloured arrows pointing at the branch cuts.'----
op <- par(mfrow = c(1, 3), mar = c(0.4, 0.2, 0.2, 0.2))
for(k in 0:2) {
  FUNstring <- paste0("z^(1/3) * exp(1i * 2*pi/3 * ", k, ")")
  phasePortrait(FUN = FUNstring, 
                xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), pi2Div = 12, 
                axes = FALSE, nCores = 2)
  title(sub = paste0("k = ", k), line = -1)
  # emphasize branch cut with a dashed line segment
  segments(-1.5, 0, 0, 0, lwd = 2, lty = "dashed") 
  # draw colored arrows
  upperCol <- switch(as.character(k),
                     "0" = "black", "1" = "red", "2" = "green")
  lowerCol <- switch(as.character(k),
                     "0" = "green", "1" = "black", "2" = "red")
  arrows(x0 = c(-1.2), y0 = c(1, -1), y1 = c(0.2, -0.2), 
         lwd = 2, length = 0.1, col = c(upperCol, lowerCol))
}
par(op)

## ----figure_11, fig.width=7, fig.height=2.7, fig.align='center', results='hide', fig.show='hold', cache=TRUE, fig.cap='Three branches of $\\log z=\\log r+\\mathrm{i}\\cdot(\\varphi + k\\cdot2\\pi), r>0, \\varphi\\in\\left[0,2\\pi\\right[$, with $k=-1,0,1$. The branch cuts are marked with dashed white lines.'----
op <- par(mfrow = c(1, 3), mar = c(0.4, 0.2, 0.2, 0.2))
for(k in -1:1) {
  FUNstring <- paste0("log(Mod(z)) + 1i * (Arg(z) + 2 * pi * ", k, ")") 
  phasePortrait(FUN = FUNstring, pi2Div = 36,
                xlim = c(-2, 2), ylim = c(-2, 2), axes = FALSE, nCores = 2)
  segments(-2, 0, 0, 0, col = "white", lwd = 1, lty = "dashed")
  title(sub = paste0("k = ", k), line = -1)
}  
par(op)

## ----figure_12, fig.width=7, fig.height=3.5, fig.align='center', results='hide', fig.show='hold', cache=TRUE, fig.cap='Mapping the complex number plane on the Riemann sphere. Left: lower (southern) hemisphere; right upper (northern hemisphere). Folding both figures face to face along a vertical line in the middle between them can be imagined as closing the Riemann sphere.'----
op <- par(mfrow = c(1, 2), mar = rep(0.1, 4))
# Southern hemisphere
phasePortrait("z", xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4), 
              pi2Div = 12, axes = FALSE, nCores = 2)
riemannMask(annotSouth = TRUE)
# Northern hemisphere
phasePortrait("z", xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4), 
              pi2Div = 12, axes = FALSE, invertFlip = TRUE, nCores = 2)
riemannMask(annotNorth = TRUE)
par(op)

## ----figure_13, fig.width=7, fig.height=3.7, fig.align='center', results='hide', fig.show='hold', cache=TRUE, fig.cap='Riemann sphere plot of the function $f(z)=\\frac{(z^{2}+\\frac{1}{\\sqrt{2}}+\\frac{\\mathrm{i}}{\\sqrt{2}})\\cdot(z+\\frac{1}{2}+\\frac{\\mathrm{i}}{2})}{z-1}$. Annotated are the zeroes $z_1$, $z_2$, $z_3$, and the poles $z_4$, $z_5$.'----
op <- par(mfrow = c(1, 2), mar = c(0.1, 0.1, 1.4, 0.1))
# Define function
FUNstring <- "(z^2 + 1/sqrt(2) * (1 + 1i)) * (z + 1/2*(1 + 1i)) / (z - 1)"
# Southern hemisphere
phasePortrait(FUNstring, xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), 
              pi2Div = 12, axes = FALSE, nCores = 2)
riemannMask()
title("Southern Hemisphere", line = 0)
# - annotate zeroes and poles
text(c(cos(5/8*pi), cos(13/8*pi), cos(5/4*pi)/sqrt(2), 1),
     c(sin(5/8*pi), sin(13/8*pi), sin(5/4*pi)/sqrt(2), 0), 
     c(expression(z[1]), expression(z[2]), expression(z[3]), expression(z[4])), 
     pos = c(1, 2, 4, 2), offset = 1, col = "white")
# Northern hemisphere
phasePortrait(FUNstring, xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), 
              pi2Div = 12, axes = FALSE, invertFlip = TRUE, nCores = 2)
riemannMask()
title("Northern Hemisphere", line = 0)
# - annotate zeroes and poles
text(c(cos(5/8*pi), cos(13/8*pi), cos(5/4*pi)*sqrt(2), 1, 0),
     c(sin(5/8*pi), sin(13/8*pi), sin(5/4*pi)*sqrt(2), 0, 0), 
     c(expression(z[1]), expression(z[2]), expression(z[3]), 
       expression(z[4]), expression(z[5])), 
     pos = c(1, 4, 3, 4, 4), offset = 1, 
     col = c("white", "white", "black", "white", "white"))
par(op)

## ---- eval=FALSE--------------------------------------------------------------
#  x11(width = 8, height = 2/3 * 8)      # Open graphics window on screen
#  op <- par(mar = c(0, 0, 0, 0))        # Do not leave plot margins
#  phasePortrait(mandelbrot, moreArgs = list(itDepth = 30),
#                ncores = 1,             # Increase or leave out for higher performance
#                xlim = c(-2, 1), ylim = c(-1, 1),
#                hsvNaN = c(0, 0, 0),    # black color for points outside the set
#                axes = FALSE,           # No coordinate axes
#                xaxs = "i", yaxs = "i") # No space between plot region and plot
#  par(op)                               # Set graphics parameters to original

## ---- eval=FALSE--------------------------------------------------------------
#  res  <- 600 # set resolution to 600 dpi
#  # open png graphics device with in DIN A4 format
#  # DIN A format has an edge length ratio of sqrt(2)
#  png("Mandelbrot Example.png",
#      width = 29.7, height = 29.7/sqrt(2), # DIN A4 landscape
#      units = "cm",
#      res = res)                   # resolution is required
#  op   <- par(mar = c(0, 0, 0, 0)) # set graphics parameters - no plot margins
#  xlim <- c(-1.254, -1.248)        # horizontal (real) plot limits
#  # the function below adjusts the imaginary plot limits to the
#  #   desired ratio (sqrt(2)) centered around the desired imaginary value
#  ylim <- ylimFromXlim(xlim, centerY = 0.02, x_to_y = sqrt(2))
#  phasePortrait(mandelbrot,
#                nCores = 1,             # Increase or leave out for higher performance
#                xlim = xlim, ylim = ylim,
#                hsvNaN = c(0, 0, 0),    # Black color for NaN results
#                xaxs = "i", yaxs = "i", # suppress R's default axis margins
#                axes = FALSE,           # do not plot axes
#                res = res)              # resolution is required
#  par(op)   # reset graphics parameters
#  dev.off() # close graphics device and complete the png file

## ---- eval=FALSE--------------------------------------------------------------
#  res <- 600
#  png("Julia Example 1.png", width = 29.7, height = 29.7/sqrt(2),
#      units = "cm", res = res)
#  op <- par(mar = c(0, 0, 0, 0))
#  xlim <- c(-1.8, 1.8)
#  ylim <- ylimFromXlim(xlim, centerY = 0, x_to_y = sqrt(2))
#  phasePortrait(juliaNormal,
#                # see documentation of juliaNormal about the arguments
#                #  c and R_esc
#                moreArgs = list(c = -0.09 - 0.649i, R_esc = 2),
#                nCores = 1, # Increase or leave out for higher performance
#                xlim = xlim, ylim = ylim,
#                hsvNaN = c(0, 0, 0),
#                xaxs = "i", yaxs = "i",
#                axes = FALSE,
#                res = res)
#  par(op)
#  dev.off()

## ---- eval=FALSE--------------------------------------------------------------
#  res <- 600
#  png("Julia Example 2.png", width = 29.7, height = 29.7/sqrt(2),
#      units = "cm", res = res)
#  op <- par(mar = c(0, 0, 0, 0))
#  xlim <- c(-0.32, 0.02)
#  ylim <- ylimFromXlim(xlim, center = -0.78, x_to_y = sqrt(2))
#  phasePortrait(juliaNormal,
#                # see documentation of juliaNormal about the arguments
#                #  c and R_esc
#                moreArgs = list(c = -0.119 - 0.882i, R_esc = 2),
#                nCores = 1, # Increase or leave out for higher performance
#                xlim = xlim, ylim = ylim,
#                hsvNaN = c(0, 0, 0),
#                xaxs = "i", yaxs = "i",
#                axes = FALSE,
#                res = res)
#  par(op)
#  dev.off()

## ---- eval=FALSE--------------------------------------------------------------
#  # Map the complex plane on itself, show all bwType options
#  
#  x11(width = 8, height = 8)
#  op <- par(mfrow = c(2, 2), mar = c(4.1, 4.1, 1.1, 1.1))
#  for(bwType in c("ma", "a", "m")) {
#    phasePortraitBw("z", xlim = c(-2, 2), ylim = c(-2, 2),
#                    bwType = bwType,
#                    xlab = "real", ylab = "imaginary",
#                    nCores = 2) # Increase or leave out for higher performance
#  }
#  # Add normal phase portrait for comparison
#  phasePortrait("z", xlim = c(-2, 2), ylim = c(-2, 2),
#                xlab = "real", ylab = "imaginary",
#                pi2Div = 18,         # Use same angular division as default
#                                     # in phasePortraitBw
#                nCores = 2)     # Increase or leave out for higher performance
#  par(op)
#  

## ----eval=FALSE---------------------------------------------------------------
#  # A rational function, show all bwType options
#  
#  x11(width = 8, height = 8)
#  funString <- "(z + 1.4i - 1.4)^2/(z^3 + 2)"
#  op <- par(mfrow = c(2, 2), mar = c(4.1, 4.1, 1.1, 1.1))
#  for(bwType in c("ma", "a", "m")) {
#    phasePortraitBw(funString, xlim = c(-2, 2), ylim = c(-2, 2),
#                    bwType = bwType,
#                    xlab = "real", ylab = "imaginary",
#                    nCores = 2) # Increase or leave out for higher performance
#  }
#  # Add normal phase portrait for comparison
#  phasePortrait(funString, xlim = c(-2, 2), ylim = c(-2, 2),
#                xlab = "real", ylab = "imaginary",
#                pi2Div = 18,         # Use same angular division as default
#                                     # in phasePortraitBw
#                nCores = 2)     # Increase or leave out for higher performance
#  par(op)
#  

## ---- eval = FALSE------------------------------------------------------------
#  # Set the plot margins at all four sides to 1/5 inch with mai,
#  # set the background color to black with bg, and the default foreground
#  # color with fg (e.g. for axes and boxes around plots, or the color of
#  # the circle outline from the function riemannMask).
#  # We catch the previous parameter values in a variable, I called
#  # "op" ("old parameters")
#  op <- par(mai = c(1/5, 1/5, 1/5, 1/5), bg = "black", fg = "white")
#  
#  # Make any phase portraits and/or other graphics of your interest
#  # ...
#  
#  # Set the graphical parameters back to the values previously stored in op
#  par(op)

## ---- eval = FALSE------------------------------------------------------------
#  phasePortrait("tan(z^3 + 1/2 - 2i)/(1 - 1i - z)",
#                xlim = c(-6, 6), ylim = c(-3, 3),
#                axes = FALSE,
#                nCores = 2) # Increase or leave out for higher performance

## ---- eval=FALSE--------------------------------------------------------------
#  phasePortrait("tan(z^3 + 1/2 - 2i)/(1 - 1i - z)",
#                xlim = c(-6, 6), ylim = c(-3, 3),
#                axes = FALSE,
#                nCores = 2) # Increase or leave out for higher performance
#  box()

## ---- eval=FALSE--------------------------------------------------------------
#  # set background and foreground colors
#  op <- par(bg = "black", fg = "white")
#  # Setting the parameter fg has an effect on the box, the axes, and the axes'
#  # ticks, but not on the axis annotations and axis labels.
#  # Also the color of the title (main) is not affected.
#  # The colors of these elements have to be set manually and separately. While we
#  # could simply set them to "white", we set them, more flexibly, to the
#  # current foreground color (par("fg")).
#  phasePortrait("tan(z^3 + 1/2 - 2i)/(2 - 1i - z)",
#                xlim = c(-6, 6), ylim = c(-3, 3),            col.axis = par("fg"),
#                xlab = "real", ylab = "imaginary",           col.lab  = par("fg"),
#                main = "All annotation in foreground color", col.main = par("fg"),
#                # Adjust text size
#                cex.axis = 0.9, cex.lab = 0.9,
#                nCores = 2) # Increase or leave out for higher performance
#  par(op)

## ---- eval=FALSE--------------------------------------------------------------
#  # Open graphics device with 16/9 aspect ratio and 7 inch width
#  x11(width = 7, height = 9/16 * 7)
#  op <- par(mar = c(0, 0, 0, 0)) # Set plot margins to zero
#  xlim <- c(-3, 3)
#  # Calculate ylim with desired center fitting the desired aspect ratio
#  ylim <- ylimFromXlim(xlim, centerY = 0, x_to_y = 16/9)
#  phasePortrait(jacobiTheta, moreArgs = list(tau = 1i/5 + 1/5), pType = "p",
#                xlim = xlim, ylim = ylim,
#                xaxs = "i", yaxs = "i",
#                axes = FALSE,
#                nCores = 2) # Increase or leave out for higher performance
#  par(op)

## ---- eval=FALSE--------------------------------------------------------------
#  # Open graphics device with 16/9 aspect ratio and a width of 7 inches
#  x11(width = 7, height = 9/16 * 7)
#  # Set plot margins to zero, outer margins to 1/7 inch,
#  #   and background color to black
#  outerMar <- 1/7 # outer margin width in inches
#  op <- par(mar = c(0, 0, 0, 0), omi = rep(outerMar, 4), bg = "black")
#  xlim <- c(-1.5, 0.5)
#  # Calculate ylim with desired center fitting the desired aspect ratio;
#  #   however, the omi settings slightly change the required
#  #   ratio of xlim and ylim
#  ratio <- (7 - 2*outerMar) / (7 * 9/16 - 2*outerMar)
#  ylim  <- ylimFromXlim(xlim, centerY = 0, x_to_y = ratio)
#  phasePortrait("sin(jacobiTheta(z, tau))/z", moreArgs = list(tau = 1i/5 + 1/5),
#                pType = "p",
#                xlim = xlim, ylim = ylim,
#                xaxs = "i", yaxs = "i",
#                axes = FALSE,
#                nCores = 1) # Increase or leave out for higher performance
#  par(op)

## ---- eval = FALSE------------------------------------------------------------
#  # Note that 'FUN =' is not required if the argument to FUN is handed to
#  # phasePortrait in the first position
#  phasePortrait(FUN = "1/(1 - z^2)", xlim = c(-5, 5), ylim = c(-5, 5), nCores = 2)
#  phasePortrait("sin((z - 2)/(z + 2))", xlim = c(-5, 5), ylim = c(-5, 5), nCores = 2)
#  phasePortrait("tan(z)", xlim = c(-5, 5), ylim = c(-5, 5), nCores = 2)

## ---- eval = FALSE------------------------------------------------------------
#  phasePortrait("-1 * sum(z^c(-k:k))", moreArgs = list(k = 11),
#                xlim = c(-2, 2), ylim = c(-1.5, 1.5),
#                pType = "p",
#                nCores = 2) # Increase or leave out for higher performance

## ---- eval = FALSE------------------------------------------------------------
#  funString <- "vapply(z, FUN = function(z) {
#                    n <- 9
#                    k <- z^(c(1:n))
#                    rslt <- sum(sin(k))
#                    return(rslt)
#                  },
#                  FUN.VALUE = complex(1))"
#  phasePortrait(funString, xlim = c(-2, 2), ylim = c(-2, 2),
#                nCores = 2) # Increase or leave out for higher performance

## ---- eval = FALSE------------------------------------------------------------
#  funString <- "vapply(z, FUN = function(z, fct) {
#                    n <- 9
#                    k <- z^(fct * c(1:n))
#                    rslt <- sum(sin(k))
#                    return(rslt)
#                  },
#                  fct = -1,
#                  FUN.VALUE = complex(1))"
#  phasePortrait(funString, xlim = c(-2, 2), ylim = c(-2, 2),
#                nCores = 2) # Increase or leave out for higher performance

## ---- eval = FALSE------------------------------------------------------------
#  # Define function
#  tryThisOne <- function(z, fct, n) {
#    k <- z^(fct * c(1:n))
#    rslt <- prod(cos(k))
#    return(rslt)
#  }
#  
#  # Call function by its name only, provide additional arguments via "moreArgs"
#  phasePortrait("tryThisOne", moreArgs = list(fct = 1, n = 5),
#    xlim = c(-2.5, 2.5), ylim = c(-2, 2),
#    nCores = 2) # Increase or leave out for higher performance

## ---- eval = FALSE------------------------------------------------------------
#  # Use argument "hsvNaN = c(0, 0, 0)" if you want the grey area black
#  phasePortrait(function(z) {
#      for(j in 1:20) {
#        z <- z * sin(z) - 1 + 1/2i
#      }
#      return(z)
#    },
#    xlim = c(-3, 3), ylim = c(-2, 2),
#    nCores = 2) # Increase or leave out for higher performance

## ---- eval = FALSE------------------------------------------------------------
#  # Use argument "hsvNaN = c(0, 0, 0)" if you want the grey area black
#  phasePortrait(function(z, n) {
#      for(j in 1:n) {
#        z <- z * cos(z)
#      }
#      return(z)
#    },
#    moreArgs = list(n = 27),
#    xlim = c(-3, 3), ylim = c(-2, 2),
#    nCores = 2) # Increase or leave out for higher performance

## ---- eval = FALSE------------------------------------------------------------
#  # atan from package base
#  phasePortrait(atan, xlim = c(-pi, pi), ylim = c(-pi, pi),
#                nCores = 2)
#  
#  # gammaz from package pracma (the package must be installed on your machine
#  # if you want this example to be working)
#  phasePortrait(pracma::gammaz, xlim = c(-9, 9), ylim = c(-5, 5),
#                nCores = 2)
#  
#  # blaschkeProd from this package (moreArgs example)
#  #  make random vector of zeroes
#  n <- 12
#  a <- complex(modulus = runif(n), argument = 2 * pi * runif(n))
#  #  plot the actual phase portrait
#  phasePortrait(blaschkeProd, moreArgs = list(a = a),
#                xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3),
#                nCores = 2)
#  
#  # User function example
#  tryThisOneToo <- function(z, n, r) {
#    for(j in 1:n) {
#      z <- r * (z + z^2)
#    }
#    return(z)
#  }
#  # Use argument "hsvNaN = c(0, 0, 0)" if you want the gray areas black
#  phasePortrait(tryThisOneToo, moreArgs = list(n = 50, r = 1/2 - 1/2i),
#                xlim = c(-3, 2), ylim = c(-2.5, 2.5),
#                nCores = 2)
#  

## ---- eval = FALSE------------------------------------------------------------
#  res <- 300 # Define desired resolution in dpi
#  png("Logistic_Function.png", width = 40, height = 40 * 3/4,
#      units = "cm", res = res)
#  phasePortrait("1/(1+exp(-z))", xlim = c(-25, 25), ylim = c(-15, 15), res = res,
#                xlab = "real", ylab = "imaginary",
#                nCores = 2) # Increase or leave out for higher performance
#  dev.off()

## ---- eval=FALSE--------------------------------------------------------------
#  switch(1 + trunc(runif(1, 0, 6)),
#         "... at all?",
#         "... in a quick-and-dirty way?",
#         "... in Hadley-Wickham-style?",
#         "... without a loop?",
#         "... without nested loops?",
#         "... in a way somebody can understand?")

## ---- include = FALSE---------------------------------------------------------
foreach::registerDoSEQ()

