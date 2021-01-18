# image.CornerDetectionF9 - Detect Corners in images 

The  **image.CornerDetectionF9** R package detects corners in images, namely:

- An implementation of the "FAST-9" corner detection algorithm explained at <http://www.edwardrosten.com/work/fast.html>. The package allows to detect corners in digital images. 
- More details in the paper at <https://arxiv.org/pdf/0810.2434.pdf>

## Examples

Read in an image with values in the 0-255 range (pgm image: http://netpbm.sourceforge.net/doc/pgm.html)

```r
library(pixmap)
library(image.CornerDetectionF9)
imagelocation <- system.file("extdata", "chairs.pgm", package = "image.CornerDetectionF9")
image   <- read.pnm(file = imagelocation, cellres = 1)
x       <- image@grey * 255
corners <- image_detect_corners(x, 100)
plot(image)
points(corners$x, corners$y, col = "red", pch = 20, lwd = 0.5)
```

![](https://raw.githubusercontent.com/bnosac/image/master/image.CornerDetectionF9/inst/extdata/result-f9.png?raw=true)


## Support in image recognition

Need support in image recognition?
Contact BNOSAC: http://www.bnosac.be