## RHTDP [![Build Status](https://travis-ci.org/jbuonagurio/RHTDP.svg)](https://travis-ci.org/jbuonagurio/RHTDP)

RHTDP provides R bindings to the NGS [Horizontal Time-Dependent Positioning](https://www.ngs.noaa.gov/TOOLS/Htdp/Htdp.shtml) utility, v3.2.5, written by Richard Snay, Chris Pearson, and Jarir Saleh of NGS. HTDP allows users to transform positional coordinates across time and between spatial reference frames.

RHTDP is a work-in-progress; currently only crustal displacement and velocity calculations are implemented. Transformations between reference frames are planned, as well as integration with the [sp](https://cran.r-project.org/package=sp) package.

### Installation

```R
# Development version from GitHub:
# install.packages("devtools")
devtools::install_github("jbuonagurio/RHTDP")
```

### Usage

Following exercises 1 and 2 in the HTDP User's Guide, the following example estimates the displacement that occurred at point 'beta' during the [Loma Prieta earthquake](https://en.wikipedia.org/wiki/1989_Loma_Prieta_earthquake) (M=7.1) of October 18, 1989. The NAD_83(2011/CORS96/2007) reference frame is used (iopt=1).

```R
beta <- matrix(data=c(-121.7722, 36.6698), nrow=1)
t0 <- as.Date("1989-10-16")
t1 <- as.Date("1989-10-18")
displace(beta, t0, t1, 1)
```

R output:
```
            de         dn           du        ve       vn        vu
1 -0.001449719 0.07401263 -0.004262457 -25.82612 37.15177 -1.329391
```

HTDP displacement:
```
 HTDP (VERSION v3.2.5    ) OUTPUT

 DISPLACEMENTS IN METERS RELATIVE TO NAD_83(2011/CORS96/2007)
 FROM 10-16-1989 TO 10-18-1989 (month-day-year)
 FROM 1989.789 TO 1989.795 (decimal years)

NAME OF SITE             LATITUDE          LONGITUDE            NORTH    EAST    UP 
beta                     36 40 11.28000 N  121 46 19.92000 W    0.074  -0.001  -0.004
```

HTDP velocity:
```
 HTDP (VERSION v3.2.5    ) OUTPUT

VELOCITIES IN MM/YR RELATIVE TO NAD_83(2011/CORS96/2007)


          beta                    
LATITUDE   =  36 40 11.28000 N  NORTH VELOCITY =  37.15 mm/yr
LONGITUDE  = 121 46 19.92000 W  EAST VELOCITY  = -25.83 mm/yr
ELLIPS. HT. =      0.000 m      UP VELOCITY    =  -1.33 mm/yr
X = -2696934.816 m              X VELOCITY     =  -9.71 mm/yr
Y = -4354426.684 m              Y VELOCITY     =  33.37 mm/yr
Z =  3788064.740 m              Z VELOCITY     =  29.01 mm/yr
```

Supported values for iopt are below; see the User's Guide for more information. Alternatively, use the `iopt()` function to list these within R. Note that option 4 (WGS_72) has been removed.

| Key | Reference Frame |
| ---: | :--- |
| 1 | NAD_83(2011/CORS96/2007) (North American plate fixed) |
| 2 | NAD_83(PA11/PACP00) (Pacific plate fixed) |
| 3 | NAD_83(MA11/MARP00) (Mariana plate fixed) |
| 5 | WGS_84(transit) (NAD_83(2011) used) |
| 6 | WGS_84(G730) (ITRF91 used) |
| 7 | WGS_84(G873) (ITRF94 used) |
| 8 | WGS_84(G1150) (ITRF2000 used) |
| 9 | WGS_84(G1674) (ITRF2008 used) |
| 10 | WGS_84(G1762) (IGb08 used) |
| 11 | SIO/MIT_92 (ITRF91 used) |
| 12 | ITRF88 |
| 13 | ITRF89 |
| 14 | ITRF90 or (PNEOS90/NEOS90) |
| 15 | ITRF91 |
| 16 | ITRF92 |
| 17 | ITRF93 |
| 18 | ITRF94 |
| 19 | ITRF96 |
| 20 | ITRF97 or IGS97 |
| 21 | ITRF2000 or IGS00/IGb00 |
| 22 | ITRF2005 or IGS05 |
| 23 | ITRF2008 or IGS08/IGb08 |

### Warning

HTDP exclusively uses 14-parameter transformations between reference frames. If you require rigorous transformations between various realizations of NAD83, use NADCON/GEOCON/GEOCON11.

Chuck Ghilani of Penn State has an [excellent presentation](http://www.personal.psu.edu/cdg3/PSLS.pdf) covering the recommended transformation process between different epochs and datums.
