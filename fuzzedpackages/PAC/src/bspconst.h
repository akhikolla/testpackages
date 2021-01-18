#ifndef BSP_CONST_H
#define BSP_CONST_H

static const uint MAX_LEVEL = 10;
static const int NUM_CUT = 4;
static const uint MIN_PTS = 50; 	/* FOR SLT */
static const uint NUM_SPLITS = 10;  /* for splitting history */
static const double LAPLACE = 10.0;

static const double ALPHA = 2;
static const double SIGMA = 1;
static const double GAMMA = 1;
static const double THRESHOLD = 4;

static const double THETA = 1;

static const bool RESAMPLE = true;
static const uint RESAMPLE_STEPS = 4;
static const uint NUM_SIS_SAMPLES = 100;
static const double MIThreshold = 0.01;
static const uint MI_NUM_SAMPLE = 5E4;
static const uint MAX_NPTS = 100;


#endif /* BSP_CONST_H */
