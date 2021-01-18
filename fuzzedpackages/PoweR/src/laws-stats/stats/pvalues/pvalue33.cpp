tmp = Rf_pnorm5(statistic[0], 1.0, sqrt((M_PI / 2.0 - 1.5) / ((double)n)), 1, 0);
if (tmp > 0.5) tmp = 1.0 - tmp;
if (alter[0] == 0) pvalue[0] = 2 * tmp; else pvalue[0] = tmp; 
