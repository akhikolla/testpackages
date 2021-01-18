if (alter[0] == 0) pvalue[0] = 2.0 * Rf_pnorm5(fabs(statistic[0]), 0.0, sqrt(165.0 / 149.0), 0, 0);
if (alter[0] == 1) pvalue[0] = Rf_pnorm5(statistic[0], 0.0, sqrt(165.0 / 149.0), 1, 0);
if (alter[0] == 2) pvalue[0] = Rf_pnorm5(statistic[0], 0.0, sqrt(165.0 / 149.0), 0, 0);
