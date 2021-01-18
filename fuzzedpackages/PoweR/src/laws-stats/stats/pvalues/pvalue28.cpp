double q2 = 4913.0 / 1458.0; // Median
if (alter[0] == 0) {
  if (statistic[0] >= q2) {
    pvalue[0] = 2.0 * pchisq(statistic[0], 4.0, 0, 0);
  } else {
    pvalue[0] = 2.0 * pchisq(statistic[0], 4.0, 1, 0);
  }
 }
if (alter[0] == 1) {
  pvalue[0] =  pchisq(statistic[0], 4.0, 1, 0); 
 }
if (alter[0] == 2) {
  pvalue[0] =  pchisq(statistic[0], 4.0, 0, 0); 
 }
