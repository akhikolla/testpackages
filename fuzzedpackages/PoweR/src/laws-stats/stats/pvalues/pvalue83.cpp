if (alter[0] == 0) {
  pvaltmp = 2.0 * pt(fabs(statistic[0]), n - 1, 1, 0);
 }
 else if (alter[0] == 1) {
   pvaltmp = pt(statistic[0], n - 1, 1, 0);
 }
 else if (alter[0] == 2) {
   pvaltmp = pt(statistic[0], n - 1, 0, 0);
 }

pvalue[0] = pvaltmp;
