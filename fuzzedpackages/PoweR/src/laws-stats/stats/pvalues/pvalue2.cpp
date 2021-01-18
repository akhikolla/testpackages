statAD = statAD * (1.0 + 0.75/(double)n + 2.25 / (double)(n * n));
if (statAD < 0.2) {
  pval = 1 - exp(-13.436 + 101.14 * statAD - 223.73 * R_pow(statAD, 2.0));
 }
 else if (statAD < 0.34) {
   pval = 1 - exp(-8.318 + 42.796 * statAD - 59.938 * R_pow(statAD, 2.0));
 }
 else if (statAD < 0.6) {
   pval = exp(0.9177 - 4.279 * statAD - 1.38 * R_pow(statAD, 2.0));
 }
 else {
   pval = exp(1.2937 - 5.709 * statAD + 0.0186 * R_pow(statAD, 2.0));
 }

pvalue[0] = pval;
