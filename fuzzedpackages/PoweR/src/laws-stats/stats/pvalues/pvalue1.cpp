pval = exp(-7.01256 * Kd * Kd * ((double)nd + 2.78019) + 2.99587 * Kd * sqrt((double)nd + 2.78019) -
	   0.122119 + 0.974598 / sqrt((double)nd) + 1.67997 / (double)nd);

if (pval > 0.1) {
  KK = (sqrt((double)n) - 0.01 + 0.85 / sqrt((double)n)) * statKS;
  if (KK <= 0.302) {
    pval = 1.0;
  } else if (KK <= 0.5) {
    KK2 = KK * KK;
    KK3 = KK2 * KK;
    KK4 = KK3 * KK;
    pval = 2.76773 - 19.828315 * KK + 80.709644 * KK2 - 138.55152 * KK3 + 81.218052 * KK4;
  } else if (KK <= 0.9) {
    KK2 = KK * KK;
    KK3 = KK2 * KK;
    KK4 = KK3 * KK;
    pval = -4.901232 + 40.662806 * KK - 97.490286 * KK2 + 94.029866 * KK3 - 32.355711 * KK4;
  } else if (KK <= 1.31) {
    KK2 = KK * KK;
    KK3 = KK2 * KK;
    KK4 = KK3 * KK;
    pval = 6.198765 - 19.558097 * KK + 23.186922 * KK2 - 12.234627 * KK3 + 2.423045 * KK4;
  } else {
    pval = 0.0;
  }
}

pvalue[0] = pval;
