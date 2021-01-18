// on doit pouvoir utiliser Cornish-Fisher pour calculer la p-value tmp (voir l'article de D'Agostino)
// Il faut aussi calculer la mediane q2
/*
if (alter[0] == 0) {
  if (statistic[0] >= q2) {
    pvalue[0] = 2 * tmp1; // upper tail
  } else {
    pvalue[0] = 2 * tmp2; // lower tail
  }
 }
if (alter[0] == 1) {
  pvalue[0] = tmp; // lower tail
 }
if (alter[0] == 2) {
  pvalue[0] = tmp; // upper tail
 }
*/
pvalcomp[0] = 0; // a virer quand je saurais comment calculer la p-valeur
