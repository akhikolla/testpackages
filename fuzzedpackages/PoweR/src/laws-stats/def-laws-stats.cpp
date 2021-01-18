// Below, add an: #include "laws/lawj.cpp"

#include "laws/law1.cpp"
#include "laws/law2.cpp"
#include "laws/law3.cpp"
#include "laws/law4.cpp"
#include "laws/law5.cpp"
#include "laws/law6.cpp"
#include "laws/law7.cpp"
#include "laws/law8.cpp"
#include "laws/law9.cpp"
#include "laws/law10.cpp"
#include "laws/law11.cpp"
#include "laws/law12.cpp"
#include "laws/law13.cpp"
#include "laws/law14.cpp"
#include "laws/law15.cpp"
#include "laws/law16.cpp"
#include "laws/law17.cpp"
#include "laws/law18.cpp"
#include "laws/law19.cpp"
#include "laws/law20.cpp"
#include "laws/law21.cpp"
#include "laws/law22.cpp"
#include "laws/law23.cpp"
#include "laws/law24.cpp"
#include "laws/law25.cpp"
#include "laws/law26.cpp"
#include "laws/law27.cpp"
#include "laws/law28.cpp"
#include "laws/law29.cpp"
#include "laws/law30.cpp"
#include "laws/law31.cpp"
#include "laws/law32.cpp"
#include "laws/law33.cpp"
#include "laws/law34.cpp"
#include "laws/law35.cpp"
#include "laws/law36.cpp"
#include "laws/law37.cpp"
#include "laws/law38.cpp"
#include "laws/law39.cpp"
#include "laws/law40.cpp"

// Modify also below.
// Change:
//   40 to 41 
// and add:
//  ,law41

void (*lawfunc[40])(int *xlen, double *x, char **name, int *trtname, double *params, int *nbparams, int *setseed) = {
  law1,law2,law3,law4,law5,law6,law7,law8,law9,law10,
  law11,law12,law13,law14,law15,law16,law17,law18,law19,law20,
  law21,law22,law23,law24,law25,law26,law27,law28,law29,law30,
  law31,law32,law33,law34,law35,law36,law37,law38,law39,law40
};


// Below, add an: #include "stats/statj.cpp"

/*
Note: stat84 are available and can be used to incorporate new tests
 */
#include "stats/stat1.cpp"
#include "stats/stat2.cpp"
#include "stats/stat3.cpp"
#include "stats/stat4.cpp"
#include "stats/stat5.cpp"
#include "stats/stat6.cpp"
#include "stats/stat7.cpp"
#include "stats/stat8.cpp"
#include "stats/stat9.cpp"
#include "stats/stat10.cpp"
#include "stats/stat11.cpp"
#include "stats/stat12.cpp"
#include "stats/stat13.cpp"
#include "stats/stat14.cpp"
#include "stats/stat15.cpp"
#include "stats/stat16.cpp"
#include "stats/stat17.cpp"
#include "stats/stat18.cpp"
#include "stats/stat19.cpp"
#include "stats/stat20.cpp"
#include "stats/stat21.cpp"
#include "stats/stat22.cpp"
#include "stats/stat23.cpp"
#include "stats/stat24.cpp"
#include "stats/stat25.cpp"
#include "stats/stat26.cpp"
#include "stats/stat27.cpp"
#include "stats/stat28.cpp"
#include "stats/stat29.cpp"
#include "stats/stat30.cpp"
#include "stats/stat31.cpp"
#include "stats/stat32.cpp"
#include "stats/stat33.cpp"
#include "stats/stat34.cpp"
#include "stats/stat35.cpp"
#include "stats/stat36.cpp"
#include "stats/stat37.cpp"
#include "stats/stat38.cpp"
#include "stats/stat39.cpp"
#include "stats/stat40.cpp"
#include "stats/stat41.cpp"
#include "stats/stat42.cpp"
#include "stats/stat43.cpp"
#include "stats/stat44.cpp"
#include "stats/stat45.cpp"
#include "stats/stat46.cpp"
#include "stats/stat47.cpp"
#include "stats/stat48.cpp"
#include "stats/stat49.cpp"
#include "stats/stat50.cpp"
#include "stats/stat51.cpp"
#include "stats/stat52.cpp"
#include "stats/stat53.cpp"
#include "stats/stat54.cpp"
#include "stats/stat55.cpp"
#include "stats/stat56.cpp"
#include "stats/stat57.cpp"
#include "stats/stat58.cpp"
#include "stats/stat59.cpp"
#include "stats/stat60.cpp"
#include "stats/stat61.cpp"
#include "stats/stat62.cpp"
#include "stats/stat63.cpp"
#include "stats/stat64.cpp"
#include "stats/stat65.cpp"
#include "stats/stat66.cpp"
#include "stats/stat67.cpp"
#include "stats/stat68.cpp"
#include "stats/stat69.cpp"
#include "stats/stat70.cpp"
#include "stats/stat71.cpp"
#include "stats/stat72.cpp"
#include "stats/stat73.cpp"
#include "stats/stat74.cpp"
#include "stats/stat75.cpp"
#include "stats/stat76.cpp"
#include "stats/stat77.cpp"
#include "stats/stat78.cpp"
#include "stats/stat79.cpp"
#include "stats/stat80.cpp"
#include "stats/stat81.cpp"
#include "stats/stat82.cpp"
#include "stats/stat83.cpp"
// #include "stats/stat84.cpp"
#include "stats/stat85.cpp"
#include "stats/stat86.cpp"
#include "stats/stat87.cpp"
#include "stats/stat88.cpp"
#include "stats/stat89.cpp"
#include "stats/stat90.cpp"
#include "stats/stat91.cpp"
#include "stats/stat92.cpp"
#include "stats/stat93.cpp"
#include "stats/stat94.cpp"
#include "stats/stat95.cpp"
#include "stats/stat96.cpp"
#include "stats/stat97.cpp"


// I added this to correct a bug that occured when calling (*statfunc[statindex-1]) from calcpuiss.cpp. This bug arrived when I removed stat84, etc below
void nothing(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {
  return;
}

// Modify also below.
// Change:
//   97 to 98
// and add:
//  ,stat98

void (*statfunc[97])(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) = {
  stat1,stat2,stat3,stat4,stat5,stat6,stat7,stat8,stat9,stat10,
  stat11,stat12,stat13,stat14,stat15,stat16,stat17,stat18,stat19,stat20,
  stat21,stat22,stat23,stat24,stat25,stat26,stat27,stat28,stat29,stat30,
    stat31,stat32,stat33,stat34,stat35,stat36,stat37,stat38,stat39,stat40,
  stat41,stat42,stat43,stat44,stat45,stat46,stat47,stat48,stat49,stat50,
  stat51,stat52,stat53,stat54,stat55,stat56,stat57,stat58,stat59,stat60,
  stat61,stat62,stat63,stat64,stat65,stat66,stat67,stat68,stat69,stat70,
  stat71,stat72,stat73,stat74,stat75,stat76,stat77,stat78,stat79,stat80,
    stat81,stat82,stat83,nothing,//stat84,
    stat85,stat86,stat87,stat88,stat89,stat90,
    stat91,stat92,stat93,stat94,stat95,stat96,stat97
};
