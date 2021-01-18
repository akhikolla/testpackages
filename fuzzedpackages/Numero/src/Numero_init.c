#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP nro_circus_paint(SEXP offsets_R, SEXP topo_R, SEXP ccodes_R, SEXP key_R, SEXP title_R);
SEXP nro_circus_show(SEXP offsets_R, SEXP topo_R, SEXP ccodes_R, SEXP labels_R, SEXP key_R);
SEXP nro_circus_write(SEXP offsets_R, SEXP topo_R, SEXP labels_R, SEXP visible_R, SEXP contrast_R, SEXP key_R, SEXP font_R);
SEXP nro_colorize(SEXP zvals_R, SEXP name_R);
SEXP nro_destratify(SEXP data_R, SEXP strata_R);
SEXP nro_diffuse(SEXP topo_R, SEXP sigma_R, SEXP bmus_R, SEXP data_R);
SEXP nro_figure(SEXP fname_R, SEXP data_R, SEXP bbox_R, SEXP script_R);
SEXP nro_impute(SEXP xdata_R, SEXP nsub_R, SEXP lag_R);
SEXP nro_kohonen(SEXP seeds_R, SEXP rho_R, SEXP sigma_R);
SEXP nro_label(SEXP topo_R, SEXP data_R, SEXP binflags_R, SEXP sigma_R);
SEXP nro_match(SEXP codebook_R, SEXP data_R);
SEXP nro_permute(SEXP topo_R, SEXP sigma_R, SEXP bmus_R, SEXP data_R, SEXP numcycl_R, SEXP lag_R);
SEXP nro_train(SEXP topo_R, SEXP sigma_R, SEXP codebook_R, SEXP data_R, SEXP nsub_R, SEXP eq_R, SEXP lag_R);
SEXP nro_webpage(SEXP fname_R, SEXP bytes_R);

R_CallMethodDef callMethods[]  = {
  {"nro_circus_paint", (DL_FUNC) &nro_circus_paint, 5},
  {"nro_circus_show", (DL_FUNC) &nro_circus_show, 5},
  {"nro_circus_write", (DL_FUNC) &nro_circus_write, 7},
  {"nro_colorize", (DL_FUNC) &nro_colorize, 2},
  {"nro_destratify", (DL_FUNC) &nro_destratify, 2},
  {"nro_diffuse", (DL_FUNC) &nro_diffuse, 4},
  {"nro_figure", (DL_FUNC) &nro_figure, 4},
  {"nro_impute", (DL_FUNC) &nro_impute, 3},
  {"nro_kohonen", (DL_FUNC) &nro_kohonen, 3},
  {"nro_label", (DL_FUNC) &nro_label, 4},
  {"nro_match", (DL_FUNC) &nro_match, 2},
  {"nro_permute", (DL_FUNC) &nro_permute, 6},
  {"nro_train", (DL_FUNC) &nro_train, 7},
  {"nro_webpage", (DL_FUNC) &nro_webpage, 2},
  {NULL, NULL, 0}
};

void R_init_Numero(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
