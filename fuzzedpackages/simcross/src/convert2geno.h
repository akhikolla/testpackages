
IntegerMatrix convert2geno(const List xodat, const NumericVector map);

IntegerVector convertchr2geno(const List chr, const NumericVector map);

IntegerMatrix combine_mat_and_pat_geno(const IntegerMatrix matmatrix, const IntegerMatrix patmatrix,
                                       const int max_geno);

IntegerMatrix combine_mat_and_pat_geno_wfounders(const IntegerMatrix matmatrix, const IntegerMatrix patmatrix,
                                                 const IntegerMatrix founder_geno);

IntegerVector convert2genoarray(const List xodat, const NumericVector map);

CharacterMatrix convert2geno_char(const List xodat, const NumericVector map, const CharacterMatrix founder_geno);
CharacterMatrix combine_mat_and_pat_geno_wfounders_char(const IntegerMatrix matmatrix, const IntegerMatrix patmatrix, 
							const CharacterMatrix founder_geno);
CharacterMatrix convert2geno_char_paste(const List xodat, const NumericVector map, const CharacterMatrix founder_geno);
CharacterMatrix combine_mat_and_pat_geno_wfounders_char_paste(const IntegerMatrix matmatrix, const IntegerMatrix patmatrix, 
							      const CharacterMatrix founder_geno);
