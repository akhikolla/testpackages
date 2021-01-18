// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif





// recode PLINK as ASCII with no spaces
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool  CreateASCIInospace_PLINK(std::string fname, std::string asciifname, std::vector<long> dims,
                         bool quiet, Rcpp::Function message)
{

int
  n_of_cols_in_geno = (dims[1] -6)/2.0;

std::vector<char> alleles0( n_of_cols_in_geno );
std::vector<char> alleles1( n_of_cols_in_geno );


std::vector<char>
     rowvec( dims[1] - 6 );  // holds allelic information from PLINK file


std::string
   tmp,
   token,
   line;

 std::ostringstream
      os;

// open PLINK ped  file
std::ifstream fileIN(fname.c_str());
if(!fileIN.good()) {
  message("ERROR: PLINK ped file could not be opened with filename  ",   fname );
  message("ERROR: ReadMarkerData has terminated with errors.  ");
  return false;
}

// open ascii file that is to hold no-spaces genotype data
std::ofstream fileOUT(asciifname.c_str(), std::ios::out );
long  counter = 0;

// initializing input line 
std::string rowinfile(n_of_cols_in_geno, '0'); // s == "000000"

int printOnlyOnce = 0;  // flag for printing warning message about missing data

while(getline(fileIN, line ))
{
  std::istringstream streamLine(line);

 // check number of columns for each line
 std::string rowinfile(n_of_cols_in_geno, '0'); // s == "000000"


   std::istringstream check_number_row_elements(line);

   long numcols = std::distance(std::istream_iterator<std::string>(check_number_row_elements),
              std::istream_iterator<std::string>()) ;

   if (numcols  != dims[1] ){
         message("\n");
         message( "Error:  PLINK file contains an unequal number of columns per row.  " );
         message( "        The error has occurred at row " ,  counter+1 , " which contains " ,  numcols  ,  " but " );
         message( "        it should contain " , dims[1] , " columns of data. " );
         message("\n");
         message(" ReadMarkerData has terminated with errors");
         return false;
   }  // end  if (number_of_columns != dims[1] )





   std::istringstream streamA(line);



   // tokenized row and placed it in std::vector rowvec
   for(int i=0; i <= 5; i++){
     streamA >> tmp;
   }
   for(long i=6; i < dims[1] ; i++){
            streamA >> rowvec[i-6];
   }  // end  for(long i=0; i < dims[1] ; i++)



   // initialize alleles structure to first row of PLINK info
   if (counter == 0) {
         for(long i=0; i < n_of_cols_in_geno ; i++){
            if( rowvec[ (2*i ) ] == '0' ||  rowvec[ (2*i + 1) ] == '0' || rowvec[ (2*i ) ] == '-' ||  rowvec[ (2*i + 1) ] == '-'){
               // missing allele
               alleles0[ i ] = 'I';
               alleles1[ i ] = 'I';
            } else {
               alleles0[ i ] =  rowvec[ (2*i ) ];
               alleles1[ i ] =  rowvec[ (2*i + 1) ];
            } //end  if (rowvec 
         }
   }

   // turn allelic info from PLINK into genotype 0,1,2 data
     // also do some checks for more than 2 alleles, and 0 and - for missing data
     for(long i=0; i <  n_of_cols_in_geno; i++){
        // Checking for missing allelic information in PLINK file

        if( rowvec[ (2*i ) ] == '0' ||  rowvec[ (2*i + 1) ] == '0' || rowvec[ (2*i ) ] == '-' ||  rowvec[ (2*i + 1) ] == '-'){
           if (printOnlyOnce == 0){
                 message("\n");
                 message(" Warning:  PLINK file contains missing alleles (i.e. 0 or - ) " );
                 message("           These missing genotypes should be imputed before running Eagle." );
                 message("           As an approximation, AMpus has set these missing genotypes to heterozygotes. " );
                 message("           Since Eagle assumes an additive model, heterozygote genotypes do not contribute to the estimation of " );
                 message("           the additive effects.  " );
                 message("\n");
                 printOnlyOnce = 1;
            } // if printOnlyOnce`
            rowvec[ (2*i) ] = 'I';      // impute
            rowvec[ (2*i + 1) ] = 'I';  // impute
        }

        // Check if allele has been seen before in allele file. 
        // If so, make sure alleles doesn't already  contain two alleles - otherwise generate error message
        for(int j = 1; j >= 0; --j){ // looping over the two alleles with indexes 0 and 1


           if (rowvec[ (2*i + j) ] != alleles0[ i ] && rowvec[ (2*i + j) ] != alleles1[ i ]){
              // situation 1: rowvec contains missing values ie 'I' then do nothing



              if (rowvec[ (2*i + j) ] == 'I' ){
                 // do nothing here

              } else {
                // situation 2: alleles contain missing values I
                if (alleles0[i] == 'I'){
                     alleles0[i] = rowvec[ (2*i + j) ];
                } else {
                      if (alleles1[i] == 'I'){
                           alleles1[i] = rowvec[ (2*i + j) ];
                       } else {
                         if (alleles0[ i ] == alleles1[ i ] ){
                             // this is okay. alleles only contains a single allele at the moment. Re-initialise alleles
                             alleles1[ i ] = rowvec[ (2*i + j) ];
                           } else {
                              // Error - we have more than two alleles segregating at a locus
                            message("\n");
                            message("Error:  PLINK file cannot contain more than two alleles at a locus.");
                            message("        The error has occurred at snp locus " , i  + 1 , " for individual " , counter+1 );
                            message("\n");
                            message(" ReadMarkerData has terminated with errors");
                            return false;
                          } // end inner if else

                       } // end if (alleles1[i] == 'I')
                } // end  if (alleles0[i] == 'I')

              } // end if (rowvec[ (2*i + j) ] == 'I' )


           }  // end if (rowvec[ (2*i + j) ] != alleles0[ i ] && rowvec[ (2*i + j) ] != alleles1[ i ])




    // set rowinfile
    if (rowvec[ (2*i) ] == 'I' || rowvec[ (2*i + 1) ] == 'I' ){
        rowinfile[i] = '1' ; // AB geno no additive effect
    } else {
        if (rowvec[ (2*i + 1) ] !=   rowvec[ (2*i) ] ){
          rowinfile[i] = '1' ;  // AB
        } else {
          if (rowvec[ (2*i ) ] == alleles0[ i ] ){  // matches first allele
               rowinfile[i] = '0';  // AA
          }  else {
               rowinfile[i] = '2';  // BB
          }
        }  // end outer if else rowvec
    } // end  if (rowvec[ (2*i) ] == 'I' || rowvec[ (2*i + 1) ] == 'I' )
 } // end  for(long i=0; i < n_of_cols_in_geno; i++)




  }  // end for(long i=0; i< n_of_cols_in_geno ; i++)

  fileOUT << rowinfile;
  fileOUT << "\n";
  counter++;



  }  // end while(getline(fileIN, line ))



// write out a few lines of the file if quiet
// open PLINK ped  file
std::ifstream fileIN_backtobeginning(fname.c_str());
counter = 0;


int nrowsp =  5;
int ncolsp = 24;
if(dims[0] < 5)
nrowsp = dims[0];
if(dims[1] < 25)
ncolsp = dims[1];

message(" First ", nrowsp, " lines and ", ncolsp, " columns of the PLINK ped file. ");


std::string rowline;

while(getline(fileIN_backtobeginning, line ) && counter < nrowsp)
{
       std::ostringstream oss;
       std::istringstream streamB(line);
       for(int i=0; i < ncolsp ; i++){
           streamB >> tmp;
           oss << tmp << " " ;
        }
        std::string rowline = oss.str();
        message(rowline);
        counter++;
}  // end  while(getline(fileIN, line ))




// close files
fileIN.close();
fileOUT.close();

return true;

}


