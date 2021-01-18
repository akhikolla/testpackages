// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include  "createM_BIN_rcpp.h"


#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif

const size_t bits_in_int = std::numeric_limits<int>::digits;


 Rcpp::DataFrame createMapDF(     std::string fname,
                                  Rcpp::IntegerVector  shouldsnpbremoved,
                                  Rcpp::Function message)
 {

  // Internal function
  // Purpose: forms a map 
  // Return:  list object 


  std::string token, line;


  // open vcf marker text  file
  std::ifstream fileIN(fname.c_str());

  if(!fileIN.good()) {
      message("ERROR: Vcf file could not be opened with filename  " , fname , "\n" );
      return 0;
  }

  // check that vcf file is a vcf file
    getline(fileIN, line );  // read first line of file
    std::istringstream streamA(line);
    streamA >> token;  // first field of first line


  if (token.rfind("##fileformat=VC", 0) != 0) {
     message("ERROR: This does not appear to be a vcf file because the first row is not beginning with  ##fileformat=VCF... ");
     return 0;
  }

  // form map vectors
  int nsnps =  shouldsnpbremoved.size() - sum( shouldsnpbremoved );


  Rcpp::CharacterVector chrm(nsnps);
  Rcpp::NumericVector pos(nsnps);
  Rcpp::CharacterVector snp(nsnps);


std::string::size_type sz;     // alias of size_t

    long row_count = 0;
    long indx = 0;
  while(getline(fileIN, line)){
     std::istringstream streamA(line);
     streamA >> token;  // first field 
     if (token.rfind("##", 0) != 0){
        // line is not preamble
        if (token.rfind("#CHROM", 0) != 0){
            // data line
            if ( shouldsnpbremoved[row_count] == 0)
            {
                chrm[indx] = token;
                streamA >> token;
                pos[indx] = std::stod (token, &sz);
                streamA >> token;
                snp[indx] = token;

                indx++;

            }



            row_count++;  // number of snps
        } // end if
     } // end if
  }  // end while getline


  // closing file
  fileIN.close();


   // return dim of M (not Mt)
   return  Rcpp::DataFrame::create(Rcpp::Named("snp")=snp, Rcpp::Named("chrm")=chrm, Rcpp::Named("pos")=pos);
   


  }


 Rcpp::IntegerVector dimOfFile(     std::string fname,
                                  Rcpp::Function message)
 {

  // Internal function
  // Purpose: check that vcf has correct first row and 
  //          calculate the dimension of M (not Mt) 
  // Return:  array[2] with col, row 


  std::string token, line;
  Rcpp::IntegerVector dim_of_M(2) ;  // row x col of M matrix


  // open vcf marker text  file
  std::ifstream fileIN(fname.c_str());

  if(!fileIN.good()) {
      message("ERROR: Vcf file could not be opened with filename  " , fname , "\n" );
      return 0;
  }
 
  message(" Checking that vcf file is of correct format.");
  // check that vcf file is a vcf file
    getline(fileIN, line );  // read first line of file
    std::istringstream streamA(line);
    streamA >> token;  // first field of first line


  if (token.rfind("##fileformat=VC", 0) != 0) {
     message("ERROR: This does not appear to be a vcf file because the first row is not beginning with  ##fileformat=VCF... ");
     return 0;
  }

  // Determine number of individuals and snps 
    long row_count = 0;
    long col_count = 0;

  message(" Determining the number of individuals and snp contained in the file. ");
  while(getline(fileIN, line)){
     std::istringstream streamA(line);
     streamA >> token;  // first field 
     if (token.rfind("##", 0) != 0){
        // line is not preamble
        if (token.rfind("#CHROM", 0) != 0){
            // data line
            row_count++;  // number of snps
            if (row_count == 1){
               while(streamA >> token){
                 col_count++;  // number of individuals
               } // end while
            } // if row_count
        } // end if
     } // end if
  }  // end while getline

   col_count = col_count - 9 + 1 ; // removing cols that do not contain genogtype data but adjusting for
                                 // counts begin from 1, not 0. 
   //Rcpp::Rcout << "number rows " << row_count << " number of cols " << col_count << std::endl;
   dim_of_M[0] = col_count;  // number of individuals
   dim_of_M[1] = row_count;  // number of columns

   message("   Number of individuals: ", col_count);
   message("   Number of snp:         ", row_count);


  // closing file
  fileIN.close();


   // return dim of M (not Mt)
   return dim_of_M;

  }

 int  createBINfiles_outofmemory_M(std::string fnamebinMt,  std::string fnamebinM,  double  max_memory_in_Gbytes, 
Rcpp::IntegerVector dim_of_M,  Rcpp::Function message, bool quiet)
{
// Purpose:  create the M (transpose of Mt) file of the recoded vcf file. 
//           This is done when there is not enough memory to hold M in  memory.
// Return:   boolean vector of whether the snp was removed due to being monomorphic, of low MAF, or multiallelic. 


//  std::ifstream fileIN(fnamebinMt.c_str(), std::ios::in  );
 std::ifstream fileIN(fnamebinMt.c_str(), std::ios::in | std::ios::binary  );
  if(!fileIN.good()) {
      message("ERROR: Temporary binary file could not be opened with filename  " , fnamebinMt.c_str()  , "\n" );
      return 0;
  }


// std::ofstream fileOUT(fnamebinM.c_str(), std::ios::out  );
 std::ofstream fileOUT(fnamebinM.c_str(), std::ios::out | std::ios::binary );
  if(!fileOUT.good()) {
      message("ERROR: Temporary binary file could not be cretaed with filename  " , fnamebinM.c_str()  , "\n" );
      return 0;
  }



std::string token, line;

double
   max_mem_in_bytes  =  max_memory_in_Gbytes * 1000000000;



  //  Block approach needed due to lack of memory

   // Calculate number of rows that can be read into available memory
    //long n_of_rows_to_be_read = max_mem_in_bytes * 1.0 / (1.5 * dim_of_M[1] * (bits_in_int/8.0)); 
    long n_of_rows_to_be_read = max_mem_in_bytes * 1.0 / (1.2 * dim_of_M[1] * sizeof (char) ); 
    // Calculate number of blocks needed
    int n_blocks = dim_of_M[0]/n_of_rows_to_be_read;
    if (dim_of_M[0] % n_of_rows_to_be_read != 0)
         n_blocks++;
 
    message(" Not enough memory to perform transpose.");
    message(" Marker data being divided into blocks for out-of-memory transpose. ");
    message(" Number of blocks is ", n_blocks );

    // Block read and transpose - requires n_blocks passes through the 
    // binary Mt file which could be slow if file is large and memory low
    for(int b=0; b < n_blocks; b++){
              message( " Processing block " , b+1 , " of a total number of blocks of " , n_blocks );


         long
              start_val = b * n_of_rows_to_be_read,
              end_val   = (b+1) * n_of_rows_to_be_read;

         if (end_val > dim_of_M[0])
              end_val = dim_of_M[0];

         long nrows = end_val - start_val ;


          long   numelem = dim_of_M[1] * nrows;

          char * line ;
          line = new char[ dim_of_M[0] ];  // entire col of data to fit in here

          char * bufferT;
          bufferT = new char[numelem];  // contains transposed block of data


         //    long counter = 0;
         for(long coli=0; coli<dim_of_M[1]; coli++){
          // read a line of data from binary Mt file
          // read block into buffer
          fileIN.read( line, dim_of_M[0] );

         

              long rowi=0 ;
              for(long ii=start_val; ii < end_val ; ii++){
                bufferT[ rowi * dim_of_M[1] + coli]  =  line[  ii ]  ;
                rowi++;
              }
          }  // end for coli



       // write transposed block to binary file
       fileOUT.write( (char *) &bufferT[0], dim_of_M[1] * nrows  * sizeof(char));


      // return to the beginning of the input file
      fileIN.clear();
      fileIN.seekg(0, std::ios::beg);

      delete [] bufferT;
      bufferT = NULL;

      delete [] line;
      line = NULL;

     }     // end for blocks
     fileIN.close();
     fileOUT.close();



     return 1;   // successful completion ... 



}



 Rcpp::IntegerVector   createBINfiles_outofmemory_Mt(std::string fname, std::string fnamebinMt,  
                                                  Rcpp::IntegerVector dim_of_M,  Rcpp::Function message, bool quiet)
 {
  // Purpose:  create the Mt file of the recoded vcf file. 
  //           This is done when there is not enough memory to hold M in  memory.
  // Return:   boolean vector of whether the snp was removed due to being monomorphic, of low MAF, or multiallelic. 

  std::string token, line;



 // open output bin file that is to hold  Mt vcf genotype data
  std::ofstream fileOUT(fnamebinMt.c_str(), std::ios::out | std::ios::binary );
 // std::ofstream fileOUT(fnamebinMt.c_str(), std::ios::out  );
  if(!fileOUT.good()) {
      message("ERROR: Temporary binary file could not be opened with filename  " , fname , "\n" );
      return 0;
  }




  if (!quiet ){
      message("");
      message(" Reading vcf File  ");
      message("");
  }




  // open vcf marker text  file
  std::ifstream fileIN(fname.c_str());
  if(!fileIN.good()) {
      message("ERROR: Vcf file could not be opened with filename  " , fname , "\n" );
      return 0;
  }

  // skip over preamble 
  bool go = true;
  while(go){
     getline(fileIN, line);
     std::istringstream streamA(line);
     streamA >> token;  // first field of first line  
     if (token.rfind("##", 0) != 0) 
         go = false;  // finish passing over preamble
  }



  // check that this vcf file is a genotype file by looking for the FORMAT column name in the 9th column of the header row
  //getline(fileIN, line);
  std::istringstream streamB(line);
  for(int i=0; i<9; i++){
     streamB >> token; 
   }


  if (token.rfind("FORMAT", 0) != 0){
       message("ERROR: This does not appear to be a vcf file containing genotype calls because the 9th header column did not contain ");
       message("       the FORMAT field (which is required for genotype data). ");
       return 0;
   }



// Reading in the genotype calls and writing to Mt.bin line by line 

 // initializing input line 
std::vector<char> rowinfile(  dim_of_M[0] );
Rcpp::IntegerVector  shouldsnpbremoved(  dim_of_M[1]  ) ;



 int counter_col,  counter_row = 0,
    al=0; 

 bool geno=false;

 char previous;
 message(" Forming recoded marker file (this may take some time if the file is large).");

 while( getline(fileIN, line))
 {
   previous = line[0];  // first character in line
 
   counter_col = 0;  // genotype counter_col
    shouldsnpbremoved[counter_row] = 0;  // set to 0 initially

  // find starting value for i - need to skip over first 8 columns of data 
  // based on \t (tabs)
  long starting_i = 0;
  int number_of_tabs = 0;
  while(number_of_tabs < 9){
    if( line[starting_i] == '\t')
       number_of_tabs++;
    starting_i++;

  }

   for(int i=starting_i ; (unsigned)i < line.size(); i++)
   {

     if (previous == '\t' && line[i] == '0') {
        al = 0;
     }  else if (previous == '\t' && line[i] == '1') {
        al = 1;
     }  else if (previous == '\t' && ( line[i] == '2' ||  line[i] == '3' ||  line[i] == '4'||  line[i] == '5'||  line[i] == '6'||  line[i] == '7'||  line[i] == '8'||  line[i] == '9' )) {
        shouldsnpbremoved[counter_row] = 1;
       if (!quiet){
          message(" Snp number ", counter_row+1, " has been removed due to having more than two alleles. Allele ",  line[i], " has been found.");
        }

    }  else if (previous == '\t' && line[i] == '.') {
        al = -9;
        geno=true;
    }  else if ( al != -9 && shouldsnpbremoved[counter_row]==0  && (previous == '/' || previous == '|') && line[i] == '0' ){
      al = al + 0;
        geno=true;
    } else if ( al != -9 && shouldsnpbremoved[counter_row]==0  && (previous == '/' || previous == '|') && line[i] == '1' ){
      al = al + 1;
        geno=true;
    } else if  ( al != -9 && shouldsnpbremoved[counter_row]==0  && (previous == '/' || previous == '|') &&  ( line[i] == '2' ||  line[i] == '3' ||  line[i] == '4'||  line[i] == '5'||  line[i] == '6'||  line[i] == '7'||  line[i] == '8'||  line[i] == '9' )  ){
        shouldsnpbremoved[counter_row] = 1;
       if (!quiet){
          message(" Snp number ", counter_row+1, " has been removed due to having more than two alleles. Allele ",  line[i], " has been found.");
        }
   }  else if (  al != -9 && shouldsnpbremoved[counter_row]==0  && (previous == '/' || previous == '|') && line[i] == '.' ){
      al = -9;
        geno=true;
   } 
   if ( shouldsnpbremoved[counter_row] ){
      geno = false;
  }   



  if (geno ){
     // check that counter_row is not greater than 
     if (counter_col ==  dim_of_M[0]  ){
      message("\n");
      message("Error:  a problem has occurred when reading in the data from the vcf file. ");
      message("        The error has occurred for snp number " , counter_row+1, ".");
      message("        This snp may contain more than ", dim_of_M[0] , " columns of data or the snp may contain tabs at the end of the line. \n ");
      message(" ReadVCF has terminated with errors");
      message("\n");
      message("\n");
      return 0 ;
     }


   geno = false;
   if (al == 0)
      rowinfile[counter_col] = '0';
   if (al == 1 || al == -9)
      rowinfile[counter_col] = '1';
   if (al == 2)
      rowinfile[counter_col] = '2';


   counter_col++; 
  }

  previous = line[i];

  }  // end for



 if (counter_col != dim_of_M[0] &&  shouldsnpbremoved[counter_row]==0  ){
      message("\n");
      message("Error:  The vcf file contains an unequal number of snp (or columns) per individual (or row).  ");
      message("        The error has occurred for snp number " , counter_row+1 , ". It contains " , counter_col , " columns of data.");
      message("        It should contain " , dim_of_M[0] , " columns of data. \n");
      message(" ReadVCF has terminated with errors");
      message("\n");
      message("\n");
      return 0 ;
}
 
   if ( shouldsnpbremoved[counter_row]==0 ){
      // writing vector to binary file
       fileOUT.write( (char *) &rowinfile[0], rowinfile.size() * sizeof(char));
   }
   counter_row++;

} // end while getline

fileOUT.close();
fileIN.close();


  return shouldsnpbremoved;


}



 Rcpp::IntegerVector   createBINfiles_withinmemory(std::string fname, std::string fnamebinMt,  std::string fnamebinM,
                                                  Rcpp::IntegerVector dim_of_M,  Rcpp::Function message, bool quiet)
 {
  // Purpose:  create the Mt and M binary files of the recoded vcf file. 
  //           This is done within memory.  
  // Return:   boolean vector of whether the snp was removed due to being monomorphic, of low MAF, or multiallelic. 

  std::string token, line;




 // open output bin file that is to hold  Mt vcf genotype data
  std::ofstream fileOUT(fnamebinMt.c_str(), std::ios::out | std::ios::binary );
 // std::ofstream fileOUT(fnamebinMt.c_str(), std::ios::out  );
  if(!fileOUT.good()) {
      message("ERROR: Temporary binary file could not be opened with filename  " , fname , "\n" );
      return 0;
  }

 if(!quiet){
      message("");
      message("Name of temporary Mt file: ", fnamebinMt.c_str(), "\n");
      message("");
      message(" Reading vcf File  ");
      message("");
  }

    // genotypesTrans will hold recoded genotypes (0,1,2) 
    //Create an array of pointers that points to more arrays
    // Note - if some SNP are removed, then this array will contain non-assigned
    //        values.
    char** genotypesTrans = new char*[  dim_of_M[0] ];
    for (int i = 0; i <  dim_of_M[0] ; ++i) {
        genotypesTrans[i] = new char[  dim_of_M[1] ];
    }




  // open vcf marker text  file
  std::ifstream fileIN(fname.c_str());
  if(!fileIN.good()) {
      message("ERROR: Vcf file could not be opened with filename  " , fname , "\n" );
      return 0;
  }

  // skip over preamble 
  bool go = true;
  while(go){
     getline(fileIN, line);
     std::istringstream streamA(line);
     streamA >> token;  // first field of first line  
     if (token.rfind("##", 0) != 0) 
         go = false;  // finish passing over preamble
  }



  // check that this vcf file is a genotype file by looking for the FORMAT column name in the 9th column of the header row
  //getline(fileIN, line);
  std::istringstream streamB(line);
  for(int i=0; i<9; i++){
     streamB >> token; 
   }


  if (token.rfind("FORMAT", 0) != 0){
       message("ERROR: This does not appear to be a vcf file containing genotype calls because the 9th header column did not contain ");
       message("       the FORMAT field (which is required for genotype data). ");
       return 0;
   }



// Reading in the genotype calls and writing to Mt.bin line by line 
long snp_count = 0;  // row count

 // initializing input line 
std::vector<char> rowinfile(  dim_of_M[0] );
Rcpp::IntegerVector  shouldsnpbremoved(  dim_of_M[1]  ) ;



 int counter_col,  counter_row = 0,
    al=0;
 
 bool geno=false;
 char previous;
 message(" Forming recoded marker file (this may take some time if the file is large).");
 while( getline(fileIN, line))
 {
   previous = line[0];  // first character in line
 
   counter_col = 0;  // genotype counter_col
    shouldsnpbremoved[counter_row] = 0;  // set to 0 initially

  // find starting value for i - need to skip over first 8 columns of data 
  // based on \t (tabs)
  long starting_i = 0;
  int number_of_tabs = 0;
  while(number_of_tabs < 9){
    if( line[starting_i] == '\t')
       number_of_tabs++;
    starting_i++;

  }

   for(int i=starting_i ; (unsigned)i < line.size(); i++)
   {

     if (previous == '\t' && line[i] == '0') {
        al = 0;
     }  else if (previous == '\t' && line[i] == '1') {
        al = 1;
     }  else if (previous == '\t' && ( line[i] == '2' ||  line[i] == '3' ||  line[i] == '4'||  line[i] == '5'||  line[i] == '6'||  line[i] == '7'||  line[i] == '8'||  line[i] == '9' )) {
        shouldsnpbremoved[counter_row] = 1;
        if (!quiet){
          message(" Snp number ", counter_row+1, " has been removed due to having more than two alleles. Allele ",  line[i], " has been found.");
        }
    }  else if (previous == '\t' && line[i] == '.') {
        al = -9;
        geno=true;
    }  else if ( al != -9 && shouldsnpbremoved[counter_row]==0  && (previous == '/' || previous == '|') && line[i] == '0' ){
      al = al + 0;
        geno=true;
    } else if ( al != -9 && shouldsnpbremoved[counter_row]==0  && (previous == '/' || previous == '|') && line[i] == '1' ){
      al = al + 1;
        geno=true;
    } else if  ( al != -9 && shouldsnpbremoved[counter_row]==0  && (previous == '/' || previous == '|') &&  ( line[i] == '2' ||  line[i] == '3' ||  line[i] == '4'||  line[i] == '5'||  line[i] == '6'||  line[i] == '7'||  line[i] == '8'||  line[i] == '9' )  ){
        shouldsnpbremoved[counter_row] = 1;
        if (!quiet){
          message(" Snp number ", counter_row+1, " has been removed due to having more than two alleles. Allele ",  line[i], " has been found.");
        }
   }  else if (  al != -9 && shouldsnpbremoved[counter_row]==0  && (previous == '/' || previous == '|') && line[i] == '.' ){
      al = -9;
        geno=true;
   } 
   if ( shouldsnpbremoved[counter_row] ){
      geno = false;
  }   



  if (geno ){
     // check that counter_row is not greater than 
     if (counter_col ==  dim_of_M[0]  ){
      message("\n");
      message("Error:  a problem has occurred when reading in the data from the vcf file. ");
      message("        The error has occurred for snp number " , counter_row+1, ".");
      message("        This snp may contain more than ", dim_of_M[0] , " columns of data or the snp may contain tabs at the end of the line. \n ");
      message(" ReadVCF has terminated with errors");
      message("\n");
      message("\n");
      return 0 ;
     }


   geno = false;
   if (al == 0)
      rowinfile[counter_col] = '0';
   if (al == 1 || al == -9)
      rowinfile[counter_col] = '1';
   if (al == 2)
      rowinfile[counter_col] = '2';


   // genotypes[counter_row][counter_col] = rowinfile[counter_col];
   genotypesTrans[counter_col][snp_count] = rowinfile[counter_col];
   counter_col++; 


  }

  previous = line[i];

  }  // end for



 if (counter_col != dim_of_M[0] &&  shouldsnpbremoved[counter_row]==0  ){
      message("\n");
      message("Error:  The vcf file contains an unequal number of snp (or columns) per individual (or row).  ");
      message("        The error has occurred for snp number " , counter_row+1 , ". It contains " , counter_col , " columns of data.");
      message("        It should contain " , dim_of_M[0] , " columns of data. \n");
      message(" ReadVCF has terminated with errors");
      message("\n");
      message("\n");
      return 0 ;
}

  
 
   if ( shouldsnpbremoved[counter_row]==0 ){
      // writing vector to binary file
       fileOUT.write( (char *) &rowinfile[0], rowinfile.size() * sizeof(char));
   }

   // check that snp is not being removed
   if ( shouldsnpbremoved[counter_row]==0 ){ 
     snp_count++;
   }

   counter_row++;


} // end while getline


fileOUT.close();
fileIN.close();

// WRite out transpose of Mt (M) to binary file

// open output bin file that is to hold  Mt vcf genotype data
 std::ofstream fileOUT_M(fnamebinM.c_str(), std::ios::out | std::ios::binary );
 // std::ofstream fileOUT_M(fnamebinM.c_str(), std::ios::out  );
 if(!fileOUT_M.good()) {
      message("ERROR: Temporary binary file could not be opened with filename  " , fname , "\n" );
      return 0;
  }


 message(" Forming transpose of recoded marker file.  ");
 for(int i=0; i< dim_of_M[0]; i++){
       //fileOUT_M.write( (char *) &genotypesTrans[i][0] , dim_of_M[1] * sizeof(char));
       fileOUT_M.write( (char *) &genotypesTrans[i][0] , snp_count * sizeof(char));
 }

fileOUT_M.close();





 //Free each sub-array
  for (int i = 0; i < (int) dim_of_M[0] ; ++i) {
        delete[] genotypesTrans[i];   
    }
 //Free the array of pointers
    delete[] genotypesTrans;

 

  return shouldsnpbremoved;



}




// [[Rcpp::export]]
Rcpp::List   create_vcf_BIN_rcpp(Rcpp::CharacterVector f_name, Rcpp::CharacterVector f_name_bin_M,
                           Rcpp::CharacterVector f_name_bin_Mt, double  max_memory_in_Gbytes,  bool quiet, Rcpp::Function message) 
{
  // Rcpp function to create binary file from vcf formatted input files
 




  double
     max_mem_in_bytes  =  max_memory_in_Gbytes * 1000000000;

  std::string token, line;

  std::string
       fname = Rcpp::as<std::string>(f_name),
       fnamebinM = Rcpp::as<std::string>(f_name_bin_M),
       fnamebinMt = Rcpp::as<std::string>(f_name_bin_Mt);

  Rcpp::IntegerVector dim_of_M(2) ;  // row x col of M matrix





  // Check if file exists and find its dimension
  dim_of_M = dimOfFile(fname, message);
  if(dim_of_M.length() == 1){
    // problem with file
    return 0;
  }


Rcpp::IntegerVector shouldsnpbremoved;

  // check if data can be held in memory as char matrix 
  // double mem_bytes = (1.2 * dim_of_M[1] * dim_of_M[0] * (sizeof(char) * CHAR_BIT) )/8.0;  
  double mem_bytes = 1.2 * dim_of_M[1] * dim_of_M[0] * sizeof(char) ;  
if ( max_mem_in_bytes  >  mem_bytes){ 
    shouldsnpbremoved  = createBINfiles_withinmemory(fname, fnamebinMt,  fnamebinM,  dim_of_M,  message, quiet);
    // adjust dim_of_M for removed snp
    if( sum(shouldsnpbremoved) > 0)
      dim_of_M[1] = shouldsnpbremoved.size() - sum(shouldsnpbremoved);


} else {
   // create Mt first when we dont have enough memory to hold all the genotype data
    shouldsnpbremoved  = createBINfiles_outofmemory_Mt(fname, fnamebinMt,  dim_of_M,  message, quiet);

    // adjust dim_of_M for removed snp
    if( sum(shouldsnpbremoved) > 0)
      dim_of_M[1] = shouldsnpbremoved.size() - sum(shouldsnpbremoved);




   // create M by taking a blocking approach to the M binary file data that has been recorded already
    createBINfiles_outofmemory_M(fnamebinMt, fnamebinM,  max_memory_in_Gbytes ,  dim_of_M,  message, quiet);


}

 

  // create map data frame for R



 Rcpp::DataFrame map = createMapDF(fname, shouldsnpbremoved,  message);



//--------------------------------------
// Summary of Genotype File
//--------------------------------------
message( "\n\n                    Summary of Marker File  " );
message( "                   ~~~~~~~~~~~~~~~~~~~~~~~~   " );
message(" File type:                    " , "vcf"  );
message(" Number of individuals:        "     , dim_of_M[0] );
message(" Total number of snp:          "     ,  shouldsnpbremoved.size()  );
message(" Final number of snp:          "     , dim_of_M[1]   );
message(" Number of snp removed:        "     ,  shouldsnpbremoved.size() -  dim_of_M[1]  );
message(" File size (gigabytes):        "  , mem_bytes/(1.2*1000000000) ); // 1.2 because that is the factor I am using to give me some breathing room 
message(" Available memory (gigabytes): " , max_memory_in_Gbytes  );
message("\n\n" );


return Rcpp::List::create(Rcpp::Named("map") = map , Rcpp::Named("dim_of_M") = dim_of_M);







}




