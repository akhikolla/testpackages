// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif


const size_t bits_in_int = std::numeric_limits<int>::digits;



// [[Rcpp::export]]
void  createMt_BIN_rcpp(Rcpp::CharacterVector f_name_in,
                          Rcpp::CharacterVector f_name_out,
                          Rcpp::CharacterVector  type,
                              double  max_memory_in_Gbytes,  std::vector <long> dims,
                              bool  quiet, Rcpp::Function message )
{

// read data from M.bin that has already been created and transpose this file

std::string
   token;

std::ostringstream
      os;

long
 n_of_cols_to_be_read;


std::string
     ftype = Rcpp::as<std::string>(type),
     fname_in = Rcpp::as<std::string>(f_name_in),
     fname_out = Rcpp::as<std::string>(f_name_out);




double
   max_mem_in_bytes  =  max_memory_in_Gbytes * 1000000000;


// Calculate number of columns that can be read in as a block with XGb of
// memory. 
   // how much memory will be needed to  store M, takes it transpose M.transpose, and 
   // store its answer in Mt (+ .5 for a buffer)
double mem_bytes = 3.5 * dims[0] * dims[1] * (bits_in_int/8);  // assumes a 64 bit system


 // open  files
std::ifstream fileIN(fname_in.c_str(),   std::ios::in | std::ios::binary );
std::ofstream fileOUT(fname_out.c_str(), std::ios::out | std::ios::binary );

//-----------------------------------------------------------------------
//  Two situations
//   1.  memory X is sufficient to read all data into memory and transpose
//   2.  memory X is insufficient to read all data into memory. 
//------------------------------------------------------------------------




if(mem_bytes < max_mem_in_bytes){
     // Situation 1
     //-------------


     // open bin file and check for its existence. 
     if(!fileIN.good()) {
        os << "\n\nERROR: Could not open  " << fname_in << "\n\n" << std::endl;
        Rcpp::stop(os.str() );

     }

    // create matrix structure to hold genotype data
//    Eigen::MatrixXi
//          M(dims[0], dims[1]) ;


    long   numelem = dims[0] * dims[1];

    char * buffer ;
    buffer = new char[numelem];

    char * bufferT;
    bufferT = new char[numelem];


   // read block into buffer
   fileIN.read( buffer, numelem );


   // Take transpose
  #if defined(_OPENMP)
     #pragma omp for  
  #endif
 for(long ii=0; ii < dims[0]; ii++){
  for(long jj=0; jj < dims[1] ; jj++){
    bufferT[ jj * dims[0] + ii]  =  buffer[ ii*dims[1] + jj ]  ;  
  }
 }


 fileOUT.write( (char *) &bufferT[0], dims[0] * dims[1]  * sizeof(char));


 fileIN.close();
 fileOUT.close();
} else {
   //  Situation 2 
   //  Block approach needed due to lack of memory

   // if (!quiet ){
        message( " A block transpose is being performed due to lack of memory.  ");
        message( " Memory parameter availmemGb is set to " , max_memory_in_Gbytes , " gigabytes" );
        message( " If possible, increase availmemGb parameter. \n\n " );
   //  }

    // Calculate number of columns that can be read into available memory
    n_of_cols_to_be_read = max_mem_in_bytes * 1.0 / (3.5 * dims[0] * (bits_in_int/8.0)); //64 bit system
    // Calculate number of blocks needed
    long n_blocks = dims[1]/n_of_cols_to_be_read;
    if (dims[1] % n_of_cols_to_be_read != 0)
         n_blocks++;

  //  if (!quiet ) {
         message(" Number of blocks being processed is ", n_blocks, "\n\n");
   //  }
    // Block read and transpose - requires n_blocks passes through the 
    // ASCII input file which could be slow if file is large and memory low
    for(long b=0; b < n_blocks; b++){
     //    if (!quiet )
              message( " Processing block " , b+1 , " of a total number of blocks of " , n_blocks );


         long
              start_val = b * n_of_cols_to_be_read,
              end_val   = (b+1) * n_of_cols_to_be_read;

         if (end_val > dims[1])
              end_val = dims[1];

         long ncols = end_val - start_val ;


          long   numelem = dims[0] * ncols;

          char * line ;
          line = new char[ dims[1] ];  // entire row of data to fit in here

          char * bufferT;
          bufferT = new char[numelem];  // contains transposed block of data





         //    long counter = 0;
         for(long rowi=0; rowi<dims[0]; rowi++){
          // read a line of data from binary M file
          // read block into buffer
          fileIN.read( line, dims[1] );

         // open binary file and check for its existence. 
         if(!fileIN.good()) {
              os << "ERROR: Could not open  " << fname_in << std::endl;
              Rcpp::stop(os.str() );
         }

              long coli=0 ;
              for(long ii=start_val; ii < end_val ; ii++){
                bufferT[ coli * dims[0] + rowi]  =  line[  ii ]  ;
                coli++;
              }
          }  // end for rowi




       // write transposed block to binary file
       fileOUT.write( (char *) &bufferT[0], dims[0] * ncols  * sizeof(char));


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


}  // end if else situation 




//--------------------------------------
// Summary of Genotype File
//--------------------------------------

message( "\n\n                    Summary of Marker File  " );
message( "                   ~~~~~~~~~~~~~~~~~~~~~~~~   " );
message( " File type:                   " , type  );
message(" Reformatted ASCII file name:  " , fname_in  );
message(" Number of individuals:        "     , dims[0] );
if (ftype == "PLINK"  ){
// message(" Number of loci:           "  , (dims[1] -6)/2.0   );
message(" Number of loci:               "  , dims[1]   );
} else {
message(" Number of loci:               "  , dims[1] );
}
message( " File size (gigabytes):       "  , mem_bytes/(3.5*1000000000) ); // 3.5 because that is the factor I am using to give me some breathing room 
message(" Available memory (gigabytes): " , max_memory_in_Gbytes  );
message("\n\n" );


}  // end function 


