// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif


const size_t bits_in_int = std::numeric_limits<int>::digits;



// [[Rcpp::export]]
void  createMt_ASCII_rcpp(Rcpp::CharacterVector f_name,
                          Rcpp::CharacterVector f_name_ascii,
                          Rcpp::CharacterVector  type,
                              double  max_memory_in_Gbytes,  std::vector <long> dims,
                              bool  quiet, Rcpp::Function message )
{

// read data from M.ascii that has already been created and transpose this file

std::string
   token,
   line;

std::ostringstream
      os;

long
 n_of_cols_to_be_read;


std::string
     ftype = Rcpp::as<std::string>(type),
     fname = Rcpp::as<std::string>(f_name),
     fnameascii = Rcpp::as<std::string>(f_name_ascii);




double
   max_mem_in_bytes  =  max_memory_in_Gbytes * 1000000000;


// Calculate number of columns that can be read in as a block with XGb of
// memory. 
   // how much memory will be needed to  store M, takes it transpose M.transpose, and 
   // store its answer in Mt (+ .5 for a buffer)
double mem_bytes = 3.5 * dims[0] * dims[1] * (bits_in_int/8);  // assumes a 64 bit system


 // open  files
std::ofstream fileOUT(fnameascii.c_str(), std::ios::out);
std::ifstream fileIN(fname.c_str());

//-----------------------------------------------------------------------
//  Two situations
//   1.  memory X is sufficient to read all data into memory and transpose
//   2.  memory X is insufficient to read all data into memory. 
//------------------------------------------------------------------------




if(mem_bytes < max_mem_in_bytes){
     // Situation 1
     //-------------

     // open ASCII file and check for its existence. 
     if(!fileIN.good()) {
        os << "\n\nERROR: Could not open  " << fname << "\n\n" << std::endl;
        Rcpp::stop(os.str() );

     }

    // create matrix structure to hold genotype data
    Eigen::MatrixXi
          M(dims[0], dims[1]) ;


    // reset position in data file
    fileIN.clear();
    fileIN.seekg(0, std::ios::beg);

   // read values into matrix
   long rowi=0;
   while(getline(fileIN, line ))
   {

       // read a line of data from ASCII file
       for(long coli=0; coli < dims[1]; coli++){
           M(rowi, coli)  = line[coli] - '0'; // trick to removes ASCII character offset for numbers
       }
       rowi++;
   }  // end while getline

  // take transpose of matrix M
  Eigen::MatrixXi Mt = M.transpose();

  // write out contents fo Mt to file (no spaces)a
 std::vector<char> rowinfile( Mt.cols() );
 for(long i=0; i < Mt.cols() ; i++)
    rowinfile[i] = '0';



  for(long rowi=0; rowi<Mt.rows(); rowi++){
     for(long coli=0; coli<Mt.cols(); coli++){
         // fileOUT << Mt(rowi, coli);
        rowinfile[coli] =  Mt(rowi, coli) + '0'; // forming string row before writing to file
     }
 for(long ii=0; ii< Mt.cols() ; ii++){
     fileOUT << rowinfile[ii];
  }

     fileOUT << "\n";
  }

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
         Eigen::MatrixXi
              M(dims[0], ncols );


         // open ASCII file and check for its existence. 
         if(!fileIN.good()) {
              os << "ERROR: Could not open  " << fname << std::endl;
              Rcpp::stop(os.str() );
         }
     //    long counter = 0;
         if (!quiet ) {
              message("\n\n");
         }
         for(long rowi=0; rowi<dims[0]; rowi++){
               // read a line of data from ASCII file
              getline(fileIN, line);
              std::istringstream streamA(line);

             long coli=0 ;
             for(long ii=start_val; ii < end_val ; ii++){
                M(rowi, coli)  = line[ii] - '0'; // trick to removes ASCII character offset for numbers
                coli++;
             }
        } // end for(rowi=0; rowi<dims[0]; rowi++)
       // transpose M


       Eigen::MatrixXi Mt = M.transpose();


      // write out contents fo Mt to file (no spaces)a
std::vector<char> rowinfile( Mt.cols() ) ;


 for(long i=0; i < Mt.cols(); i++)
    rowinfile[i] = '0';


      for(long rowi=0; rowi<Mt.rows(); rowi++){
         for(long coli=0; coli<Mt.cols(); coli++){
             // fileOUT << Mt(rowi, coli);
            rowinfile[coli] =  Mt(rowi, coli) + '0'; // forming string row before writing to file
         }
         for(long ii=0; ii < Mt.cols(); ii++)
            fileOUT << rowinfile[ii]; // writing entire row of data
         fileOUT << "\n";
      }




      // return to the beginning of the input file
      fileIN.clear();
      fileIN.seekg(0, std::ios::beg);

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
message(" Reformatted ASCII file name:  " , fname  );
message(" Number of individuals:        "     , dims[0] );
if (ftype == "PLINK"  ){
// message(" Number of loci:           "  , (dims[1] -6)/2.0   );
message(" Number of loci:               "  , dims[1]   );
} else {
message(" Number of loci:               "  , dims[1] );
}
message( " File size (gigabytes):       "  , mem_bytes/1000000000 );
message(" Available memory (gigabytes): " , max_memory_in_Gbytes  );
message("\n\n" );


}  // end function 


