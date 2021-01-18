//readblock.h

#ifndef readblock_h_INCLUDED   // if x.h hasn't been included yet...
#define readblock_h_INCLUDED   //   #define this so the compiler knows it has been included


Eigen::MatrixXd  ReadBlock(std::string asciifname,
                           long start_row,
                           long numcols,
                           long numrows_in_block);

Eigen::MatrixXd  ReadBlockBin(std::string binfname,
                           long start_row,
                           long numcols,
                           long numrows_in_block);



#endif 
