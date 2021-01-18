//
//  debug.h
//  ccdr_proj
//
//  Created by Bryon Aragam on 3/30/15.
//  Copyright (c) 2014-2015 Bryon Aragam. All rights reserved.
//

#ifndef debug_h
#define debug_h

#include <iostream>

#include "defines.h"

//------------------------------------------------------------------------------/
//   HELPER FUNCTIONS FOR DEBUGGING
//------------------------------------------------------------------------------/

#ifdef _DEBUG_ON_
    // print only upper rxr principal submatrix
    std::string printToFile(SparseBlockMatrix& betas, int r){
        std::ostringstream matrix_out;
        matrix_out << "\n";
        
        r = std::min(betas.dim(), r); // if r > dimension then just print the whole thing
        int field_width = 10;
        
        for(int i = 0; i < r; ++i){
            for(int j = 0; j < r; ++j){
                int found = -1;
                for(int k = 0; k < betas.rowsizes(j); ++k){
                    if(betas.row(j, k) == i){
                        found = k;
                        break;
                    }
                }
                
                if(found >= 0){
                    matrix_out << std::setw(field_width) << std::setprecision(4);
                    if(fabs(betas.value(j, found)) > ZERO_THRESH){
                        matrix_out << betas.value(j, found);
                    } else{
                        matrix_out << "*";
                    }
                } else{
                    matrix_out << std::setw(field_width) << 0;
                }
            }
            
            matrix_out << std::endl;
        }
        
        return matrix_out.str();
    }

    // print a "moving" window snapshot of the beta matrix centred at (row, col)
    std::string printToFile(SparseBlockMatrix& betas, int row, int col){
        std::ostringstream matrix_out;
        matrix_out << "\n";
        
        int r = 3;
        int field_width = 8;
        int maxRow = row + r, minRow = row - r;
        int maxCol = col + r, minCol = col - r;
        
        if(minRow < 0){
            minRow = 0;
            maxRow += r;
        }
        if(maxRow >= betas.dim()){
            maxRow = betas.dim() - 1;
            minRow -= r;
        }
        if(minCol < 0){
            minCol = 0;
            maxCol += r;
        }
        if(maxCol >= betas.dim()){
            maxCol = betas.dim() - 1;
            minCol -= r;
        }
        if(minRow < 0 || minCol < 0) matrix_out << "Issue printing betas matrix!" << std::endl;
        
        matrix_out << std::endl << std::setw(field_width) << "";
        for(int j = minCol; j <= maxCol; ++j){
            matrix_out << std::setw(field_width) << j;
        }
        
        matrix_out << std::endl;
        for(int i = minRow; i <= maxRow; ++i){
            matrix_out << std::setw(field_width) << i;
            for(int j = minCol; j <= maxCol; ++j){
                int found = -1;
                for(int k = 0; k < betas.rowsizes(j); ++k){
                    if(betas.row(j, k) == i){
                        found = k;
                        break;
                    }
                }
                
                if(found >= 0){
                    matrix_out << std::setw(field_width) << std::setprecision(4);
                    if(fabs(betas.value(j, found)) > ZERO_THRESH){
                        matrix_out << betas.value(j, found);
                    } else{
                        matrix_out << "*";
                    }
                } else{
                    matrix_out << std::setw(field_width) << 0;
                }
                
            }
            
            matrix_out << std::endl;
        }
        
        return matrix_out.str();
    }
#endif

#endif
    
