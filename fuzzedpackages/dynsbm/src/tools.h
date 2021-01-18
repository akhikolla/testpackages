/*
    This file is part of dynsbm.

    dysbm is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    dynsbm is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with dynsbm.  If not, see <http://www.gnu.org/licenses/>
 */
#ifndef DYNSBM_tOOLS
#define DYNSBM_tOOLS

namespace dynsbm{
  const double precision = 1e-10;
}

template<typename Ttype>
void allocate2D(Ttype**& ptr, int d1, int d2){
  ptr = new Ttype*[d1];
  ptr[0] = new Ttype[d1*d2];
  for (int p=0;p<d1*d2;p++) ptr[0][p] = Ttype(0);
  for (int i=1;i<d1;i++){
    ptr[i] = ptr[i-1] + d2;
  }
}
template<typename Ttype>
void deallocate2D(Ttype**& ptr, int d1, int d2){
  delete[] ptr[0];
  delete[] ptr;
}
template<typename Ttype>
void allocate3D(Ttype***& ptr, int d1, int d2, int d3){
  ptr = new Ttype**[d1];
  for (int i=0;i<d1;i++){
    ptr[i] = new Ttype*[d2];
    for (int j=0;j<d2;j++){
      ptr[i][j] = new Ttype[d3];
      for(int k=0;k<d3;k++){
	ptr[i][j][k] = Ttype(0);
      }
    }
  }
  /*
  for (int i=0;i<d1;i++){
    ptr[i] = new Ttype*[d2];
  }
  ptr[0][0] = new Ttype[d1*d2*d3];
  for (int p=0;p<d1*d2*d3;p++) ptr[0][0][p] = Ttype(0);
  int i=0, j=1;
  while(i<d1){
    if(j>0)
      ptr[i][j] = ptr[i][j-1]+d2;
    else
      ptr[i][j] = ptr[i-1][d2-1]+d2;
    j++;
    if(j==d2){
      j=0;i++;
    }
  }
  */
}
template<typename Ttype>
void deallocate3D(Ttype***& ptr, int d1, int d2, int d3){
  /*
    delete[] ptr[0][0];
    for (int i=0;i<d1;i++)
    delete[] ptr[i];
  */
  for (int i=0;i<d1;i++){
    for (int j=0;j<d2;j++){
      delete[] ptr[i][j];
    }
    delete[] ptr[i];
  }
  delete[] ptr;
}
template<typename Ttype>
void allocate4D(Ttype****& ptr, int d1, int d2, int d3, int d4){
  ptr = new Ttype***[d1];
  for (int i=0;i<d1;i++){
    ptr[i] = new Ttype**[d2];
    for (int j=0;j<d2;j++){
      ptr[i][j] = new Ttype*[d3];
      for(int k=0;k<d3;k++){
	ptr[i][j][k] = new Ttype[d4];
	 for(int l=0;l<d4;l++){
	   ptr[i][j][k][l] =Ttype(0);
	 }
      }
    }
  }
  /*
    ptr = new Ttype***[d1];
    for(int k=0;k<d1;k++) allocate3D<Ttype>(ptr[k], d2, d3, d4);
  */
}
template<typename Ttype>
void deallocate4D(Ttype****& ptr, int d1, int d2, int d3, int d4){
  for (int i=0;i<d1;i++){
    for (int j=0;j<d2;j++){
      for(int k=0;k<d3;k++){
	delete[] ptr[i][j][k];
      }
      delete[] ptr[i][j];
    }
    delete[] ptr[i];
  }
  delete[] ptr;
}
#endif
