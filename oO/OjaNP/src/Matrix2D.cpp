#include "stdafx.h"
#include "Matrix2D.h"

#include <stdarg.h>
#include <math.h>
#include <limits>
#include <iostream>

#define DETERMINANT_FAILURE 0.01



using namespace std;

Matrix2D::Matrix2D()
{
  this->m = 0;
  this->n = 0;
}

Matrix2D::Matrix2D(int m, int n)
{
  this->m = m;
  this->n = n;
  values = new double*[this->m];
  for(int i = 0; i < m; i++)
  {
    values[i] = new double[n];
    for(int j = 0; j < n; j++)
      values[i][j] = 0;
  }
}
/*
Matrix2D::Matrix2D(double m, double n, ...)

{
  this->m = (int)m;
  this->n = (int)n;

  values = new double*[this->m];
  for(int i = 0; i < m; i++)
    values[i] = new double[this->n];

  va_list params;
  va_start(params, n);              // Aufruf an Initialisierungmakro

  for (int j = 0 ; j < this->n ; j++)
  {
    Vector v = va_arg(params, Vector); // Extrahiere einen Parameter
    for(int i = 0; i < this->m; i++)
      values[i][j] = v.getValue(i);
  }
  va_end(params);                        // Schliessmakro;
}

*/
Matrix2D::Matrix2D(const Matrix2D & reducedMat,int* validRows, int reducedM, int reducedN)
{
	this->m = reducedM;
	this->n = reducedN;
	
	this->values = new double *[this->m];     // speicher f�r zeilen anfordern
	for (int i = 0; i < m; i++)
	{
		this->values[i] = new double[this->n];  // speicher f�r spalten anfordern
	}
	int red_n_helper= reducedMat.n - reducedN;
	for (int r = 0; r < this->m; r++)
	{
		for (int c = 0; c < this->n; c++)
		{
			this->values[r][c] = reducedMat.values[validRows[r]][red_n_helper+c];    // daten kopieren
		}
	}
} 	

Matrix2D::~Matrix2D()
{
  for(int i = 0; i < m; i++)
    delete[] values[i];
  delete[] values;
}

bool Matrix2D::setRow(int i, Vector& v)
{
  if(v.getSize()!= n)
  {
    //cout << "Matrix2D::setRow - sizes do not match!" << endl;
    return false;
  }
  
  for(int j = 0; j < n; j++)
    values[i][j] = v.getValue(j);

  return true;
}

bool Matrix2D::setColumn(int j, Vector& v)
{
  if(v.getSize() != m)
  {
    //cout << "Matrix2D::setColumn - sizes do not match!" << endl;
    return false;
  }
  
  for(int i = 0; i < m; i++)
    values[i][j] = v.getValue(i);

  return true;
}

bool Matrix2D::setRow(int i, double* v)
{
  for(int j = 0; j < n; j++)
    values[i][j] = v[j];

  return true;
}

bool Matrix2D::setColumn(int j, double* v)
{
  for(int i = 0; i < m; i++)
    values[i][j] = v[i];

  return true;
}

bool Matrix2D::setColumn(int j, double p, Vector& v)
{  
  values[0][j] = p;
  for(int i = 1; i < m; i++)
    values[i][j] = v.getValue(i-1);

  return true;
}

bool Matrix2D::setValue(int i, int j, double v)
{
  if((i >= m) || (j >= n))
  {
    //cout << "Matrix2D::setValue - index too large!" << endl;
    return false;
  }

  values[i][j] = v;
  return true;
}


Vector Matrix2D::getRow(int i)
{
  if(i >= m)
  {
    //cout << "Matrix2D::Row - index too large!" << endl;
  }

  Vector v(n);
  for(int j = 0; j < n; j++)
    v.setValue(j, values[i][j]);
  return v;
}

Vector * Matrix2D::getColumn(int j)
{
  if(j >= n)
  {
    //cout << "Matrix2D::Column - index too large!" << endl;
  }

  Vector *v=new Vector(m);
  for(int i = 0; i < m; i++)
    v->setValue(i, values[i][j]);
  return v;
}

double Matrix2D::getValue(int i, int j)
{
  if((i >= m) || (j >= n))
  {
    //cout << "Matrix2D::Value - index too large!" << endl;
  }

  return values[i][j];
}

void Matrix2D::print() const
{
  //cout << "m = (" << endl;
  //for(int i = 0; i < m; i++)
  //{
    //cout << "  ";
    //for(int j  = 0; j < n-1; j++)
    //  cout << values[i][j] << ", ";
    //cout << values[i][n-1] << endl;
  //}
  //cout << ")" << endl;
}

Matrix2D & Matrix2D::transpose() const 
{
	Matrix2D * res = new Matrix2D(); // empty Matrix
	res->m = this->n;     // set Number of Rows to Number of Cols 
	res->n = this->m;     // set Number of Cols to Number of Rows
	res->values = new double *[res->m];   // allocate space for rows of Matrix
	for (int r = 0; r < res->m; r++)
	{
		res->values[r] = new double [res->n];   // allocate space for columns of Matrix
	}
	/* transpose data */ 
	for (int r = 0; r < res->m; r++)
	{
		for(int c=0; c < res->n; c++)
		{
			res->values[r][c] = this->values[c][r];   // switch by diagonal 
		}
	}
	return *res;
}

// generate a balanced Matrix with identical eigenvalues
Matrix2D & Matrix2D::balance() const
{
	if(this->n != this->m)
  	{
    //cout << "Matrix2D::determinant - not a quadratic matrix!" << endl;
    }
    Matrix2D D(this->n,this->n);
    for (int i = 0; i < this->n; i++)
    {
    	double sumA_ij = 0;
    	for (int j = 0; j< this->n; j++)
    	    sumA_ij += fabs(this->values[i][j]);
    	    
    	if (sumA_ij != 0)
    	   D.values[i][i] = double(1) / double (sumA_ij);  
    }
    Matrix2D & balanced = D.mult(*this);
    D.~Matrix2D();
     
	
	
	return balanced;
}

Matrix2D & Matrix2D::getACopy()const
{
	Matrix2D * mat_copy = new Matrix2D();
	
	mat_copy->n = this->n;
	mat_copy->m = this->m;
	
	mat_copy->values = new double*[this->m];
	for (int i = 0; i < this->m; i++)
	{
		mat_copy->values[i] = new double[this->n];
	}
	
	for (int r = 0; r < this->m; r++)
	{
		for (int c = 0; c < this->n; c++)
		{
			mat_copy->values[r][c] = this->values[r][c];
		}
	}
	
	
	return *mat_copy;
}

Matrix2D & Matrix2D::reduceToHessenbergForm(bool doBalancing) const
{
	Matrix2D & hessenberg = this->getACopy();
	
	double w;
	double c,s,h;
	double delta = numeric_limits<double>::epsilon();
	
	for (int j = 0; j < hessenberg.n-2; j++)
	{
		for (int i = j+2 ; i < hessenberg.n; i++)
		{
			if (hessenberg.values[i][j] != 0) 
			{
				if (hessenberg.values[j+1][j] < delta * fabs(hessenberg.values[i][j]))
				{
					w = - hessenberg.values[i][j];
					c = 0;
					s = 1;
				} else 
				{
					w = sqrt(hessenberg.values[j+1][j]*hessenberg.values[j+1][j] + hessenberg.values[i][j]*hessenberg.values[i][j]);
					w = (hessenberg.values[j+1][j] < 0)?(-w):(w);
					c = hessenberg.values[j+1][j] / w;
					s = -hessenberg.values[i][j] / w;
				}
				hessenberg.values[j+1][j] = w;
				hessenberg.values[i][j]   = 0;
				for (int k = j+1; k < hessenberg.n; k++)
				{
					h = c * hessenberg.values[j+1][k] - s * hessenberg.values[i][k];
					hessenberg.values[i][k] = s * hessenberg.values[j+1][k] + c * hessenberg.values[i][k];
					hessenberg.values[j+1][k] = h;
				}
				for (int k = 0; k < hessenberg.n; k++)
				{
					h = c * hessenberg.values[k][j+1] - s * hessenberg.values[k][i];
					hessenberg.values[k][i] = s * hessenberg.values[k][j+1] + c * hessenberg.values[k][i];
					hessenberg.values[k][j+1] = h;
				}
			}
		}
	}
	
		
	return hessenberg;
} 
  

void Matrix2D::printReducedToMapleFile (char * filename, char * matrix_variable_name , int * validRows, int reducedM, int reducedN)
{
	Matrix2D reduced(*this,validRows,reducedM,reducedN);
	reduced.printToMapleFile(filename,matrix_variable_name);
}

void Matrix2D::printToMapleFile(char * filename, char * matrix_variable_name)
{
	 /*FILE *outputDatei=fopen(filename,"a");	
     fprintf(outputDatei,matrix_variable_name); 
   	 
     fprintf(outputDatei," := matrix(%i,%i,[",this->m,this->n); 
   	 
   	for (int i = 0; i < this->m; i++){
	 fprintf(outputDatei,"["); 
	 for (int j = 0; j < this->n-1; j++)
	 {
	  fprintf(outputDatei,"%.15e,",this->values[i][j]);
	 }
	 if (i != this->m-1) 
	   fprintf(outputDatei,"%.15e],",this->values[i][this->n-1]);
	 else 
	   fprintf(outputDatei,"%.15e]]):;",this->values[i][this->n-1]);
	}
   	 
	fclose(outputDatei);*/
}



double Matrix2D::determinant(int* validRows, int reducedM, int reducedN)
{
  double det = 0;
  if(reducedN != reducedM)
  {
    //cout << "Matrix2D::determinant - not a quadratic matrix!" << endl;
    return det;
  }
  if(reducedN < 1)
  {
    //cout << "Matrix2D::determinant - empty matrix!" << endl;
    return det;
  }
  if(reducedN == 1)
  {
    det = values[validRows[0]][n-1];
    return det;
  }
  if(reducedN == 2)
  {
    det = values[validRows[0]][n-2] * values[validRows[1]][n-1] - values[validRows[0]][n-1] * values[validRows[1]][n-2];
    return det;
  }
  
  if(reducedN == 3)
  {
    det = values[validRows[0]][n-3] * values[validRows[1]][n-2] * values[validRows[2]][n-1]
        + values[validRows[0]][n-2] * values[validRows[1]][n-1] * values[validRows[2]][n-3]
        + values[validRows[0]][n-1] * values[validRows[1]][n-3] * values[validRows[2]][n-2]
        - values[validRows[2]][n-3] * values[validRows[1]][n-2] * values[validRows[0]][n-1]
        - values[validRows[2]][n-2] * values[validRows[1]][n-1] * values[validRows[0]][n-3]
        - values[validRows[2]][n-1] * values[validRows[1]][n-3] * values[validRows[0]][n-2];
    return det;
  }

  for(int i = 0; i < reducedM; i++)
  {
    int* newValidRows = new int[reducedM - 1];
    int k = 0;
    for(int j = 0; j < reducedM; j++)
    {
      if(i != j)
      {
        newValidRows[k] = validRows[j];
        k++;
      }
    }
    det += (pow(-1.0, 1.0 + i + 1.0) * values[validRows[i]][n - reducedN] * determinant(newValidRows, reducedM - 1, reducedN - 1));
    delete[] newValidRows;
  }

  return det;
}

Matrix2D & Matrix2D::mult(const Matrix2D & other_mat) const
{
        /*
	if (this->n != other_mat.m) 
	{
		cout << "Matrix2D::mult - the Matrices cannot be multiplied because the row number and the column number does not match" << endl;
	}
	*/
	Matrix2D *res = new Matrix2D();
	res->m = this->m;
	res->n = other_mat.n;
	res->values = new double *[res->m];   // allocate space for rows of Matrix
	for (int r = 0; r < res->m; r++)
	{
		res->values[r] = new double [res->n];   // allocate space for columns of Matrix
	}
	double temp_scalar;
	for (int r = 0; r < res->m; r++)
	{
		for (int c = 0; c < res->n; c++)
		{
			temp_scalar = 0;
			for (int k = 0; k < this->n; k++)
			{
				temp_scalar+=this->values[r][k] * other_mat.values[k][c];
			}
			res->values[r][c] = temp_scalar;
		}
	}	
	return *res; 
}

void Matrix2D::loadIdentity()
{
	if (this->m != this->n)
	{
		//cout << " Cannot load identity - matrix must be square " << endl;
		return;
	}
	for (int i = 0; i < this->n; i++)
	{
		for (int j = 0; j < this->m ; j++)
		{
			this->values[j][i] = 0;
		}
		this->values[i][i] = 1;
	}
}

double Matrix2D::determinant2(int * validRows,int reducedM, int reducedN)
{
	// m Zeilen , n Spalten 
	// n - reducedN           // der spalten Index der den Anfang der redurzierten Spalten angibt
	// ValidRows[k]           // Zeilen Indices der reduzierten Matrix
	if(reducedN != reducedM)
  	{
    	//cout << "Matrix2D::determinant - not a quadratic matrix!" << endl;
    	return 0;
  	}
	//DETERMINANT_FAILURE = 0.01;
    double sqrdet = 0;
    Matrix2D reduced(*this,validRows,reducedM,reducedN);
    Matrix2D & transp_reduced = reduced.transpose();
    Matrix2D & mtn = reduced.mult(transp_reduced); // mult_reduced_transporsed
 	
    // now mult_trans_norm symetric!    det(A*A^T) = det(A) * det(A^T)

    double lowst_dbl_greater_zero = std::numeric_limits<double>::epsilon();
    double sum;
    do
    {
        sum = 0;
        for (int i = 1; i < mtn.getNumberRows(); i++)
        {
            for(int j = 0; j < i; j++)
            {
                sum += (mtn.values[i][j]*mtn.values[i][j]);
            }
        }
        if (2 * sum < (DETERMINANT_FAILURE*DETERMINANT_FAILURE))
        {
            break;
        }
        for (int p = 0; p < mtn.getNumberRows()-1; p++)
        {
            for(int q = p+1; q < mtn.getNumberRows(); q++)
            {
                if (fabs(mtn.values[q][p]) >= (DETERMINANT_FAILURE*DETERMINANT_FAILURE))
                {
                     double delta = (mtn.values[q][q] - mtn.values[p][p]) / (2* mtn.values[q][p]);
                     double t = 1;
                     if (fabs(delta) > lowst_dbl_greater_zero)
                     {
                         double sign = (delta < 0)?(-1.0f):(1.0f);
                         t = 1 / (delta +  sign * sqrt(delta*delta + 1));
                     }
                         double c = 1 / sqrt(1+t*t);
                         double s = c*t;
                         double r = s/(1+c);
                         mtn.values[p][p] -= t * mtn.values[q][p];
                         mtn.values[q][q] += t * mtn.values[q][p];
                         mtn.values[q][p] = 0;
                         for (int j = 0; j < p; j++)
                         {
                             double g = mtn.values[q][j] + r * mtn.values[p][j];
                             double h = mtn.values[p][j] - r * mtn.values[q][j];
                             mtn.values[p][j] -= s * g;
                             mtn.values[q][j] += s * h;
                         }
                         for (int i = p+1; i < q; i++)
                         {
                             double g = mtn.values[q][i] + r * mtn.values[i][p];
                             double h = mtn.values[i][p] - r * mtn.values[q][i];
                             mtn.values[i][p] -= s * g;
                             mtn.values[q][i] += s * h;
                         }
                         for (int i = q+1; i < mtn.getNumberRows(); i++) 
                         {
                             double g = mtn.values[i][q] + r * mtn.values[i][p];
                             double h = mtn.values[i][p] - r * mtn.values[i][q];
                             mtn.values[i][p] -= s * g;
                             mtn.values[i][q] += s * h;
                         }
                } // end if
            } // endfor q
        } // endfor p
    } while(true);

    sqrdet = 1;
    
    for (int i = 0; i < mtn.getNumberRows(); i++)
    {
       sqrdet *= mtn.values[i][i];
    }
    
    return (sqrdet < 0)?(sqrt(-sqrdet)):(sqrt(sqrdet));
}

double Matrix2D::reduced_determinant_pivoting(int * validRows,int reducedM, int reducedN,int option)
{
 	Matrix2D reduced(*this,validRows,reducedM,reducedN);
 	double det = reduced.determinant_pivoting(false,option);
 	//reduced.~Matrix2D();
 	return det;
}

void Matrix2D::swap_rows(int row_1 , int row_2, int start_col_index , int end_col_index )
{
	for (int c = start_col_index; c < end_col_index; c++)
	{ 
		double temp = this->values[row_1][c];
		this->values[row_1][c] = this->values[row_2][c];
		this->values[row_2][c] = temp;
	}
}

void Matrix2D::swap_columns(int col_1, int col_2, int start_row_index, int end_row_index)
{
	for (int r = start_row_index; r < end_row_index; r++)
	{
		double temp = this->values[r][col_1];
		this->values[r][col_1] = this->values[r][col_2];
		this->values[r][col_2] = temp;
	}
}

void Matrix2D::do_pivoting( int r , int c)
//double Matrix2D::do_pivoting( int r , int c)
{
	for (int piv_row = r+1; (piv_row < this->m) ; piv_row++ )
	{
		for (int piv_col = c+1; piv_col < this->n; piv_col++)
		{
			this->values[piv_row][piv_col] = this->values[piv_row][piv_col] - this->values[piv_row][c]* this->values[r][piv_col] / this->values[r][c];
		}
	}
}

double Matrix2D::determinant_pivoting(bool copydata, int piv_option)
{
	/* option: 2: total pivoting 	  expensive	
	 *         1: column pivoting      
	 *         0: normal pivoting     fast but not very stable
	 * */
	 int pivot_row = 0;
	 int pivot_col = 0; 
	 int count_exchanges = 0;     
	 
	 Matrix2D * piv_Mat;
	 if (copydata)
	 	piv_Mat = & this->getACopy();//this->balance(); //this->getACopy();
	 else
	 	piv_Mat = this;
	 
	 
	 for (; (pivot_row < piv_Mat->m) && (pivot_col < piv_Mat->n); pivot_row++,pivot_col++)
	 {
	 	int next_piv_element_row = pivot_row;
	 	int next_piv_element_col = pivot_col;
	 	// find the next pivot element according to the pivot option
	 	switch (piv_option)  // s.o. 
	 	{
	 		case 0: // we take the A(pivot_row,pivot_col) element as our pivot element
	 		{
	 			// if this element is 0 than we have to search for an element != 0 
	 			// we do that in the A(...,pivot_col) column
	 			bool found_element = false;
	 			// next_piv_element_row is the first index of an element in this column neq 0
	 			for ( int r = pivot_row; (r < piv_Mat->m) && (!found_element) ; r++)
	 			{
	 				if (piv_Mat->values[r][pivot_col] != 0)
	 				{
	 					found_element = true;
	 					next_piv_element_row = r;
	 				}
	 			}
	 			break;
	 		}
	 		case 1: // we take the maximum element of the column A(...,pivot_col) as out next piv_element
	 		{
	 			double max_col = piv_Mat->values[pivot_row][pivot_col];
	 			for (int r = pivot_row; (r < piv_Mat->m); r++)
	 			{
	 				if (max_col < piv_Mat->values[r][pivot_col])
	 				{
	 					max_col = piv_Mat->values[r][pivot_col];
	 					next_piv_element_row = r;
	 				}
	 			}
	 			break;
	 		}
	 		case 2: // we take the maximum element of the submatrix (pivot_row ... m) X (pivot_col ... n)
	 		{
	 			double max = piv_Mat->values[pivot_row][pivot_col];
	 			for (int r = pivot_row; r < piv_Mat->m; r++)
	 			{
	 				for (int c = pivot_col; c < piv_Mat->n; c++)
	 				{
	 					if (max < piv_Mat->values[r][c])
	 					{
	 						max = piv_Mat->values[r][c];
	 						next_piv_element_row = r;
	 						next_piv_element_col = c;
	 					}
	 				}
	 			} 
	 			break;
	 		}
	 	}
	 	
	 	if (piv_Mat->values[next_piv_element_row][next_piv_element_col] == 0) 
	 	{
	 		// if the pivot element is 0 we can stop here and return 0 
	 		// because at least one of the diag. elments will be 0
	 		// so det(A) = det(diag(A)) = a11 * a22 * ... * ann = 0
	 		return 0;
	 	}	 	
	 	
	 	// if next_piv_element_row != pivot_row we have to swap these rows
	 	if (next_piv_element_row != pivot_row)
	 	{
	 		piv_Mat->swap_rows(pivot_row,next_piv_element_row,
	 		                pivot_col,piv_Mat->n); 
	 		count_exchanges++;
	 		// we have to swap only the elements pivot_col ... n-1
	 		// because the elements 0 ... pivot_col-1 will not effect 
	 		// our det(A) 		
	 	} 
	 	if (next_piv_element_col != pivot_col)
	 	{
	 		piv_Mat->swap_columns(pivot_col,next_piv_element_col,
	 						pivot_row,piv_Mat->m);
	 		count_exchanges++;
	 	}
	 	// AND FINALLY we do the pivoting for our submatrix 
	 	// (pivot_row+1 ... m) X (pivot_col+1 ... n)
	 	
	 	piv_Mat->do_pivoting(pivot_row,pivot_col); 
	 } 
	 double det = 1;
	 
	 for (int i = 0; i < piv_Mat->m; i++)
	 	det *= piv_Mat->values[i][i];
	 return (count_exchanges % 2 == 0)?(det):(-det);
}


int Matrix2D::getNumberColumns()
{
  return n;
}

int Matrix2D::getNumberRows()
{
  return m;
}

void Matrix2D::deleteColumn(int j, Matrix2D* matrix)
{
  for(int k = 0; k < m; k++)
  {
    for(int l = 0; l < n; l++)
    {
      if(l < j)
      {
        matrix->setValue(k, l, values[k][l]);
      }
      else if (l > j)
      {
        matrix->setValue(k, l-1, values[k][l]);
      }
    }
  }
}

void Matrix2D::deleteRow(int i, Matrix2D* matrix)
{
  for(int k = 0; k < m; k++)
  {
    for(int l = 0; l < n; l++)
    {
      if(k < i)
      {
        matrix->setValue(k, l, values[k][l]);
      }
      else if (k > i)
      {
        matrix->setValue(k-1, l, values[k][l]);
      }
    }
  }
}

void Matrix2D::deleteRowAndColumn(int i, int j, Matrix2D* matrix)
{
  for(int k = 0; k < m; k++)
  {
    for(int l = 0; l < n; l++)
    {
      if(k < i)
      {
        if(l < j)
          matrix->setValue(k, l, values[k][l]);
        else if (l > j)
          matrix->setValue(k, l-1, values[k][l]);
      }
      else if (k > i)
      {
        if(l < j)
          matrix->setValue(k-1, l, values[k][l]);
        else if (l > j)
          matrix->setValue(k-1, l-1, values[k][l]);
      }
    }
  }
}

double Matrix2D::getVariance()
{
  // expected value e
  double e = 0;
  double* length = new double[n];
  for(int j = 0; j < n; j++)
  {
    double sum = 0;
    for(int i = 0; i < m; i++)
    {
      sum += values[i][j] * values[i][j];
    }
    length[j] = sqrt(sum);
    e += length[j];
  }
  e /= n;

  // calculate variance
  double var = 0;
  for(int i = 0; i < n; i++)
  {
    var += (length[i] - e) * (length[i] - e);
  }
  var /= n;

  delete[] length;

  return var;
}
