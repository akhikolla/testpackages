
#include "ForDebug.h"

/*Define the namespace*/
namespace ROPTLIB{

	void ForDebug::Print(const char *name, const double *M, integer row, integer col, integer num)
	{
		OUTSTREAM << "=============" << name << "============" << std::endl;
		if (col == 1 && num == 1)
		{
			for (integer i = 0; i < row; i++)
				OUTSTREAM << M[i] << std::endl;
		}
		else
			if (num == 1)
			{
				for (integer j = 0; j < row; j++)
				{
					for (integer k = 0; k < col; k++)
					{
						OUTSTREAM << M[j + row * k] << "\t";
					}
					OUTSTREAM << std::endl;
				}
			}
			else
			{
				for (integer i = 0; i < num; i++)
				{
					OUTSTREAM << "(:, :, " << i << ")" << std::endl;
					for (integer j = 0; j < row; j++)
					{
						for (integer k = 0; k < col; k++)
						{
							OUTSTREAM << M[i * row * col + j + row * k] << "\t";
						}
						OUTSTREAM << std::endl;
					}
				}
			}
	}
} /*end of ROPTLIB namespace*/
