#include "Functions.h"


namespace anomalymv
{
  
int cmpfunc_sorting (const void * a, const void * b) 
{

	int out = 1;

	if (((struct position_saving *)a)->saving < ((struct position_saving *)b)->saving )
	{
		out = -1;
	}

	return(out);
}

} // namespace anomalymv
