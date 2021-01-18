
namespace anomalymv
{

int cmpfunc_nosorting (const void * a, const void * b) 
{

	int out = 1;

	if (*(double*)a < *(double*)b )
	{
		out = -1;
	}

	return(out);

}

} // namespace anomalymv
