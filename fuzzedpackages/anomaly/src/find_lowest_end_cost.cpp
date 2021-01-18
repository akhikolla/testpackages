namespace anomalymv
{

double find_lowest_end_cost(double* segmentcosts, int jj, int p, int l)
{

	int ii = 0, kk = jj;
	double min = 200;

	for (ii = 0; ii < (l+1); ii++)
	{

		if (min > segmentcosts[kk])
		{

			min = segmentcosts[kk];
			
		}

		kk = kk + p;

	}

	return(min);

}

} // namespace anomalymv
