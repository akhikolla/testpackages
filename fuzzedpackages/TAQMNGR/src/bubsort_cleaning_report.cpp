#include <taq.h>

void bubsort_cleaning_report(vector<int> &date, vector<double> &total, vector<double> &notcorrected_delayed, vector<double> BrownGallo, int N)
{
	int a, b, c;
	double x, y, z;
	
	for(a = 1; a < N; a++)
	{
		for(b = N - 1; b >= a; b--)
		{
			if(date[b-1] > date[b])
			{
				c = date[b-1];
				x = total[b-1];
				y = notcorrected_delayed[b-1];
				z = BrownGallo[b-1];
				
				date[b-1] = date[b];
				total[b-1] = total[b];
				notcorrected_delayed[b-1] = notcorrected_delayed[b];
				BrownGallo[b-1] = BrownGallo[b];
				
				date[b] = c;
				total[b] = x;
				notcorrected_delayed[b] = y;
				BrownGallo[b] = z;
			}
		}
	}
}
				
	
	
