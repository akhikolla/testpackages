#include <taq.h>


void bubsort(vector<double> &nums, vector<int> &ordr, int size)
{
	int a, b, c;
	double t;
	for(a = 1; a < size; a++)
	{
		for(b = size - 1; b >= a; b--)
		{
			if(nums[b-1] > nums[b])
			{
				t = nums[b-1];
				c = ordr[b-1];
				nums[b-1] = nums[b];
				ordr[b-1] = ordr[b];
				nums[b] = t;
				ordr[b] = c;
			}
		}
	}
}
