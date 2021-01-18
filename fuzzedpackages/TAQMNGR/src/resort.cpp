#include <taq.h>


void resort(vector<double> &nums, vector<int> &ordr, double prezzo, int size)
{
	int t = 0, y, UpB = size-1;
	double x;
	
	while(ordr[t] != 1) t++;
	nums[t] = prezzo;
	ordr[t] = size + 1;
	for(int j = 0; j < size; j++) ordr[j] -= 1;
	
	while(t < UpB && nums[t] > nums[t+1])
	{
		x = nums[t];
		nums[t] = nums[t+1];
		nums[t+1] = x;
		y = ordr[t];
		ordr[t] = ordr[t+1];
		ordr[t+1] = y;
		t++;
	}
	
	while(t > 0 && nums[t] < nums[t-1])
	{
		x = nums[t];
		nums[t] = nums[t-1];
		nums[t-1] = x;
		y = ordr[t];
		ordr[t] = ordr[t-1];
		ordr[t-1] = y;
		t--;
	}
}
