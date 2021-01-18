#include <taq_import.h>

void BubsortDelDouble(vector<int> &nums, vector<int> &numsOut)
{
	int a, b, t, size;
	size = nums.size();
	for(a=1; a<size; a++) {
		for(b=size-1; b>=a; b--) {
			if(nums[b-1] > nums[b]) {
				t=nums[b-1];
				nums[b-1]=nums[b];
				nums[b]=t;
			}
		}
	}
	int r, c, Ind;
	for(r = 0; r < (size-1); r++)
	{
		Ind = 1;
		for(c = (r + 1); c < size; c++)
		{
			if(nums[r] == nums[c])
			{
				Ind = 0;
				break;
			}
		}
		if(Ind) numsOut.push_back(nums[r]);
	}
	numsOut.push_back(nums[size-1]);
}
