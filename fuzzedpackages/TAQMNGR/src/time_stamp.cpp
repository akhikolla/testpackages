#include <taq_aggregate.h>

int TimeStamp(int interval, vector<int> &tempo)
{
	int Open = 34200, Close = 57900, ComTime, nStamp = 0;
	
	ComTime = Open;
	while(ComTime <= Close)
	{
		tempo.push_back(ComTime);
		nStamp++;
		ComTime += interval;
	}
	if(tempo.back() != Close)
	{
		tempo.pop_back();
		tempo.push_back(Close);
	}
	return nStamp;
}
