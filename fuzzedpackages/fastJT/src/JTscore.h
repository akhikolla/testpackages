#ifndef __JTSCORE_H__
#define __JTSCORE_H__

// class to store the results

using namespace std;

class JTscore{

private:
	double 		jStar;
	double 	nonEqCount;
	long    	gSnipID;
public:
	JTscore(){};
	JTscore(double JStar, double NonEqCount, long gSnipID);

	double 		getJStar();
	double 	getNonEqCount();
	long 		getGSnipID();
};

JTscore::JTscore(double JStar, double NonEqCount, long GSnipID)
{
	this->jStar 	= JStar;
	this->nonEqCount= NonEqCount;
	this->gSnipID   = GSnipID;
}

double JTscore::getJStar()
{
	return jStar;	
}

double JTscore::getNonEqCount()
{
	return nonEqCount;
}

long JTscore::getGSnipID()
{
	return gSnipID;
}

#endif 
