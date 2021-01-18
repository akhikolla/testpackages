#ifndef __JTSCOREMP_H__
#define __JTSCOREMP_H__

using namespace std;

class JTscoreMP{

private:
	double 		jStar;
	double 	nonEqCount;
	long 		gSnipID;

public:
	JTscoreMP(){};
	JTscoreMP(double JStar, double NonEqCount, long gSnipID);

	double 		getJStar();
	double 	getNonEqCount();
	long 		getGSnipID();
};

JTscoreMP::JTscoreMP(double JStar, double NonEqCount, long GSnipID)
{
	this->jStar 	= JStar;
	this->nonEqCount= NonEqCount;
	this->gSnipID	= GSnipID;
}

double JTscoreMP::getJStar()
{
	return jStar;	
}

double JTscoreMP::getNonEqCount()
{
	return nonEqCount;
}

long JTscoreMP::getGSnipID()
{
	return gSnipID;
}

#endif 
