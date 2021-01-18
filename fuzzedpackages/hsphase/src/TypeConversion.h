// Copyright (C) 2014 Mohammad H. Ferdosi

#ifndef TYPECONVERSION_H
#define	TYPECONVERSION_H

#include <iostream>
#include <sstream>

using namespace std;

class typeConversion
{
public:
	typeConversion(const int initialIntiger, string &initialString);
	typeConversion(const unsigned int initialIntiger, string &initialString);
	typeConversion(const double initialIntiger, string &initialString);
	~typeConversion();
private:
	int Intiger;
	double Double;
	string String;

};

#endif	/* TYPECONVERSION_H */

