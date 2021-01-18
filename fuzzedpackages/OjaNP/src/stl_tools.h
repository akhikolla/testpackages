#include <iostream>
#include <list>
#include <set>

//#include <iterator> //DF

//using namespace std;   //von df eingefügt

template<class T> ostream& operator<<(ostream& O,const list<T>& S)
{
	typename list<T>::const_iterator i;

	O << '<';
	for(i=S.begin(); i!=S.end(); i++)
	{
		if(i!=S.begin())
			O << ',';
		O << *i;
	}
	O << '>';
	
	return O;
}

template<class T> ostream& operator<<(ostream& O,const vector<T>& S)
{
	typename vector<T>::const_iterator i;

	O << '(';
	for(i=S.begin(); i!=S.end(); i++)
	{
		if(i!=S.begin())
			O << ',';
		O << *i;
	}
	O << ')';
	
	return O;
}

template<class T> ostream& operator<<(ostream& O,const set<T>& S)
{
	typename set<T>::const_iterator i;

	O << '{';
	for(i=S.begin(); i!=S.end(); i++)
	{
		if(i!=S.begin())
			O << ',';
		O << *i;
	}
	O << '}';
	
	return O;
}
