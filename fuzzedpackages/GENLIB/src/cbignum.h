/*! \file cbignum.h
\brief Library Genlib: Classe et Enumation de base

	Declaration de toute les classes,enumeration pour tous les autres fichiers

	Ce fichier prend en charge toute les modificatiosn du code pour rendre le tous compatible entre win32 et unix

\author Sébastien Leclerc
\contributor Jean-François Lefebvre

*/

/* Calcul de de gros entier*/


#ifndef SBIGNUM
#define SBIGNUM

typedef unsigned char  UCHAR;
typedef unsigned char* PCHAR;

const short  CBignumBase = 256;
const unsigned char CBignumMax = 255;

///Structure qui represente un entier a precision multiple
class CBignum
{
	UCHAR* buffer;
	int  tailleBuffer;

	UCHAR* num;   //Null terminated string qui represente un gros nombre
	int  taillenum;

	UCHAR opShiftopMul(char direction);


	public:
	//FONCTION STANDARD

	CBignum();
	~CBignum();	

	void reset();
	void compact();
	int settaille(int Newtaille);
	int gettaille();
	UCHAR* getnum();
	void printstatus();

	//FONCTION DIMPORTATION EXPORTATION
	void ConvInputLong(unsigned int val); 
	unsigned int ConvOutputLong(); 

	void copy(const CBignum& Original); 

	//FONCTION MATHEMATIQUE

	void opMod(CBignum& num,CBignum& den); 

	//void opDiv(CBignum& num,CBignum& den);
	//void opAdd(CBignum& num,CBignum& den);
	//void opSub(CBignum& num,CBignum& den);
	//void opMul(CBignum& num,CBignum& den);

};


#endif


