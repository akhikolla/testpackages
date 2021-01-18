#include "base.h"

#ifndef GENOUTILSANAL
#define GENOUTILSANAL

#include "map"
typedef std::map<int, int> liste;

/**

	Cette classe permet de lier un individu à une série de nombres répétitifs non successifs
	Outil de la fonction compareGen

*/
class Tuple 
{
private:
	CIndSimul *ind;     ///< l'individu auquel il se rapporte
	liste tab;			///< map indiquant combien d'enfants ont telle ou telle valeur

public:
	Tuple();
	Tuple( const Tuple& );			///< constructeur de copie
	~Tuple(); 
	void addtab(int val);			///< incrémente d'un le nombre d'enfant ayant la valeur val
	const liste& gettab() const;	///< accesseur à l'attribut tab
	void setNoeud(CIndSimul* ref);	///< mutateur de l'attribut ind
	CIndSimul* getNoeud();			///< accesseur à l'attribut ind
	void clear(CIndSimul* ref=NULL);///< fonction de remise à zéro 
	Tuple& operator=(const Tuple& original);///< opérateur de copie

// Declaration de fonction Friend
friend int operator<(const Tuple& lhs, const Tuple& rhs);
friend int operator==(const Tuple& t2,const Tuple& t1);

}; 



int WINCDECL QSORTtuple(const void* e1,const void* e2);

int operator<(const Tuple& lhs, const Tuple& rhs);
int operator==(const Tuple& t2,const Tuple& t1);


#endif

