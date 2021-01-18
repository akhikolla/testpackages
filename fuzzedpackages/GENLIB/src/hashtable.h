/*! \file hashtable.h
\brief Implémentation d'une table de hashage parametre

\author Claire Gires
\contributor Jean-François Lefebvre

*/ 


#ifndef HASHTABLE
#define HASHTABLE

#include "base.h"
#include "outils.h"

#ifdef NEG
	#undef NEG
#endif
extern "C" {
	#include "mpi.h"
	#include "mplogic.h"
}

const int VIDE = 0;
const double HASH_TAUX_REMPLISSAGE = 0.8;

// ********************************************************************
//
//			IMPLEMENTATION : HASHAGE DOUBLE
//
// ********************************************************************
template<class T=int>
class HashDouble
{
private:
	unsigned int m_taille;
	T* m_tab;
	char m_isError;
	unsigned int m_clef1;
	unsigned int m_clef2;
public:

	HashDouble(unsigned int nbElement, double ratio);
	~HashDouble();
	char Error();
	int Add(T& newElement);  ///Ajoute un nouvel element a la liste 
							//retourne GTRUE si le rajout a été fait (pas dupplicata)
							//retourne GFALSE si le rajout n'a pas été fait (dupplicata)
};

template<class T> inline HashDouble<T>::HashDouble(unsigned int nbElement, double ratio)
{

	m_taille = (unsigned int)(nbElement/ratio);

	if (m_taille%2 == 0)
		m_taille++;

	while(millerRabin(m_taille, 40))	
		m_taille+=2;

	//Fonction de hashage
	m_clef1=m_taille;
	m_clef2=m_taille-1;

	m_tab = new T[m_taille];
	if(m_tab==NULL)
	{
		m_isError=(m_tab==NULL);
		return;
	}
	//Remise à zero
	for(unsigned int i=0; i<m_taille; i++)
		m_tab[i] = VIDE;
	
}

template<class T> inline HashDouble<T>::~HashDouble()
{
	delete [] m_tab;
}

template<class T> inline char HashDouble<T>::Error() {return m_isError;}

template<class T> inline int HashDouble<T>::Add(T& newElement)
{
	unsigned int clef1, clef2,i;
	
	clef1 = newElement%m_clef1;
	
	if(m_tab[clef1]==VIDE)
		{
			m_tab[clef1]=newElement;
			return GTRUE;
		}

	if(m_tab[clef1]==newElement)
		return GFALSE;

	clef2 = (newElement%m_clef2)+1;
	
	i=0;
	const unsigned int limit = m_taille - clef2;		
	while(i<m_taille)
	{
		//On Avance
		if (clef1>=limit)
			clef1 -= limit; //Depassement
		else
			clef1 += clef2;			
		i++;

		if(m_tab[clef1]==VIDE)
		{
			m_tab[clef1]=newElement;
			return GTRUE;
		}

		if(m_tab[clef1]==newElement)
		return GFALSE;
	}
	m_isError=GTRUE; //SURCHARGE
	return GFALSE;
}


//DECLARATION DE SPECIALISATION MP_INT
template<> inline HashDouble<mp_int>::HashDouble(unsigned int nbElement, double ratio) 
{
	//Nombre d'élément dans la table de hashage
	m_taille = (unsigned int)(nbElement/ratio);

	if (m_taille%2 == 0)
		m_taille++;
	
	//Ajuste la taille du tableau au premier nombre impair suppérieur
	while(millerRabin(m_taille, 40))	
		m_taille+=2;
	
	//Fonction de hashage
	m_clef1=m_taille;
	m_clef2=m_taille-1;
	
	//Creation du tableau de MP_INT
	m_tab = new mp_int[m_taille];
	if(m_tab==NULL )
	{
		m_isError=(m_tab==NULL);
		return;
	}
	//Creation du tableau de valeur pour la table de hashage
	m_isError= (mp_init_array(m_tab, m_taille)!= MP_OKAY  );   	
}

template<> inline HashDouble<mp_int>::~HashDouble()
{
	mp_clear_array(m_tab, m_taille);
	delete [] m_tab;	
}

template<> inline int HashDouble<mp_int>::Add(mp_int& newElement)
{	
	unsigned int pos, interval;
	
	mp_mod_d(&newElement,m_clef1,&pos); //pos = newElement%m_clef1;
		
	//Est-ce que la case est vide
	if(mp_cmp_d(&m_tab[pos],VIDE)==0) //0 donne égal
	{
			mp_copy(&newElement,&m_tab[pos]);	
			return GTRUE;
	}
	
	//Est-ce un dupplicata direct
	if (mp_cmp(&newElement,&m_tab[pos])==0)
		return GFALSE;

	mp_mod_d(&newElement,m_clef2,&interval); //clef2 = newElement%m_clef2;
	interval+=1;

	unsigned int i=0;
	const unsigned int limit = m_taille - interval;		
	while(i++<m_taille) 
	{
		//Avance dans la table		
		if (pos>=limit)
			pos -= limit; //Depassement
		else
			pos += interval;			
	
		//Est-ce que la case est vide
		if(mp_cmp_d(&m_tab[pos],VIDE)==0) //0 donne égal
		{
				mp_copy(&newElement,&m_tab[pos]);		
				return GTRUE;
		}
		
		//Est-ce un dupplicata direct
		if (mp_cmp(&newElement,&m_tab[pos])==0) //0 donne égal
			return GFALSE;
	}
	m_isError=GTRUE;
	return GFALSE; //SURCHARGE...
}

// ********************************************************************
//
//			IMPLEMENTATION : HASHAGE SIMPLE
//
// ********************************************************************

template<class T=int>
class HashTable
{
private:
	unsigned int m_taille;

	//Classe pour le contenu...
	struct InterE
	{
		T elem;
		InterE* next;
	} *m_tab;  //Table de hashage...

public:
	HashTable();
	HashTable(unsigned int nbElement, double ratio=HASH_TAUX_REMPLISSAGE);
	~HashTable();

	void Initialise(unsigned int nbElement, double ratio=HASH_TAUX_REMPLISSAGE);
 	
	int add(T& newElement); ///Ajoute un nouvel element a la table de hashage
							//retourne GTRUE si le rajout a été fait (pas dupplicata)
							//retourne GFALSE si le rajout n'a pas été fait (dupplicata)
	int del(T& newElement);
							//retourne GTRUE si l'effacage est fait
							//retourne GFALSE si newElement n'as pas été trouvé
};

//FONCTION
template<class T> inline  HashTable<T>::HashTable() {m_taille=-1;m_tab=NULL;}
template<class T> inline  HashTable<T>::~HashTable() {if (m_tab) delete m_tab;}

template<class T> inline  HashTable<T>::HashTable(unsigned int nbElement, double ratio)
{Initialise(nbElement,ratio);}

template<class T> inline void HashTable<T>::Initialise(unsigned int nbElement, double ratio)
{
	m_taille = (unsigned int)(nbElement/ratio);

	if (m_taille%2 == 0)
		m_taille++;

	while(millerRabin(m_taille, 40))	
		m_taille+=2;
	
	//Creation du tableau
	m_tab = new InterE[m_taille];

	//REMISE A ZERO..(remet le contenue et Next a zero)
	for(unsigned int i=0; i<m_taille; i++)
	{
		m_tab[i].elem = VIDE;
		m_tab[i].next = NULL;
	}		
}

template<class T> inline int HashTable<T>::add(T& newElement)
{
	//Trouve la position du truc
	const unsigned int pos = (unsigned int)(newElement) % m_taille;
	
	//Est-ce qu'il est là....
	if (m_tab[pos].elem==VIDE)
	{
		//CASE VIDE DONC LE RAJOUTE...
		m_tab[pos].elem=newElement;
		return GTRUE;
	}
	else
		if (m_tab[pos].elem==newElement)
			return GFALSE; //Oui...Dupplicata
		else
		{
			//ici sa se complique....
			InterE *tmp=&m_tab[pos];
			while(tmp->next!=NULL)
			{
				tmp = tmp->next;
				if (tmp->elem==newElement)
					return GFALSE; //Oui...Dupplicata							
			}

			//PAS DANS LA LISTE CHAINE DONC INSERTION
			tmp->next	= new InterE;
			tmp			= tmp->next; //avance curseur
			tmp->elem	= newElement;
			tmp->next	= NULL;
			printf("\nOVERFLOW");
			return GTRUE; //Ajoute pas de dupplicata	
		}

}
template<class T> inline int HashTable<T>::del(T& newElement)
{

	//TROUVE LA POSITION DU TRUC
	const unsigned int pos = (unsigned int)(newElement) % m_taille;
	
	//Est-ce qu'il est là....
	if (m_tab[pos].elem==VIDE)
	{
		//CASE VIDE PEUT DONC PAS L'EFFACE
		return GFALSE;
	}
	else
		if (m_tab[pos].elem==newElement)
		{
			//BON, ON EFFACE
			if (m_tab[pos].next==NULL)			
				//Pas d'élément dans la liste chaine
				m_tab[pos].elem=VIDE;		
			else
			{ //Il y a un élément dans la liste chainé
				InterE* next=m_tab[pos].next;
				memcpy(&m_tab[pos],next,sizeof(InterE));
				printf("\nUNDERFLOW");
				delete next; //Efface le bloc initule
			}
			return GTRUE; //Element efface
		}
		else
		{
			//ici sa se complique....
			InterE *cur=m_tab[pos].next;
			InterE *prev=&m_tab[pos];
			while(cur!=NULL)
			{
				if (cur->elem==newElement)
				{
					//TROUVE, ON EFFACE
					prev->next=cur->next;
					delete cur;
					printf("\nUNDERFLOW");
					return GTRUE;					
				}
				
				//avamce
				prev	= cur;
				cur		= cur->next;
			}
			//RIEN TROUVE DONC PAS D'EFFACAGE			
			return GFALSE; //PAS TROUVE ON EFFACE PAS
		}
}



#endif

