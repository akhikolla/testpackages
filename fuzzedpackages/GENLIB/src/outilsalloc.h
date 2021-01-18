/*! \file blockalloc.h
\brief Implémentation d'une table d'allocation par block

\author Sébastien Leclerc
\contributor Jean-Francois Lefebvre

*/ 
#include "base.h"
#include "hal.h"

#include <Rcpp.h>
#define R_NO_REMAP

#ifndef BLOCKALLOC
#define BLOCKALLOC


template<class P> 
class BlockAlloc
{
private:
	static const int NBBLOCKALLOUER;
	GestionMemoire memcheck; //Pour la gestion de memoire par bloc

	int m_taille; //Taille de chaque block
	int m_count; //Nombre d'élément restant dans le block courant;
	P* current; //Pointeur vers le prochain élément a retourné

	int m_used;

public:
	BlockAlloc(size_t taille) 
	{
		m_used=1;m_taille=taille;m_count=0;current=NULL;
	}

	BlockAlloc()			  
	{
		m_used=0;m_taille=0;m_count=0;current=NULL;
	}
	
	inline void setTaille(size_t taille) 
	{
		//if (m_used==1)
		//	GENError("Utilisation de BlockAlloc invalide, l'objet est déjà initialisé");
		m_used=1;m_taille=taille;
	}
	
	inline P* Alloc()
	{
		try{
		if (!m_used){
			//GENError("Utilisation de BlockAlloc invalide, doit-être initialisé avant");
//			GENError("Invalid use of BlockAlloc, must be initialized first.");
			throw std::range_error("Invalid use of BlockAlloc, must be initialized first.");
		}
		if (m_count==0)
		{
			//Besoin d'un nouveau block
			const int tailleblock=m_taille*NBBLOCKALLOUER;
			current = (P*) memcheck.alloc(tailleblock,sizeof(P));
			m_count = m_taille;

			//RESET LE BLOCK
			for(int i=0;i<tailleblock;i++)
				current[i]=0;
		}
		P* tmp=current;
		current += m_taille;
		--m_count;

		return tmp;
 		} catch(std::exception &ex) {
 			forward_exception_to_r(ex);
 		}catch(...){
 			::Rf_error("c++ exception (unknown reason)"); 
 		} 
 		return 0;
	}

};


//CONSTANTE
template<class P> const int BlockAlloc<P>::NBBLOCKALLOUER=2000;

//PETITE PILE TRES RAPIDE...
template<class P, int SIZE> class SMPile
{
private:
	P array[SIZE];
	P* current;
public:
	inline SMPile() {current=array;}

	inline P top() {return *current;}
	inline P pop() {return *(current--);}
	inline void push(P elem) {*(++current)=elem;}
};

#endif



