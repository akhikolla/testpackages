#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include"Species.h"

class IndividualStatus {
	public:
	int initialSp;
	unsigned long id;
	double x, y;
	
	double creationTime, deathTime;
	std::list<double> growthTimes;

	IndividualStatus(int sp, unsigned long id, double x, double y, double ctime);
	
	void setGrowth(double time);
	void setDeath(double time);

	void addToHistory(Arena *ar);
};

/* INDIVIDUAL */

class Individual {
	private:
	static unsigned long id_MAX;
	Position p;
	const unsigned long id;
	int spnum;
	double baseD, baseG, baseR, Rad, SqRad; /* these are set at birth/growth */
	double D, G, R; /* these are the current rates */
	Species *species, *seedStage;
	Arena *arena;
	std::list<Individual*>::iterator ref;
	IndividualStatus *info;
	/* array of lists of neighbours by species */
	std::vector<std::list<Individual*> > affectingMeNeighbours;
	std::vector<std::list<Individual*> > affectedByMeNeighbours;
	void initNeighbours();
    void updateRates();

	public:
	Individual(Arena *ar, Species *sp, double x, double y);
	Individual(Arena *ar, Species *sp, Position p);
    Individual(Arena *ar, Species *sp, Position pos, unsigned long restored_id, double bTime);
	/* general action function */
	void act();

	/* GETS */
	double getTotalRate();
	int getSpeciesId();
	const unsigned long getId();
	Position getPosition();
	double getRadius();
	bool isPresent(Position p, double radius = 0);



	/* INTERACTIONS */
	/** adds a neighbour list and cross-adds yourself to everyone in that list */
	void addAffectingMeNeighbourList(std::list<Individual*> neighList);
	/** adds a neighbour to the list. Don't forget to add cross-reference to neighbour's list! */
	void addAffectingMeNeighbour(Individual *i);
	/** removes neighbour from list. Doesn't remove cross-reference from neighbour's list */
	void removeAffectingMeNeighbour(Individual *i);
	/** adds a neighbour list and cross-adds yourself to everyone in that list */
	void addAffectedByMeNeighbourList(std::list<Individual*> neighList);
	/** adds a neighbour to the list. Don't forget to add cross-reference to neighbour's list! */
	void addAffectedByMeNeighbour(Individual *i);
	/** removes neighbour from list. Doesn't remove cross-reference from neighbour's list */
	void removeAffectedByMeNeighbour(Individual *i);
	/** removes all neighbours */
	void clearNeighbours();
	bool noAffectingMeNeighbours(int i);
	
	~Individual();

	private:
	void setSpecies(Species *sp);
	void grow();
	void reproduce();
	void die();

};


#endif
