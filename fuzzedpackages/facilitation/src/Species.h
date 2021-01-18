#ifndef SPECIES_H
#define SPECIES_H

#include"Arena.h"

class Species {
	private:
	int id;
	double D, G, R, dispersalRadius, Rad, maxStressEffect;
	int spnum, kernelType;
	double totalRate;

    Arena *arena;
	std::list<Individual*> population;
	Species *nextStage, *seedStage;
	/* array of interaction coeficients (affecting D,G,R) */
	double *interactionsD,*interactionsG,*interactionsR;

	public:
	Species(Arena *ar,int id, double *par);
	Species(Arena *ar,int id, double D, double G, double R, double dispersal, double Rad, double maxStressEffect, int dkernel);
	~Species();
	/* BASIC RUN ACTION */
	void act();

	/* INTERACTIONS */
	/* note: for the following functions, if radius is unspecified (=0), the radius used is the species own radius */
	bool isPresent(Position p, double radius = 0);
	std::list<Individual*> getPresent(Position p, double radius = 0);

	/* REPRODUCTION AND DEATH */
	std::list<Individual*>::iterator add(Individual *i);
	void remove(std::list<Individual*>::iterator i);
	void addIndividual(double x, double y);
	void addIndividual(Position p);	
	void disperseIndividual(double x, double y);
	void disperseIndividual(Position p);	
	Position dispersalKernel();


	/* SETS */
	void setNextStage(Species *st);
	void setSeedStage(Species *st);
	void setInteractionD(int s, double effect);
	void setInteractionG(int s, double effect);
	void setInteractionR(int s, double effect);

    void updateTotalRate(double change);

	/* GETS */
	double getTotalRate();
	double getG();
	double getR();
	double getD(Position p);
	double getRad();
	double getInteractionD(int species_id);
	double getInteractionG(int species_id);
	double getInteractionR(int species_id);
    bool   affectedBy(int species_id); /* finds if there is interaction */
	int getId();
	Species* getNextStage();
	Species* getSeedStage();

	int getAbundance();
};
#endif
