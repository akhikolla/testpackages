
#include"Individual.h"
#include"Random.h"

unsigned long Individual::id_MAX = 0;

Individual::Individual(Arena *ar, Species *sp, double x, double y):Individual(ar,sp,Position(x,y)){} 

/* This is now the secondary constructor */
Individual::Individual(Arena *ar, Species *sp, Position pos) : Individual( ar, sp, pos, id_MAX++, ar->getTotalTime() ){}

/* primary (bottom level) constructor */
Individual::Individual(Arena *ar, Species *sp, Position pos, unsigned long restored_id, double bTime) : id(restored_id),  arena(ar),
	affectingMeNeighbours(ar->getSpNum()+1), affectedByMeNeighbours(ar->getSpNum()+1) {

    if(id_MAX <= id) id_MAX = id+1;
	p = ar->boundaryCondition(pos), 
    D = G = R = 0; /* this sets totalRate to 0 so that the update works fine */
	spnum = arena->getSpNum();
	setSpecies(sp);
	info  = new IndividualStatus(sp->getId(),id,p.x,p.y,bTime);
	if(p.x==-1) die();
}



void	Individual::setSpecies(Species *sp) {
    D = G = R = 0; /* this sets totalRate to 0 so that the update works fine */

	species = sp;
	baseG = species->getG();
	baseR = species->getR();
	baseD = species->getD(p);
	Rad = species->getRad();
	SqRad = Rad*Rad;
	seedStage = species->getSeedStage();
	ref = species->add(this);

    updateRates();
	initNeighbours();
}

Individual::~Individual(){
	info->addToHistory(arena);
	clearNeighbours();
	delete(info);
}


int 		Individual::getSpeciesId(){return species->getId();}
Position 	Individual::getPosition(){return p;}
double 		Individual::getRadius(){return Rad;}
const unsigned long Individual::getId(){return id;}

void Individual::updateRates(){
	int sp, num;
	double effect, oldTotalRate = getTotalRate();
    D=baseD,G=baseG,R=baseR;
	for(sp = 1; sp <= spnum; sp++){
        num = affectingMeNeighbours[sp].size();/* note that effect is LINEAR on number of affecting neighbours */
		if(!affectingMeNeighbours[sp].empty()){
            if((effect = species->getInteractionD(sp)) != 0)
                D -= effect*num; /* note that D is a negative trait so positive effects decrease D */

            if((effect = species->getInteractionG(sp)) != 0)
                G += effect*num; 

            if((effect = species->getInteractionR(sp)) != 0)
                R += effect*num; 

		}
	}
	if(D < 0) D = 0;
	if(G < 0) G = 0;
	if(R < 0) R = 0;
    species->updateTotalRate(getTotalRate() - oldTotalRate);
}

double Individual::getTotalRate(){

	return D + G + R;
}

bool   Individual::isPresent(Position p2, double sqRadius){
	if(sqRadius == 0) sqRadius = SqRad;
	p2 -= p;
	if((p2.x)*(p2.x) + (p2.y)*(p2.y) < sqRadius) return true;
	else return false;
}

void   Individual::act(){
	double r = Random(getTotalRate());

	if(r < G) grow();
	else if (r < G+R) reproduce();
	else die();
}

void 	Individual::grow(){
	info->setGrowth(arena->getTotalTime());

	species->remove(this->ref);
	clearNeighbours();
    species->updateTotalRate(-getTotalRate());

	setSpecies(species->getNextStage()); /* note: setSpecies inits the neighbours */
}

void	Individual::reproduce(){
	species->disperseIndividual(p);
}

void 	Individual::die(){
	info->setDeath(arena->getTotalTime());

	species->remove(this->ref);
	clearNeighbours();
    species->updateTotalRate(-getTotalRate());

	delete(this);
}

void 	Individual::clearNeighbours(){
	int sp;
	std::list<Individual*>::iterator i;
	for(sp = 1; sp <= spnum; sp++){
		for(i=affectingMeNeighbours[sp].begin();i!=affectingMeNeighbours[sp].end();i = affectingMeNeighbours[sp].erase(i)){
			(*i)->removeAffectedByMeNeighbour(this);
		}

		for(i=affectedByMeNeighbours[sp].begin();i!=affectedByMeNeighbours[sp].end();i = affectedByMeNeighbours[sp].erase(i)){
			(*i)->removeAffectingMeNeighbour(this);
		}
	}
}

void Individual::initNeighbours(){
	int s;
    /* find what neighbours are affecting me */
	for(s=1;s<=spnum;s++){
		if(species->affectedBy(s) != 0){
			/* do not use radius in looking for affecting neighbours, 
			 * effect radius is the affecting neighbour's radius */
			addAffectingMeNeighbourList(arena->getPresent(s,p));
		}
	}

    /* find what neighbours are affected by me */
	/* this function automatically uses my radius */
	arena->addAffectedByMe(this);
}

void Individual::addAffectedByMeNeighbourList(std::list<Individual*> neighList){
	std::list<Individual*>::iterator i;

	for(i=neighList.begin();i!=neighList.end();i++){
		addAffectedByMeNeighbour(*i);
		/* makes sure that your neighbours adds you too */
		(*i)->addAffectingMeNeighbour(this);
	}
}

void Individual::addAffectingMeNeighbourList(std::list<Individual*> neighList){
	std::list<Individual*>::iterator i;

	for(i=neighList.begin();i!=neighList.end();i++){
		addAffectingMeNeighbour(*i);
		/* makes sure that your neighbours adds you too */
		(*i)->addAffectedByMeNeighbour(this);
	}
}

void 	Individual::addAffectedByMeNeighbour(Individual *i){
	int s = i->getSpeciesId();
	if(i==this){return;}
	affectedByMeNeighbours[s].push_back(i);
}

void 	Individual::addAffectingMeNeighbour(Individual *i){
	int s = i->getSpeciesId();
	if(i==this){return;}
	affectingMeNeighbours[s].push_back(i);
    updateRates(); /* update rates because affecting neighbour will affect me */
}

void 	Individual::removeAffectedByMeNeighbour(Individual *i){
	int s = i->getSpeciesId();
	affectedByMeNeighbours[s].remove(i);
}

void 	Individual::removeAffectingMeNeighbour(Individual *i){
	int s = i->getSpeciesId();
	affectingMeNeighbours[s].remove(i);
    updateRates(); /* update rates because affecting neighbour will stop affecting me */
}

bool 	Individual::noAffectingMeNeighbours(int i){
	return affectingMeNeighbours[i].empty();
}

IndividualStatus::IndividualStatus(int sp, unsigned long pid, double px, double py, double ctime):initialSp(sp),id(pid),x(px),y(py),creationTime(ctime),deathTime(-1){
	growthTimes = {};
}
void IndividualStatus::setGrowth(double time){ growthTimes.push_back(time); }
void IndividualStatus::setDeath(double time){ deathTime=time; }

void IndividualStatus::addToHistory(Arena *ar){
	std::list<double>::iterator i;
	double time1=creationTime,time2;
	int sp = initialSp;
	for(i = growthTimes.begin(); i!=growthTimes.end();i++){
		time2 = *i;
		ar->addToHistory(sp,id,x,y,time1,time2);
		time1 = time2;
		sp++;
	}
	ar->addToHistory(sp,id,x,y,time1,deathTime);
}
