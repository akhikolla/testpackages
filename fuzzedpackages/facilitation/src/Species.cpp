#include"Individual.h"
#include"Random.h"
#include<cstdio>

Species::Species(Arena *ar,int myid, double *par) : Species(ar,myid,par[0],par[1],par[2],par[3],par[4],par[5],par[6]){}

Species::Species(Arena *ar,int myid, double death, double growth, double rep=0, double dispersal=0, double radius=0, double maxEf=0, int dkernel=1)
    :id(myid),D(death),G(growth),R(rep),dispersalRadius(dispersal),Rad(radius),maxStressEffect(maxEf),kernelType(dkernel)
{
    int i;
    nextStage = NULL;
    seedStage = NULL;

    arena = ar;
    spnum = ar->getSpNum();

    interactionsD = (double*)malloc((spnum+1)*(sizeof(double)));
    interactionsG = (double*)malloc((spnum+1)*(sizeof(double)));
    interactionsR = (double*)malloc((spnum+1)*(sizeof(double)));
    interactionsD[0]=interactionsG[0]=interactionsR[0]=0; /* this is actually not used but let's set it to 0 just in case */
    for(i=1;i<=spnum;i++){
        interactionsD[i]=0;
        interactionsG[i]=0;
        interactionsR[i]=0;
    }
    totalRate = 0;
}

Species::~Species(){
    /* clear population */
    std::list<Individual*>::iterator i;

    for(i=population.begin();i!=population.end();i++){
        delete(*i);
    }

    free(interactionsD);
    free(interactionsG);
    free(interactionsR);
}

void Species::setInteractionD(int s, double effect){
    if(effect > D){
        Rcpp::warning("Interaction parameter set to be bigger than rate.");
    }
    interactionsD[s] = effect;
}

void Species::setInteractionG(int s, double effect){
    if(effect > G){
        Rcpp::warning("Interaction parameter set to be bigger than rate.");
    }
    interactionsG[s] = effect;
}

void Species::setInteractionR(int s, double effect){
    if(effect > R){
        Rcpp::warning("Interaction parameter set to be bigger than rate.");
    }
    interactionsR[s] = effect;
}

void Species::addIndividual(double x, double y){
    if(G > 0 && nextStage==NULL) {
        Rcpp::warning("Next stage set to NULL but G > 0. Check input data.");
        throw id;
    }
    if(R > 0 && seedStage==NULL) {
        Rcpp::warning("Seed stage set to NULL but R > 0. Check input data.");
        throw id;
    }
    /*Individual *i =*/ new Individual(arena,this,x,y);
}

void Species::addIndividual(Position p){
    addIndividual(p.x,p.y);
}

void Species::disperseIndividual(double x, double y){
    Position p(x,y);
    disperseIndividual(p);
}

void Species::disperseIndividual(Position p){
    if(kernelType==0){ /* Fully random on the arena */
        seedStage->addIndividual(Position(Random(arena->getWidth()),Random(arena->getHeight())));
    }
    else{
        seedStage->addIndividual(p + dispersalKernel());
    }
}

Position Species::dispersalKernel(){
    Position p;
    switch(kernelType){
        case 1: /* EXPONENTIAL */
        default:
            if(dispersalRadius <= 0) return Position(0,0);
            p = RandomDirection();
            return Exponential(1.0/dispersalRadius)*p;
    }
}

double Species::getTotalRate(){
    return totalRate;
}

bool Species::isPresent(Position p, double radius) {
    std::list<Individual*>::iterator i;

    for(i=population.begin();i!=population.end();i++){
        if((*i)->isPresent(p,radius*radius)) return true;
    }

    return false;
}

std::list<Individual*> Species::getPresent(Position p, double radius){
    std::list<Individual*>::iterator i;
    std::list<Individual*> list;

    for(i=population.begin();i!=population.end();i++){
        if((*i)->isPresent(p,radius*radius)) list.push_back(*i);
    }

    return list;
}

void Species::act(){
    std::list<Individual*>::iterator i;
    double r = Random(totalRate);

    for(i=population.begin();i!=population.end();i++){
        r -= (*i)->getTotalRate();
        if(r < 0) {
            (*i)->act();
            return;
        }
    }
    /* if the below code is executed, it's becase no individual was selected */
    Rcpp::warning ("No individual selected on Species::act.");
}

void Species::setNextStage(Species *st) {nextStage = st;}
void Species::setSeedStage(Species *st) {
    seedStage = st;
}

void Species::remove(std::list<Individual*>::iterator i){
    population.erase(i);
}

std::list<Individual*>::iterator Species::add(Individual *i){
    population.push_front(i);
    return population.begin();
}

Species* Species::getSeedStage() {return seedStage;}
Species* Species::getNextStage() {return nextStage;}
double Species::getG(){return G;}
double Species::getR(){return R;}
double Species::getRad(){return Rad;}
int Species::getId(){return id;}
double Species::getD(Position p){
    if(maxStressEffect == 0){
        return D;
    }
    else {
        return arena->getStressValue(p)*maxStressEffect;
    }
}

double Species::getInteractionD(int species_id){return interactionsD[species_id];}
double Species::getInteractionG(int species_id){return interactionsG[species_id];}
double Species::getInteractionR(int species_id){return interactionsR[species_id];}
bool Species::affectedBy(int s) {
    return (interactionsD[s] != 0 || interactionsG[s] != 0 || interactionsR[s] != 0);
}

int Species::getAbundance(){
    return population.size();
}


void Species::updateTotalRate(double change){
    totalRate+=change;

    /* recalculate rate if it's close to 0 */
    if (totalRate < 0.000000000000001){
        totalRate = 0;
        std::list<Individual*>::iterator i;

        for(i=population.begin();i!=population.end();i++){
            totalRate += (*i)->getTotalRate();
        }
    }
}

