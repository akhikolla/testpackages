#include"Individual.h"
#include"Random.h"

Arena::Arena(int maxspid, double *parameters, double w, double h, int bc, double starttime) :maxsp(maxspid),width(w),height(h),bcond(bc) {
	int i;
    /* allocates vector of size n+1 so that we can ignore number 0 */
	species = (Species**)malloc((1+maxsp)*(sizeof(Species*)));
	ratesList = (double*)malloc((1+maxsp)*(sizeof(double)));

	for(i=1;i<=maxsp;i++){
		species[i] = new Species(this,i,parameters+FACILITATION_NUMPARAMETERS*(i-1));
	}

	totalTime = starttime;

	history = new History();
}

void Arena::createStructuredSpecies(int minId, int maxId) {
	int i;

	for(i=minId;i<maxId;i++){
		species[i]->setNextStage(species[i+1]);
		species[i]->setSeedStage(species[minId]);
	}
	// last stage doesn't have next stage
	species[i]->setSeedStage(species[minId]);
}

void Arena::createSimpleSpecies(int id){
	species[id]->setSeedStage(species[id]);
}

History * Arena::finalStatus(){
	int i;
	for(i=1;i<=maxsp;i++){
		delete(species[i]);
	}
	free(species);
	free(ratesList);
	return history;
}

void Arena::setInteractionsD(double *interactions){
	int i,j;
	for(i=1;i<=maxsp;i++){
		for(j=1;j<=maxsp;j++){
			species[i]->setInteractionD(j,interactions[maxsp*(i-1)+(j-1)]);
		}
	}
}

void Arena::setInteractionsG(double *interactions){
	int i,j;
	for(i=1;i<=maxsp;i++){
		for(j=1;j<=maxsp;j++){
			species[i]->setInteractionG(j,interactions[maxsp*(i-1)+(j-1)]);
		}
	}
}

void Arena::setInteractionsR(double *interactions){
	int i,j;
	for(i=1;i<=maxsp;i++){
		for(j=1;j<=maxsp;j++){
			species[i]->setInteractionR(j,interactions[maxsp*(i-1)+(j-1)]);
		}
	}
}

void Arena::populate(int *speciesinit){
	int i,j;

	for(i=1;i<=maxsp;i++){
		for(j=0;j<speciesinit[i-1];j++){
            species[i]->addIndividual(Random(width),Random(height));
		}
	}
}

double History::globalEndTime(){
    return std::max(*std::max_element(beginTime_list.begin(),beginTime_list.end()),*std::max_element(endTime_list.begin(),endTime_list.end()));
}

int History::size(){ return sp_list.size(); }

void Arena::populate(Rcpp::DataFrame init){
    int i;
    double t;
    std::vector<int> sp;
    std::vector<unsigned long> id;
    std::vector<double> x,y,beginTime,endTime;

    if(init.nrows()==0){
    }
    else {
        sp = Rcpp::as<std::vector<int> >(init["sp"]);
        id = Rcpp::as<std::vector<unsigned long> >(init["id"]);
        x = Rcpp::as<std::vector<double> >(init["x"]);
        y = Rcpp::as<std::vector<double> >(init["y"]);
        beginTime = Rcpp::as<std::vector<double> >(init["begintime"]);
        endTime = Rcpp::as<std::vector<double> >(init["endtime"]);

        for(i=0;i<init.nrows();i++){
            if(Rcpp::NumericVector::is_na(endTime[i])){
                new Individual(this,species[sp[i]],Position(x[i],y[i]),id[i],beginTime[i]);
            }
        }
        t = std::max(*std::max_element(beginTime.begin(),beginTime.end()),*std::max_element(endTime.begin(),endTime.end()));
        if(t>totalTime){
            totalTime = t;
        }
    }
    history->globalBeginTime = totalTime;
}

bool Arena::turn() {
    int i;
    double r, time;

    totalRate = 0;
    for(i=1;i<=maxsp;i++){
        ratesList[i] =  species[i]->getTotalRate();
        totalRate += ratesList[i];
    }

    if(totalRate < 0) {
        Rcpp::warning("#This simulation has reached an impossible state (totalRate < 0).");
        return false;
    }

    if(totalRate == 0) {
        Rcpp::warning("#This simulation has reached a stable state (totalRate = 0).");
        return false;
    }

    time = Exponential(totalRate);
    totalTime += time;

    /* select stage to act */
    r = Random(totalRate);
    for(i=1;i<=maxsp-1;i++){
        r -= ratesList[i];
        if(r < 0){
            break;
        }
    }
    species[i]->act();
    return true;

}

/*TODO: should this array be dynamically allocated? */
int* Arena::getAbundance(){
    int i;
    int *ab;
    ab = (int*)malloc(maxsp*sizeof(int));
    for(i=1;i<=maxsp;i++){
        ab[i] = species[i]->getAbundance();
    }
    return ab;
}

int Arena::getTotalAbundance(){
    int i;
    int ab=0;
    for(i=1;i<=maxsp;i++){
        ab += species[i]->getAbundance();
    }
    return ab;
}

bool Arena::findPresent(int species_id, Position p){
    return species[species_id]->isPresent(p);
}

std::list<Individual*> Arena::getPresent(int species_id,Position p){
    return species[species_id]->getPresent(p);
}

/* This will add a List of neighbours to an individual */
void Arena::addAffectedByMe(Individual *ind){
    int j;
    int sp = ind->getSpeciesId();
    Position p = ind->getPosition();
    double radius = ind->getRadius(); /* will look for inds within this radius of me */

    for(j=1;j<=maxsp;j++){
        if(species[j]->affectedBy(sp)){
            ind->addAffectedByMeNeighbourList(species[j]->getPresent(p,radius));
        }
    }
}

double Arena::getTotalTime(){
    return totalTime;
}
double Arena::getWidth(){
    return width;
}
double Arena::getHeight(){
    return height;
}

int Arena::getSpNum(){
    return maxsp;
}

Position Arena::boundaryCondition(Position p){
    switch(bcond){
        case(1):
            /* REFLEXIVE */
            while(p.x <0 || p.x > width){
                if(p.x < 0) p.x = -p.x;
                if(p.x > width) p.x = width - (p.x - width);
            }
            while(p.y <0 || p.y > height){
                if(p.y < 0) p.y = -p.y;
                if(p.y > height) p.y = height - (p.y - height);
            }
            break;

        case(2):
            /* PERIODIC */
            while(p.x < 0) p.x += width;
            while(p.x > width) p.x -= width;;
            while(p.y < 0) p.y += height;
            while(p.y > height) p.y -= height;;
            break;

        case(0):
            /* ABSORTIVE */
            if(p.x < 0 || p.x > width || p.y < 0 || p.y > height) { p.x = -1; p.y=-1; }
            break;
        default:
            Rcpp::warning("Unsuported boundary condition");
    }
    return p;
}

void Arena::addToHistory(int sp, unsigned long id, double x, double y, double beginT, double endT){
    history->sp_list.push_back(sp);
    history->id_list.push_back(id);
    history->x_list.push_back(x);
    history->y_list.push_back(y);
    history->beginTime_list.push_back(beginT);
    history->endTime_list.push_back(endT);
}

double Arena::getStressValue(Position p){
    return p.x/width;
}
