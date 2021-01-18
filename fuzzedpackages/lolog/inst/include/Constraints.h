#ifndef CONSTRAINTS_H_
#define CONSTRAINTS_H_

#include "Constraint.h"
#include "ParamParser.h"

namespace lolog{

/*!
 * Restricts the sample space to the set of networks with all nodes having degree
 * >=lower and <=upper
 *
 */
template<class Engine>
class BoundedDegree : public BaseConstraint< Engine >{
protected:
	int upper;
	int lower;
	std::vector<int> scratch;
	double dist;
	double lastDist;
public:

	BoundedDegree() : upper(10000000), lower(0), dist(0.0), lastDist(0.0){}

	BoundedDegree(int low, int up) : dist(0.0), lastDist(0.0){
		lower=low;
		upper=up;
	}

	BoundedDegree(List params) : dist(0.0), lastDist(0.0){
		ParamParser p(name(), params);
		lower = p.parseNext<int>("lower");
		upper = p.parseNext<int>("upper");
		p.end();
	}

	std::string name(){
		return "boundedDegree";
	}

	/*!
	 * calculate how many steps away the constraint is from being satisfied
	 */
	double initialize(const BinaryNet<Engine>& net){
		dist = 0.0;
		for(int i=0;i<net.size();i++){
			int deg = net.degree(i);
			if(deg>upper)
				dist+= deg-upper;
			if(deg<lower)
				dist+= lower - deg;
		}
		return dist;
	}

	//dyad update
	double dyadUpdateDistance(const BinaryNet<Engine>& net, const int& from, const int&to){
		lastDist = dist;
		bool addingEdge = !net.hasEdge(from,to);
		int dfrom = net.degree(from);
		int dto = net.degree(to);
		//if(dfrom<lower || dto<lower || dfrom>upper || dto>upper)
		//	::Rf_error("Network degrees outside degree bounds");
		if(addingEdge){
			if(dfrom<lower)
				dist--;
			else if(dfrom>=upper)
				dist++;
			if(dto<lower)
				dist--;
			else if(dto>=upper)
				dist++;
		}else{
			if(dfrom<=lower)
				dist++;
			else if(dfrom>upper)
				dist--;
			if(dto<=lower)
				dist++;
			else if(dto>upper)
				dist--;
		}
		return dist;
	}


	void rollback(const BinaryNet<Engine>& net){
		dist = lastDist;

	}

	//vertex update
	double discreteVertexUpdateDistance(const BinaryNet<Engine>& net,
			const int& vert, const int& variable, const int& newValue){
		return dist;
	}

	double continVertexUpdateDistance(const BinaryNet<Engine>& net, const int& vert,
			const int& variable, const double& newValue){
		return dist;
	}

	bool isOrderIndependent(){
		return true;
	}

};

typedef Constraint<Directed, BoundedDegree<Directed> > DirectedBoundedDegreeConstraint;
typedef Constraint<Undirected, BoundedDegree<Undirected> > UndirectedBoundedDegreeConstraint;


}


#endif /* CONSTRAINTS_H_ */
