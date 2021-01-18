#include "VariableFactory.h"

Element* VariableFactory::GetVariable(const std::string& name, integer n, integer m, integer p)
{
	if (name == "Euclidean") {
		return new EucVariable(n, m);
	} else if (name == "Sphere") {
		return new SphereVariable(n);
	} else if (name == "Stiefel") {
		return new StieVariable(n, p);
	} else if (name == "Oblique") {
		return new ObliqueVariable(n, m);
	} else if (name == "LowRank") {
		return new LowRankVariable(n, m, p);
	} else if (name == "OrthGroup") {
		return new OrthGroupVariable(n);
	} else if (name == "L2Sphere") {
		return new L2SphereVariable(n);
	} else if (name == "SPDManifold") {
		return new SPDVariable(n);
	} else if (name == "CpxNStQOrth") {
		throw ManifoldOptimException("CSOVariable type is not currently implemented");
		// return new CSOVariable(n, p);
	} else if (name == "Grassmann") {
		return new GrassVariable(n, p);
	} else {
		throw ManifoldOptimException("Variable type is not implemented in this library");
	}
}

