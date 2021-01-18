#include "ManifoldFactory.h"

Manifold* ManifoldFactory::GetManifold(const std::string& name, integer n, integer m, integer p)
{
	if (name == "Euclidean") {
		return new Euclidean(n, m);
	} else if (name == "Sphere") {
		return new Sphere(n);
	} else if (name == "Stiefel") {
		return new Stiefel(n, p);
	} else if (name == "Oblique") {
		throw ManifoldOptimException("Oblique Manifold type is not currently implemented");
	} else if (name == "LowRank") {
		return new LowRank(n, m, p);
	} else if (name == "OrthGroup") {
		return new OrthGroup(n);
	} else if (name == "L2Sphere") {
		throw ManifoldOptimException("L2Sphere Manifold type is not currently implemented");
	} else if (name == "SPDManifold") {
		return new SPDManifold(n);
	} else if (name == "CpxNStQOrth") {
		throw ManifoldOptimException("CpxNStQOrth Manifold type is not currently implemented");
		// return new CpxNStQOrth(n, p);
	} else if (name == "Grassmann") {
		return new Grassmann(n, p);
	} else {
		throw ManifoldOptimException("Manifold type is not implemented in this library");
	}
}


