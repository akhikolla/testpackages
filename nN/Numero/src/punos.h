/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#ifndef punos_INCLUDED
#define punos_INCLUDED

#include <string>
#include <vector>
#include <map>
#include "medusa.h"
#include "scriptum.h"

namespace punos {

  /*
   * Coordinate information for visualization.
   */
  struct District {
    medusa::mdreal x;
    medusa::mdreal y;
    std::pair<medusa::mdreal, medusa::mdreal> radii;
    std::pair<medusa::mdreal, medusa::mdreal> angles;
  };

  /*
   *
   */
  class Topology {
  private:
    void* buffer;
  public:
    Topology();

    /* Create a set of unconnected centroids. */
    Topology(const medusa::mdsize);
    
    /* Create a 3D network structure. The first input sets the vertical
       positions of the planes. The second input sets the map radius
       (as the number of concentric circles of districts). Use the size()
       function to verify the total number of districts and radius() to
       check the exact space occupied by the circular plane.  */
    Topology(const std::vector<medusa::mdreal>&, const medusa::mdsize);

    /* Create a 3D network structure from existing districts. */
    Topology(const std::vector<medusa::mdreal>&,
	     const std::vector<District>&);

    /* Copy contents from the input. */
    Topology(const Topology&);
    void operator=(const Topology&);
    
    /* Free resources. */
    ~Topology();
    
    /* Return district characteristics for visualization. */
    District operator[](const medusa::mdsize) const;

    /* Number of component levels. */
    medusa::mdsize depth() const;

    /* Estimate smoothed sum of values on the network. The first input
       contains the best-matching districts, and the second contains the
       sample values. */
    std::vector<medusa::mdreal>
    diffuse(const std::vector<medusa::mdsize>&,
	    const std::vector<medusa::mdreal>&) const;

    /* Estimate smoothed sum of values on the network. The first input
       specifies the sample depths, see stratify(). The second input
       contains the best-matching districts, and the third contains the
       sample values. */
    std::vector<std::vector<medusa::mdreal> >
    diffuse(const std::vector<medusa::Site>&,
	    const std::vector<medusa::mdsize>&,
	    const std::vector<medusa::mdreal>&) const;

    /* Calculate spatial distance between two map districts. */
    medusa::mdreal distance(const medusa::mdsize,
			    const medusa::mdsize) const;    

    /* Write single-character highlight labels on the map districts.
       IMPORTANT: The first two inputs define the center point of the map
       in the graphics device units, not the native spatial units. */
    scriptum::Frame highlight(const medusa::mdreal, const medusa::mdreal,
			      const std::vector<char>&,
			      const std::vector<scriptum::Color>&,
			      const scriptum::Style&) const;

    /* Import a network structure from a file. Any previous
       contents is discarded. Returns an error message if failed. */
    std::string import(const std::string&);

    /* Interpolate district prototypes according to pivot points. The
       input contains seed profiles that will be spread over the map. */
    std::vector<std::vector<medusa::mdreal> >
    interpolate(const std::vector<std::vector<medusa::mdreal> >&) const;
    
    /* Vertical positions of level planes. */
    std::vector<medusa::mdreal> levels() const;

    /* List of connected districts, including self. */
    std::vector<medusa::mdsize> neighbors(const medusa::mdsize) const;

    /* Paint map districts.
       IMPORTANT: The first two inputs define the center point of the map
       in the graphics device units, not the native spatial units. */
    scriptum::Frame paint(const medusa::mdreal, const medusa::mdreal,
			  const std::vector<scriptum::Color>&,
			  const scriptum::Style&) const;
    
    /* Maximum extent from origin. */
    medusa::mdreal radius() const;

    /* Set up neighborhood network. The input must be positive and at most
       half of map radius. Returns true if any edges were created. */
    bool rewire(const medusa::mdreal);

    /* Save network structure in a text file. */
    unsigned long save(const std::string&) const;

    /* Current neighborhood radius. */
    medusa::mdreal sigma() const;

    /* Number of districts. */
    medusa::mdsize size() const;

    /* Determine adjacent planes (i.e. the layer in between)
       according to vertical position. */
    medusa::Site stratify(const medusa::mdreal) const;

    /* Link weight between two map districts. */
    medusa::mdreal weight(const medusa::mdsize,
			  const medusa::mdsize) const;

    /* Write text labels on map districts.
       IMPORTANT: The first two inputs define the center point of the map
       in graphics device units, not the native spatial units. */
    scriptum::Frame write(const medusa::mdreal, const medusa::mdreal,
			  const std::vector<std::string>&,
			  const std::vector<scriptum::Color>&,
			  const scriptum::Style&) const;
  };

  /* Version information. */
  extern std::string version();
}

#endif /* punos_INCLUDED */
