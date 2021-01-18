/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#ifndef koho_INCLUDED
#define koho_INCLUDED

#include <string>
#include <vector>
#include "medusa.h"
#include "punos.h"

namespace koho {

  /*
   *
   */
  struct Resident {

    /* Name of the data point. */
    std::string identity;

    /* Map district. */
    medusa::mdsize district;

    /* Distance in data space to district prototype. */
    medusa::mdreal residual;
  };

  /*
   *
   */
  class Model {
  private:
    void* buffer;
  public:
    Model();

    /* Set up a model with the given topology. The second input sets
       the maximum number of training samples per cycle. The third
       input sets the balancing coefficient for spatial point histogram:
       0.0 means no balancing and 1.0 means maximum balancing. */
    Model(const punos::Topology&, const medusa::mdsize,
	  const medusa::mdreal);

    /* Copy the contents from the input. */
    Model(const Model&);
    void operator=(const Model&);

    /* Free local resources. */
    ~Model();
    
    /* Set a prototype vector. The first input sets the district, the second
       contains the data. Any updates are incremental: valid values are
       inserted, whereas any other previous elements are left untouched. */
    std::string configure(const medusa::mdsize,
			  const std::vector<medusa::mdreal>&);
 
    /* Estimate distances to the district profiles in data space. */
    std::vector<medusa::mdreal> distance(const std::string&) const;

    /* Return data point identities. */
    std::vector<std::string> identities() const;

    /* Insert a new data point. The first input is the data point identity
       and the second input sets the values for each dimension. No missing
       values or non-uniform number of elements is allowed. Returns a non-
       empty message if failed. */
    std::string insert(const std::string&,
		       const std::vector<medusa::mdreal>&);

     /* Return the number of data dimensions. */
    medusa::mdsize order() const;
    
   /* Return the prototype vector for the specified map district. */
    std::vector<medusa::mdreal> prototype(const medusa::mdsize) const;

    /* Return the number of training data points. */
    medusa::mdsize size() const;

    /* Return current map topology. */
    punos::Topology topology() const;

    /* Train the map according to the inserted data points. The first input
       is filled with the final layout. The second input is filled with
       training errors from each cycle. The last input sets a time quota
       for training. If training error has stabilized, no further cycles are
       run. Returns a message if failed. */
    std::string train(std::vector<Resident>&, std::vector<medusa::mdreal>&,
		      const medusa::mdreal);
  };

  /*
   *
   */
  class Engine {
  private:
    void* buffer;
  public:
    Engine();

    /* Set up an engine with the given topology. */
    Engine(const punos::Topology&);

    /* Copy the contents from the input. */
    Engine(const Engine&);
    void operator=(const Engine&);

    /* Free local resources. */
    ~Engine();
    
    /* Estimate district averages using the current layout and contents. */
    std::vector<std::vector<medusa::mdreal> > average() const;

    /* Data point histograms of the current contents. */
    std::vector<std::vector<medusa::mdreal> > histograms() const;
    
    /* Insert a new one-dimensional data point. The first input is the
       point identity, the second indicates the map district, and the third
       contains the data value. Returns an error if failed. */
    std::string insert(const std::string&, const medusa::mdsize,
		       const medusa::mdreal);
 
    /* Insert a new multi-dimensional data point. The first input is the
       point identity, the second indicates the map district, and the third
       contains the data values. Returns an error if failed. */
    std::string insert(const std::string&, const medusa::mdsize,
		       const std::vector<medusa::mdreal>&);
 
    /* Shuffle values randomly. If the second input is true, shuffle with
       replacement. Returns true if values were shuffled. */
    bool shuffle(const bool);
  };

  /* Version information. */
  extern std::string version();
}

#endif /* koho_INCLUDED */
