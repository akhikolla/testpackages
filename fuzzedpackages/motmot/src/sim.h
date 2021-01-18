#ifndef __SIM_H__
#define __SIM_H__

#include <vector>
#include "tree.h"   

using std::vector;

class Sim
{
    Tree tree;

    int num_segment;
    int total_time_steps;               // Total time / dt for clade.
    int num_species;                    // How many species in a given time-segment of the tree?
    int num_traits;
    double dt;  
    double a;
    double rate;                        // BM sigma^2.
    double limit;                       // trait value limits
    double time;                        // keep track of progression thru tree
    double sumDist;                     // Distance between two species in traitspace.
    double sumSqDist;
    vector<int> segment_steps;          // Number of time steps in each segment.
    vector<double> dists;
    vector<double> sqDists;

    vector<vector<double> > s;           // matrix of times lineages come into sympatry
    vector<vector<double> > al;          // matrix of times lineages leave sympatry

    /*  Modify trait values for segment between speciation
        events. Needs number of time steps within that segment  */
    void evolve_segment(int&);

    /*  Update trait values for one time step.  */
    void step_segment();

    /*  Do one evolutionary step on one species  
        This is where all the runtime is spent!  */
    void step_species(int&);
    void interaction(int, int);
    void update_distance(int, int);

    /*  An approximation to R's pnorm function  */
    double pnorm(double);

public:
    
    // Trait values.
    vector<vector<double> > tval;    

    /*  Copy R's inputs into class data.  */
    void set_values(double&, double&, double &a, double[], Tree&, int &nt, double[], double[], double&); 

    /*  Perform simulation on whole tree, modifying and lengthening tvals.  */
    void path();        
};

#endif
