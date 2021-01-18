#include <math.h>
#include <vector>
#include <random>

#include "tree.h"
#include "sim.h"

extern "C" void pathsim(int &ntip, double &dt, double &rate, double &a, double r_intervals[], 
    int splitters[], double trait_val[], int &nt, double symp[], double allo[], double &lim)
{
    /* -------------- INITIALISE TREE ----------------- */
    Tree tree;
    tree.num_tips = ntip;
    tree.total_time = 0;
    int n_interval = ntip-1;
    tree.speciators.assign(n_interval, 0);
    for(int i=0; i<n_interval; ++i)
    {
        tree.total_time += r_intervals[i];
        tree.speciators[i] = splitters[i];  
    }
    /* ------------------------------------------------ */


    /* --------------- RUN SIMULATION ----------------- */
    Sim sim;
    sim.set_values(dt, rate, a, r_intervals, tree, nt, symp, allo, lim);        
    sim.path();
    /* ------------------------------------------------ */
    

    /* -------------- RETURN VALS TO R ---------------- */
    for(int i = 0; i < ntip; ++i)
    {
        for (int j = 0; j < nt; ++j)
        {
            trait_val[j + i*nt] = sim.tval[i][j];
        } 
    }
    /* ------------------------------------------------ */
}
