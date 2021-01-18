#include <time.h>
#include "stdio.h"
#include "string.h"
#include <stdlib.h>
#include <unistd.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>	// GSL_SUCCEnn ...

#include <gsl/gsl_types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_min.h>

#include <gsl/gsl_roots.h>

#include "math.h"

//Custom headers and functions:
#include "Nrutil.h"
#include "covid19_model.h"
#include "time_window_utils.h"



int covid19_model (
            COVID19ParamStruct *params,
            int *out_pops_seed,
            int *out_pops_pop,
            int *out_pops_time,
            int *out_pops_S_pop,
            int *out_pops_E_pop,
            int *out_pops_I_asym_pop,
            int *out_pops_I_presym_pop,
            int *out_pops_I_sym_pop,
            int *out_pops_I_home_pop,
            int *out_pops_I_hosp_pop,
            int *out_pops_I_icu1_pop,
            int *out_pops_I_icu2_pop,
            int *out_pops_R_pop,
            int *out_pops_D_pop,
            int *out_events_pos,
            int *out_events_sym,
            int *out_events_total_hosp,
            int *out_events_total_icu,
            int *out_events_n_death)
{

    // Allocate random numbers structures:
    gsl_rng_env_setup();
    const gsl_rng_type *T1 = gsl_rng_default;
    gsl_rng *rand1 = gsl_rng_alloc(T1);
    //*****************************************

    double r0 = 0, r0_slope = 0, r0_intercept = 0;
    double m_slope = 0, m_intercept = 0;
    double dist_param = 0, dist_param_low = 0, dist_param_slope = 0;
    double dist_param_intercept = 0;
    double imm_frac_slope = 0, imm_frac_intercept = 0;
    bool dist_param_changed = false;
    bool daily_mode_on = false;

    // Import the time window data into a linked list
    TimeWindow *head_node, *current_node = NULL;
    head_node = importTimeWindowData(params->total_windows,
                                     params->input_r0,
                                     params->input_dist_param,
                                     params->input_m,
                                     params->input_imm_frac,
                                     params->input_hosp_rate,
                                     params->input_icu_rate,
                                     params->input_death_rate,
                                     params->input_recov_hosp,
                                     params->input_window_length);
    current_node = head_node;
    // Set initial parameters for funcs.h
    params->hosp_rate = current_node->hosp_rate;
    params->icu_rate = current_node->icu_rate;
    params->death_rate = current_node->death_rate;
    params->recov_hosp = current_node->recov_hosp;

    int n_pop = params->n_pop;
    int n_times;
    n_times = params->t_max / params->tau;  // t_max = 250, tau = 1

    // Events and Updates:
    int n_events = 18; //15 non-migration
    int *n_occur = new int[n_events];
    int *update_vec = new int[params->n_equations];
    int *update_vec_migrants = new int[params->n_equations];
    int *update_vec_move = new int[n_pop];
    //*****************************************

    /* *************************************************************** */
    /* Read in the Pop-specific data and Distance Matrix   :           */
    /* *************************************************************** */

    float *pop_N, *prob_ColSum; // prevalence data
    float *census_area; // prevalence data
    float **dist_mat, **prob_move;

    //TODO: replace these matrix/vector allocations from Nrutil.h with
    //pre-allocated pointer arguments, so we don't have to worry about
    //memory management in here.
    dist_mat = matrix(1, n_pop, 1, n_pop);
    prob_move = matrix(1, n_pop, 1, n_pop);
    prob_ColSum = vector(1, n_pop);
    pop_N = vector(1, n_pop);
    census_area = vector(1, n_pop);

    // compute prob_ColSum and prob_move based on dist_mat and
    // dist_param_low.
    for(int i = 1; i <= n_pop; i++){
        pop_N[i] = params->input_N_pops[i-1];//TODO remove pop_N/census_area, just use input_*
        census_area[i] = params->input_census_area[i-1];
        prob_ColSum[i] = 0.0;

        for(int j = 1; j <= n_pop; j++){
            dist_mat[i][j] = params->input_dist_vec[(i-1)+(j-1)*n_pop];
            prob_move[i][j] = 0;

            if(i != j){
                prob_move[i][j] = 1 / exp(dist_mat[i][j] / dist_param_low);
            }

            prob_ColSum[j] = prob_ColSum[j] + prob_move[i][j];
            // printf("prob_move[%d][%d]: %f, ", i, j, prob_move[i][j]);
        }
         // printf("\n");
    }

    params->pop_N = pop_N;
    params->census_area = census_area;
    params->dist_mat = dist_mat;
    params->prob_move = prob_move;

    /* ***************************************************  */
    /* Allocate and store populaitons   :                   */
    /* ***************************************************  */
    // TRACKED POPULATIONS
    // # ROWS = POPS
    // # Cols = 1: Current, 2: Future
    int **S_pop, **E_pop, **I_asym_pop;
    int **I_presym_pop, **I_sym_pop, **I_home_pop;
    int **I_hosp_pop, **I_icu1_pop, **I_icu2_pop;
    int **R_pop, **D_pop;
    S_pop = imatrix(1, n_pop, 1, 2);
    E_pop = imatrix(1, n_pop, 1, 2);
    I_asym_pop = imatrix(1, n_pop, 1, 2);
    I_presym_pop = imatrix(1, n_pop, 1, 2);
    I_sym_pop = imatrix(1, n_pop, 1, 2);
    I_home_pop = imatrix(1, n_pop, 1, 2);
    I_hosp_pop = imatrix(1, n_pop, 1, 2);
    I_icu1_pop = imatrix(1, n_pop, 1, 2);
    I_icu2_pop = imatrix(1, n_pop, 1, 2);
    R_pop = imatrix(1, n_pop, 1, 2);
    D_pop = imatrix(1, n_pop, 1, 2);


    // MIGRANTS (from = cols, to = rows)

    int **S_move, **I_move;

    S_move = imatrix(1, n_pop, 1, n_pop);
    I_move = imatrix(1, n_pop, 1, n_pop);

    for(int i=1;i<=n_pop;i++){
        for(int j=1;j<=2;j++){
            S_pop[i][j] = 0.0; E_pop[i][j] = 0.0; I_asym_pop[i][j] = 0.0;
            I_presym_pop[i][j] = 0.0; I_sym_pop[i][j] = 0.0;
            I_home_pop[i][j] = 0.0; I_hosp_pop[i][j] = 0.0;
            I_icu1_pop[i][j] = 0.0; I_icu2_pop[i][j] = 0.0;
            R_pop[i][j] = 0.0; D_pop[i][j] = 0.0;
        }

        for(int j=1;j<=n_pop;j++){
            S_move[i][j] = 0;
            I_move[i][j] = 0;
        }
    }


    //*****************************************

    /* ***************************************************  */
    /* INITIAL CONDITIONS  :                                */
    /* ***************************************************  */

    //********************

    COVID19PopStruct AllPops;

    AllPops.S_pop = S_pop;
    AllPops.E_pop = E_pop;
    AllPops.I_asym_pop = I_asym_pop;
    AllPops.I_presym_pop = I_presym_pop;
    AllPops.I_sym_pop = I_sym_pop;
    AllPops.I_home_pop = I_home_pop;
    AllPops.I_hosp_pop = I_hosp_pop;
    AllPops.I_icu1_pop = I_icu1_pop;
    AllPops.I_icu2_pop = I_icu2_pop;
    AllPops.R_pop = R_pop;
    AllPops.D_pop = D_pop;

    int n_I_move_total = 0;

    COVID19MoveMatStruct MovePops;

    MovePops.S_move = S_move;
    MovePops.I_move = I_move;

    //*****************************************
    //*****************************************
    //*****************************************

    // RUN THE SIMULATION:

    int out_pops_line=0;
    int out_events_line=0;

    for(int seed_i = 0; seed_i < params->n_realz; seed_i++){ // n_realz how many times running the model
        int this_seed = params->input_realz_seeds[seed_i];
        //printf("seed_i=%d this_seed=%d\n", seed_i, this_seed);
        gsl_rng_set(rand1, this_seed);
        /* ***************************************************  */
        /* INITIAL CONDITIONS  :                                */
        /* ***************************************************  */

        for(int k=1;k<=n_pop;k++){
            for(int j=1;j<=2;j++){

                // SET AS INPUT_*_POPS
                S_pop[k][j] = params->input_S_pops[(k-1)];
                E_pop[k][j] = params->input_E_pops[(k-1)];
                I_asym_pop[k][j] = params->input_I_asym_pops[(k-1)];
                I_presym_pop[k][j] = params->input_I_presym_pops[(k-1)];
                I_sym_pop[k][j] = params->input_I_sym_pops[(k-1)];
                I_home_pop[k][j] = params->input_I_home_pops[(k-1)];
                I_hosp_pop[k][j] = params->input_I_hosp_pops[(k-1)];
                I_icu1_pop[k][j] = params->input_I_icu1_pops[(k-1)];
                I_icu2_pop[k][j] = params->input_I_icu2_pops[(k-1)];
                R_pop[k][j] = params->input_R_pops[(k-1)];
                D_pop[k][j] = params->input_D_pops[(k-1)];
            }

            for(int j=1;j<=n_pop;j++){
                S_move[k][j] = 0;
                I_move[k][j] = 0;
            }
        }

        AllPops.S_pop = S_pop;
        AllPops.E_pop = E_pop;
        AllPops.I_asym_pop = I_asym_pop;
        AllPops.I_presym_pop = I_presym_pop;
        AllPops.I_sym_pop = I_sym_pop;
        AllPops.I_home_pop = I_home_pop;
        AllPops.I_hosp_pop = I_hosp_pop;
        AllPops.I_icu1_pop = I_icu1_pop;
        AllPops.I_icu2_pop = I_icu2_pop;
        AllPops.R_pop = R_pop;
        AllPops.D_pop = D_pop;

        // TEST:
        // exit(1);

        /* ***************************************************  */
        /* BEGIN  :                                             */
        /* ***************************************************  */

        for(int t = 1; t <= n_times; t++){ //n_times

            // Check for moving to the next time window (1st window has length 0)
            if (current_node->window_length < 1)
            {
                if (current_node->next != NULL)
                {
                    // Get the next time window
                    current_node = current_node->next;
                    // Check if running in daily mode or time window mode
                    if (current_node->window_length == 1)
                    {
                        daily_mode_on = true;
                    }
                    else
                    {
                        daily_mode_on = false;
                    }
                }
                else
                {
                    /*
                     * No more time windows but loop is still going.
                     * Use initial values for remainder of the loop.
                     */
                    current_node = head_node;
                    current_node->window_length = n_times - t;
                }

                // Min/max for R0
                r0_slope = current_node->getR0Slope();
                r0_intercept = current_node->getR0Intercept(t - 1);

                if (dist_param != current_node->dist_param)
                {
                    dist_param_changed = true;

                    // dist_param
                    dist_param_low = current_node->getMinDistParam();
                    dist_param_slope = current_node->getDistParamSlope();
                    dist_param_intercept = current_node->getDistParamIntercept(t - 1);
                }
                else
                {
                    dist_param_changed = false;
                }

                // m
                m_slope = current_node->getMSlope();
                m_intercept = current_node->getMIntercept(t - 1);

                // imm_frac
                imm_frac_slope = current_node->getImmFracSlope();
                imm_frac_intercept = current_node->getImmFracIntercept(t - 1);
            }

            //*****************************************
            // Effects of time window data

            if (daily_mode_on)
            {
                params->m = current_node->m;
                params->imm_frac = current_node->imm_frac;
                r0 = current_node->r0;
                params->beta = calculateBeta(r0, params);
                dist_param = current_node->dist_param;
                params->hosp_rate = current_node->hosp_rate;
                params->icu_rate = current_node->icu_rate;
                params->death_rate = current_node->death_rate;
                params->recov_hosp = current_node->recov_hosp;
            }
            else
            {
                params->m = m_slope * t + m_intercept;
                params->imm_frac = imm_frac_slope * t + imm_frac_intercept;
                r0 = r0_slope * t + r0_intercept;
                params->beta = calculateBeta(r0, params);
                dist_param = dist_param_slope * t + dist_param_intercept;
            }

            // Only deal with prob_move if dist_param changes
            if (dist_param_changed)
            {
                for(int k = 1; k <= n_pop; k++){
                    // prob_ColSum[k] = 0.0;

                    for(int j = 1; j <= n_pop; j++){
                        prob_move[k][j] = 0.0;

                        if(k != j){
                            prob_move[k][j] = 1 / exp(dist_mat[k][j] / dist_param);
                        }
                        // prob_ColSum[j] = prob_ColSum[j] + prob_move[k][j];
                    }
                }

                params->prob_move = prob_move;
            }


            //*****************************************
            // 1 tau leap (internal)
            // printf("up to tau_leap_1step.\n");

            for(int this_pop = 1; this_pop <= n_pop; this_pop++){
                // printf("up to tau_leap_1step.\n");
                tau_leap_1step( n_occur,
                                this_pop,
                                params,
                                AllPops,
                                rand1,
                                n_events);
                update_pops(update_vec,
                            n_occur,
                            this_pop,
                            t,
                            this_seed,
                            params,
                            AllPops,
                            rand1,
                            out_events_pos       +out_events_line,
                            out_events_sym       +out_events_line,
                            out_events_total_hosp+out_events_line,
                            out_events_total_icu +out_events_line,
                            out_events_n_death   +out_events_line);

                out_events_line++;
                AllPops.S_pop[this_pop][2] = update_vec[0];
                AllPops.E_pop[this_pop][2] = update_vec[1];
                AllPops.I_asym_pop[this_pop][2] = update_vec[2];
                AllPops.I_presym_pop[this_pop][2] = update_vec[3];
                AllPops.I_sym_pop[this_pop][2] = update_vec[4];
                AllPops.I_home_pop[this_pop][2] = update_vec[5];
                AllPops.I_hosp_pop[this_pop][2] = update_vec[6];
                AllPops.I_icu1_pop[this_pop][2] = update_vec[7];
                AllPops.I_icu2_pop[this_pop][2] = update_vec[8];
                AllPops.R_pop[this_pop][2] = update_vec[9];
                AllPops.D_pop[this_pop][2] = update_vec[10];

                // printf("S_pop[%d]: %d\n", this_pop, AllPops.S_pop[this_pop][2]);
                // printf("E_pop[%d]: %d\n", this_pop, AllPops.E_pop[this_pop][2]);
                // printf("I_sym_pop[%d]: %d\n", this_pop, AllPops.I_sym_pop[this_pop][2]);
                // printf("I_asym_pop[%d]: %d\n", this_pop, AllPops.I_asym_pop[this_pop][2]);
                // printf("H_pop[%d]: %d\n", this_pop, AllPops.H_pop[this_pop][2]);
                // printf("R_pop[%d]: %d\n", this_pop, AllPops.R_pop[this_pop][2]);


                //// MIGRATION:

                // Calculate where to move folks:
                // printf("up to move_pop()\n");
                // Move S
                if(n_occur[15] > 0){
                    move_pops(update_vec_move, n_occur[15], this_pop, params, rand1);

                    for(int k=0;k<n_pop;k++){
                        // printf("update_move[%d]: %d\n", (k+1), update_move[k]);
                        MovePops.S_move[k+1][this_pop] = update_vec_move[k];
                    }
                }

                // Move I
                n_I_move_total = (n_occur[16] + n_occur[17]);
                if(n_I_move_total > 0){
                    move_pops(update_vec_move, n_I_move_total, this_pop, params, rand1);

                    for(int k=0;k<n_pop;k++){
                        // printf("update_I_move[%d][%d]: %d\n", this_pop, (k+1), update_move[k]);
                        MovePops.I_move[k+1][this_pop] = update_vec_move[k];
                    }
                }

                // for(k = 0; k < n_pop; k++){
                //   printf("S_move[%d][%d]: %d\n", this_pop, k+1, MovePops.S_move[k+1][this_pop]);
                //   printf("I_move[%d][%d]: %d\n", this_pop, k+1, MovePops.I_move[k+1][this_pop]);
                // }
            }//end j loop (pop)

            //*****************************************
            //*****************************************
            // Migrants can cause exposure (I_asym) or migrants can become exposed (S)

            // printf("up to migrants.\n");
            for(int this_pop = 1; this_pop <= n_pop; this_pop++){
                update_pop_migrants(update_vec_migrants,
                                    this_pop,
                                    params,
                                    AllPops,
                                    MovePops,
                                    rand1);

                AllPops.S_pop[this_pop][2] = update_vec_migrants[0];
                AllPops.E_pop[this_pop][2] = update_vec_migrants[1];
                AllPops.I_asym_pop[this_pop][2] = update_vec_migrants[2];
                AllPops.I_presym_pop[this_pop][2] = update_vec_migrants[3];
                AllPops.I_sym_pop[this_pop][2] = update_vec_migrants[4];
                AllPops.I_home_pop[this_pop][2] = update_vec_migrants[5];
                AllPops.I_hosp_pop[this_pop][2] = update_vec_migrants[6];
                AllPops.I_icu1_pop[this_pop][2] = update_vec_migrants[7];
                AllPops.I_icu2_pop[this_pop][2] = update_vec_migrants[8];
                AllPops.R_pop[this_pop][2] = update_vec_migrants[9];
                AllPops.D_pop[this_pop][2] = update_vec_migrants[10];

                // printf("S_pop[%d]: %d\n", this_pop, AllPops.S_pop[this_pop][2]);
                // printf("E_pop[%d]: %d\n", this_pop, AllPops.E_pop[this_pop][2]);
                // printf("I_sym_pop[%d]: %d\n", this_pop, AllPops.I_sym_pop[this_pop][2]);
                // printf("I_asym_pop[%d]: %d\n", this_pop, AllPops.I_asym_pop[this_pop][2]);
                // printf("H_pop[%d]: %d\n", this_pop, AllPops.H_pop[this_pop][2]);
                // printf("R_pop[%d]: %d\n", this_pop, AllPops.R_pop[this_pop][2]);
            }// end j loop (this_pop)


            //*****************************************
            //*****************************************
            // REDUCING LOOPING BY COMBINING STEPS:

            for(int this_pop=1;this_pop<=n_pop;this_pop++){
                // Check for negatives:
                if(AllPops.S_pop[this_pop][2] < 0) AllPops.S_pop[this_pop][2] = 0;
                if(AllPops.E_pop[this_pop][2] < 0) AllPops.E_pop[this_pop][2] = 0;
                if(AllPops.I_asym_pop[this_pop][2] < 0) AllPops.I_asym_pop[this_pop][2] = 0;
                if(AllPops.I_presym_pop[this_pop][2] < 0) AllPops.I_presym_pop[this_pop][2] = 0;
                if(AllPops.I_sym_pop[this_pop][2] < 0) AllPops.I_sym_pop[this_pop][2] = 0;
                if(AllPops.I_home_pop[this_pop][2] < 0) AllPops.I_home_pop[this_pop][2] = 0;
                if(AllPops.I_hosp_pop[this_pop][2] < 0) AllPops.I_hosp_pop[this_pop][2] = 0;
                if(AllPops.I_icu1_pop[this_pop][2] < 0) AllPops.I_icu1_pop[this_pop][2] = 0;
                if(AllPops.I_icu2_pop[this_pop][2] < 0) AllPops.I_icu2_pop[this_pop][2] = 0;
                if(AllPops.R_pop[this_pop][2] < 0) AllPops.R_pop[this_pop][2] = 0;
                if(AllPops.D_pop[this_pop][2] < 0) AllPops.D_pop[this_pop][2] = 0;

                //Recalculate pop size after all movement:
                pop_N[this_pop] =
                AllPops.S_pop[this_pop][2] + AllPops.E_pop[this_pop][2] + AllPops.I_asym_pop[this_pop][2] +
                AllPops.I_presym_pop[this_pop][2] + AllPops.I_sym_pop[this_pop][2] + AllPops.I_home_pop[this_pop][2] +
                AllPops.I_hosp_pop[this_pop][2] + AllPops.I_icu1_pop[this_pop][2] + AllPops.I_icu2_pop[this_pop][2] +
                AllPops.R_pop[this_pop][2] - AllPops.D_pop[this_pop][2];

                // write one output line.
                int one_or_two;
                if(t==1){
                    one_or_two=1;
                }
                else{
                    one_or_two=2;
                }

            	out_pops_seed[out_pops_line] = this_seed;
            	out_pops_pop[out_pops_line] = this_pop;
            	out_pops_time[out_pops_line] = t;
            	out_pops_S_pop[out_pops_line] = AllPops.S_pop[this_pop][one_or_two];
            	out_pops_E_pop[out_pops_line] = AllPops.E_pop[this_pop][one_or_two];
            	out_pops_I_asym_pop[out_pops_line] = AllPops.I_asym_pop[this_pop][one_or_two];
            	out_pops_I_presym_pop[out_pops_line] = AllPops.I_presym_pop[this_pop][one_or_two];
            	out_pops_I_sym_pop[out_pops_line] = AllPops.I_sym_pop[this_pop][one_or_two];
            	out_pops_I_home_pop[out_pops_line] = AllPops.I_home_pop[this_pop][one_or_two];
            	out_pops_I_hosp_pop[out_pops_line] = AllPops.I_hosp_pop[this_pop][one_or_two];
            	out_pops_I_icu1_pop[out_pops_line] = AllPops.I_icu1_pop[this_pop][one_or_two];
            	out_pops_I_icu2_pop[out_pops_line] = AllPops.I_icu2_pop[this_pop][one_or_two];
            	out_pops_R_pop[out_pops_line] = AllPops.R_pop[this_pop][one_or_two];
            	out_pops_D_pop[out_pops_line] = AllPops.D_pop[this_pop][one_or_two];
            	out_pops_line++;

                // Reset movement:
                for(int k=1;k<=n_pop;k++){
                    MovePops.S_move[this_pop][k] = 0;
                    MovePops.I_move[this_pop][k] = 0;
                }

                // Update Pops 1 step;
                AllPops.S_pop[this_pop][1] = AllPops.S_pop[this_pop][2];
                AllPops.E_pop[this_pop][1] = AllPops.E_pop[this_pop][2];
                AllPops.I_asym_pop[this_pop][1] = AllPops.I_asym_pop[this_pop][2];
                AllPops.I_presym_pop[this_pop][1] = AllPops.I_presym_pop[this_pop][2];
                AllPops.I_sym_pop[this_pop][1] = AllPops.I_sym_pop[this_pop][2];
                AllPops.I_home_pop[this_pop][1] = AllPops.I_home_pop[this_pop][2];
                AllPops.I_hosp_pop[this_pop][1] = AllPops.I_hosp_pop[this_pop][2];
                AllPops.I_icu1_pop[this_pop][1] = AllPops.I_icu1_pop[this_pop][2];
                AllPops.I_icu2_pop[this_pop][1] = AllPops.I_icu2_pop[this_pop][2];
                AllPops.R_pop[this_pop][1] = AllPops.R_pop[this_pop][2];
                AllPops.D_pop[this_pop][1] = AllPops.D_pop[this_pop][2];

                AllPops.S_pop[this_pop][2] = 0;
                AllPops.E_pop[this_pop][2] = 0;
                AllPops.I_asym_pop[this_pop][2] = 0;
                AllPops.I_presym_pop[this_pop][2] = 0;
                AllPops.I_sym_pop[this_pop][2] = 0;
                AllPops.I_home_pop[this_pop][2] = 0;
                AllPops.I_hosp_pop[this_pop][2] = 0;
                AllPops.I_icu1_pop[this_pop][2] = 0;
                AllPops.I_icu2_pop[this_pop][2] = 0;
                AllPops.R_pop[this_pop][2] = 0;
                AllPops.D_pop[this_pop][2] = 0;

                // printf("S_pop1[%d]: %d\n", this_pop, AllPops.S_pop[this_pop][1]);
                // printf("E_pop1[%d]: %d\n", this_pop, AllPops.E_pop[this_pop][1]);
                // printf("I_sym_pop1[%d]: %d\n", this_pop, AllPops.I_sym_pop[this_pop][1]);
                // printf("I_asym_pop1[%d]: %d\n", this_pop, AllPops.I_asym_pop[this_pop][1]);
                // printf("H_pop1[%d]: %d\n", this_pop, AllPops.H_pop[this_pop][1]);
                // printf("R_pop1[%d]: %d\n", this_pop, AllPops.R_pop[this_pop][1]);

            }
            // fclose(fp2);

            // Remove a day from the time window
            current_node->window_length--;

        }//end t loop (time)

    } // end i loop (realz)

    //*****************************************

    /* ***************************************************  */
    /* FREE MEMORY ALLOCATIONS   :                          */
    /* ***************************************************  */

    // Clean up memory used by time windows
    clearTimeWindows(head_node);

    delete[] n_occur;
    delete[] update_vec;
    delete[] update_vec_migrants;
    delete[] update_vec_move;

    free_matrix(dist_mat, 1, n_pop, 1, n_pop);
    free_matrix(prob_move, 1, n_pop, 1, n_pop);
    free_vector(prob_ColSum, 1, n_pop);
    free_vector(pop_N, 1, n_pop);
    free_vector(census_area, 1, n_pop);

    free_imatrix(S_pop, 1, n_pop, 1, 2);
    free_imatrix(E_pop, 1, n_pop, 1, 2);
    free_imatrix(I_asym_pop, 1, n_pop, 1, 2);
    free_imatrix(I_presym_pop, 1, n_pop, 1, 2);
    free_imatrix(I_sym_pop, 1, n_pop, 1, 2);
    free_imatrix(I_home_pop, 1, n_pop, 1, 2);
    free_imatrix(I_hosp_pop, 1, n_pop, 1, 2);
    free_imatrix(I_icu1_pop, 1, n_pop, 1, 2);
    free_imatrix(I_icu2_pop, 1, n_pop, 1, 2);
    free_imatrix(R_pop, 1, n_pop, 1, 2);
    free_imatrix(D_pop, 1, n_pop, 1, 2);

    free_imatrix(S_move, 1, n_pop, 1, n_pop);
    free_imatrix(I_move, 1, n_pop, 1, n_pop);
    //free_imatrix(output_pops, 1, dim_out2, 1, dim_out1);

    return 0;
}



//**************** Funcs.h ******************************

void trans_type_beta(double& beta_scaled, int this_pop, double infect_sum,
                     COVID19ParamStruct *Params, COVID19PopStruct AllPops, gsl_rng *rand1){

    double noise_temp = 0.0;
    double pop_dens = 0.0;
    double beta_dens = 0.0;

    noise_temp = gsl_ran_gaussian(rand1, Params->stoch_sd);

    switch(Params->trans_type){

    case 1:
        // FREQUENCY-DEPENDENT TRANSMISSION

        beta_scaled = fabs(Params->beta / Params->pop_N[this_pop] *
          (1 + (noise_temp) / pow(infect_sum, 0.5)));

        break;

    case 2:
        // DENSITY-DEPENDENT TRANSMISSION
        //// MONOD EQUATION (to determine rel'n btw beta and raw pop_dens)
        pop_dens = (Params->pop_N[this_pop] / Params->census_area[this_pop]);
        beta_dens = Params->beta * pop_dens / (Params->dd_trans_monod_k + pop_dens);

        beta_scaled = fabs(beta_dens / Params->pop_N[this_pop] *
          (1 + (noise_temp) /  pow(infect_sum, 0.5)));

        break;

    default:
        break;

    }// end switch

}

void tau_leap_1step(int *n_occur, int this_pop, COVID19ParamStruct *Params,
                    COVID19PopStruct AllPops, gsl_rng *rand1, int n_events){

    double infect_sum  = 0.0;
    double beta_scaled;

    double *event_prob = new double[n_events]; // size = number of events

    int i;

    // sum of infectious:
    infect_sum =
        AllPops.I_asym_pop[this_pop][1] +
        AllPops.I_presym_pop[this_pop][1] +
        AllPops.I_sym_pop[this_pop][1] +
        AllPops.I_home_pop[this_pop][1] +
        AllPops.I_hosp_pop[this_pop][1] +
        AllPops.I_icu1_pop[this_pop][1] +
        AllPops.I_icu2_pop[this_pop][1];

    if(infect_sum > 0.0){ // As long as there are some infectious hosts...

        trans_type_beta(beta_scaled,
                        this_pop,
                        infect_sum,
                        Params,
                        AllPops,
                        rand1);

    }else{
        beta_scaled = 0.0;
    }

    // printf("beta_scaled = %f\n", beta_scaled);


    //*****************************************
    // Calculate Event Probabilities:

    //// TRANSMISSION:
    // Transmission from I_asym:
    event_prob[TRANS_FROM_ASYM] =
        beta_scaled * AllPops.S_pop[this_pop][1] * AllPops.I_asym_pop[this_pop][1] * Params->frac_beta_asym;
    // Transmission from I_presym:
    event_prob[TRANS_FROM_PRESYM] =
        beta_scaled * AllPops.S_pop[this_pop][1] * AllPops.I_presym_pop[this_pop][1];
    // Transmission from I_sym:
    event_prob[TRANS_FROM_SYM] =
        beta_scaled * AllPops.S_pop[this_pop][1] * AllPops.I_sym_pop[this_pop][1];
    // Transmission from I_home:
    event_prob[TRANS_FROM_HOME] =
        beta_scaled * AllPops.S_pop[this_pop][1] * AllPops.I_home_pop[this_pop][1];
    // Transmission from I_hosp:
    event_prob[TRANS_FROM_HOSP] =
        beta_scaled * AllPops.S_pop[this_pop][1] * AllPops.I_hosp_pop[this_pop][1] * Params->frac_beta_hosp;
    // Transmission from I_icu1:
    event_prob[TRANS_FROM_ICU1] =
        beta_scaled * AllPops.S_pop[this_pop][1] * AllPops.I_icu1_pop[this_pop][1] * Params->frac_beta_hosp;
    // Transmission from I_icu2:
    event_prob[TRANS_FROM_REHAB] =
        beta_scaled * AllPops.S_pop[this_pop][1] * AllPops.I_icu2_pop[this_pop][1] * Params->frac_beta_hosp;

    //// LATENCY:
    // Exposed to infectious:
    event_prob[EXPOSED] =
        Params->delta * AllPops.E_pop[this_pop][1];

    //// RECOVERIES/TRANSITIONS:
    // Recovery asym:
    event_prob[RECOVERY_ASYM] =
        Params->recov_a * AllPops.I_asym_pop[this_pop][1];
    // Presym to sym:
    event_prob[PRESYM_TO_SYM] =
        Params->recov_p * AllPops.I_presym_pop[this_pop][1];
    // Sym to home, hosp, or ICU:
    event_prob[SYM_TO_HOME_HOSP_OR_ICU] =
        Params->recov_s * AllPops.I_sym_pop[this_pop][1];
    // Recov home
    event_prob[RECOV_HOME] =
        Params->recov_home * AllPops.I_home_pop[this_pop][1];
    // Recov hosp or move to ICU
    event_prob[RECOV_HOSP_OR_MOVE_TO_ICU] =
        Params->recov_hosp * AllPops.I_hosp_pop[this_pop][1];
    // ICU1 to ICU2(rehab) or Death
    event_prob[ICU1_TO_REHAB_OR_DEATH] =
        Params->recov_icu1 * AllPops.I_icu1_pop[this_pop][1];
    // Recov ICU2
    event_prob[RECOV_ICU2] =
        Params->recov_icu2 * AllPops.I_icu2_pop[this_pop][1];


    //// MIGRATION:
    // Migrate: S
    event_prob[MIGRATE_S] = Params->m * AllPops.S_pop[this_pop][1];
    // Migrate: I_asym
    event_prob[MIGRATE_I_ASYM] = Params->m * AllPops.I_asym_pop[this_pop][1];
    // Migrate: I_presym
    event_prob[MIGRATE_I_PRESYM] = Params->m * AllPops.I_presym_pop[this_pop][1];

    // How many events of each type will occur over time period tau?

    for(i=0;i<n_events;i++){
        // Draw Poisson random variate:
        n_occur[i] = gsl_ran_poisson(rand1, event_prob[i]*Params->tau);
    }

    delete[] event_prob;
}

/* ***************************************************  */
/* ***************************************************  */
/* ***************************************************  */

void update_pops
    (int *update_vec, int* n_occur, int this_pop,
     int this_time, int this_seed,
     COVID19ParamStruct *Params, COVID19PopStruct AllPops, gsl_rng *rand1,
     int *out_events_pos,
     int *out_events_sym,
     int *out_events_total_hosp,
     int *out_events_total_icu,
     int *out_events_n_death
    ){

    int total_exposed;
    int n_sym, n_asym, n_home, n_sym_to_icu;
    int n_hosp, n_icu1, n_hosp_recov;
    int n_icu2, n_death;
    int total_hosp, total_icu, total_pos;

    // # 1 Suscept
    total_exposed = n_occur[TRANS_FROM_ASYM] + n_occur[TRANS_FROM_PRESYM] + n_occur[TRANS_FROM_SYM] + n_occur[TRANS_FROM_HOME] +
        n_occur[TRANS_FROM_HOSP] + n_occur[TRANS_FROM_ICU1] + n_occur[TRANS_FROM_REHAB];

    update_vec[0] =
        AllPops.S_pop[this_pop][1] - total_exposed;

    if(update_vec[0] < 0){
        total_exposed = total_exposed + update_vec[0]; // adds a negative
        update_vec[0] = 0;
    }

    // # 2 Exposed
    update_vec[1] =
        AllPops.E_pop[this_pop][1] + total_exposed - n_occur[EXPOSED];
    if(update_vec[1] < 0){
        n_occur[EXPOSED] = n_occur[EXPOSED] + update_vec[1];
        update_vec[1] = 0;
    }

    // # 3 I_asym,
    n_sym = gsl_ran_binomial(rand1, (1.0 - Params->asym_rate), n_occur[EXPOSED]);
    //n_sym = (int) ( (1.0 - Params->asym_rate) * n_occur[EXPOSED] + 0.5);
    n_asym = (n_occur[EXPOSED] - n_sym);
    update_vec[2] =
        AllPops.I_asym_pop[this_pop][1] + n_asym - n_occur[RECOVERY_ASYM];
    if(update_vec[2] < 0){
        n_occur[RECOVERY_ASYM] = n_occur[RECOVERY_ASYM] + update_vec[2];
        update_vec[2] = 0;
    }

    // # 4 I_presym,
    update_vec[3] =
        AllPops.I_presym_pop[this_pop][1] + n_sym - n_occur[PRESYM_TO_SYM];
    if(update_vec[3] < 0){
        n_occur[PRESYM_TO_SYM] = n_occur[PRESYM_TO_SYM] + update_vec[3];
        update_vec[3] = 0;
    }

    // # 5 I_sym
    update_vec[4] =
        AllPops.I_sym_pop[this_pop][1] + n_occur[PRESYM_TO_SYM] - n_occur[SYM_TO_HOME_HOSP_OR_ICU];
    if(update_vec[4] < 0){
        n_occur[SYM_TO_HOME_HOSP_OR_ICU] = n_occur[SYM_TO_HOME_HOSP_OR_ICU] + update_vec[4];
        update_vec[4] = 0;
    }

    // # 6 I_home
    n_hosp = gsl_ran_binomial(rand1, Params->hosp_rate, n_occur[SYM_TO_HOME_HOSP_OR_ICU]);
    //n_hosp = (int) (Params->hosp_rate * n_occur[SYM_TO_HOME_HOSP_OR_ICU] + 0.5);
    n_sym_to_icu = gsl_ran_binomial(rand1, Params->sym_to_icu_rate, n_occur[SYM_TO_HOME_HOSP_OR_ICU]);
    //n_sym_to_icu = (int) (Params->sym_to_icu_rate * n_occur[SYM_TO_HOME_HOSP_OR_ICU] + 0.5);
    n_home = n_occur[SYM_TO_HOME_HOSP_OR_ICU] - n_hosp - n_sym_to_icu;

    update_vec[5] =
        AllPops.I_home_pop[this_pop][1] + n_home - n_occur[RECOV_HOME];
    if(update_vec[5] < 0){
        n_occur[RECOV_HOME] = n_occur[RECOV_HOME] + update_vec[5];
        update_vec[5] = 0;
    }

    // # 7 I_hosp
    update_vec[6] =
        AllPops.I_hosp_pop[this_pop][1] + n_hosp - n_occur[RECOV_HOSP_OR_MOVE_TO_ICU];
    if(update_vec[6] < 0){
        n_occur[RECOV_HOSP_OR_MOVE_TO_ICU] = n_occur[RECOV_HOSP_OR_MOVE_TO_ICU] + update_vec[6];
        update_vec[6] = 0;
    }

    // # 8 I_icu1
    n_icu1 = gsl_ran_binomial(rand1, Params->icu_rate, n_occur[RECOV_HOSP_OR_MOVE_TO_ICU]);
    //n_icu1 = (int) (Params->icu_rate * n_occur[RECOV_HOSP_OR_MOVE_TO_ICU] + 0.5);
    n_hosp_recov = n_occur[RECOV_HOSP_OR_MOVE_TO_ICU] - n_icu1;

    update_vec[7] =
        AllPops.I_icu1_pop[this_pop][1] + n_sym_to_icu + n_icu1 - n_occur[ICU1_TO_REHAB_OR_DEATH];
    if(update_vec[7] < 0){
        n_occur[ICU1_TO_REHAB_OR_DEATH] = n_occur[ICU1_TO_REHAB_OR_DEATH] + update_vec[7];
        update_vec[7] = 0;
    }

    // # 9 I_icu2
    n_icu2 = gsl_ran_binomial(rand1, (1 - Params->death_rate), n_occur[ICU1_TO_REHAB_OR_DEATH]);
    //n_icu2 = (int) ( (1 - Params->death_rate) * n_occur[ICU1_TO_REHAB_OR_DEATH] + 0.5 );
    n_death = n_occur[ICU1_TO_REHAB_OR_DEATH] - n_icu2;
    //n_icu2 = n_occur[ICU1_TO_REHAB_OR_DEATH] - n_death;

    update_vec[8] =
        AllPops.I_icu2_pop[this_pop][1] + n_icu2 - n_occur[RECOV_ICU2];
    if(update_vec[8] < 0){
        n_occur[RECOV_ICU2] = n_occur[RECOV_ICU2] + update_vec[8];
        update_vec[8] = 0;
    }

    // # 10 Recov
    update_vec[9] =
        AllPops.R_pop[this_pop][1] + n_occur[RECOVERY_ASYM] + n_occur[RECOV_HOME] + n_hosp_recov + n_occur[RECOV_ICU2];

    // # 11 Deaths
    update_vec[10] =
        AllPops.D_pop[this_pop][1] + n_death;

    total_pos = (n_sym + n_asym);
    total_hosp = (n_hosp + n_sym_to_icu);
    total_icu = (n_sym_to_icu + n_icu1);

    // SAVE NEW EVENTS
    *out_events_pos = total_pos;
    *out_events_sym = n_occur[PRESYM_TO_SYM];
    *out_events_total_hosp = total_hosp;
    *out_events_total_icu = total_icu;
    *out_events_n_death = n_death;
}

/* ***************************************************  */
/* ***************************************************  */
/* ***************************************************  */

void move_pops(int *update_vec_move, int n_occur, int this_pop,
               COVID19ParamStruct *Params, gsl_rng *rand1){

    int n_pop = Params->n_pop;
    int i;
    int this_total;
    unsigned int *temp_vec = new unsigned int[n_pop];

    //temp variables
    double *temp_prob = new double[(n_pop)];

    for(i=0;i<n_pop;i++){
        temp_prob[i] = Params->prob_move[(i+1)][this_pop];
        // temp_prob[i] = Params->dist_mat[(i+1)][this_pop];

        //printf("temp_prob[%d][%d]: %f\n", i+1, this_pop, temp_prob[i]);
    }

    // How many will travel to each population:
    //*****************************************

    //:
    this_total = n_occur;
    // printf("this_total: %d\n", this_total);
    // Choose which pops...

    // printf("Up to multinom\n");
    gsl_ran_multinomial(rand1, n_pop, this_total, temp_prob, temp_vec);

    //*****************************************
    for(i=0;i<n_pop;i++){
        update_vec_move[i] = temp_vec[i];
    }

    delete[] temp_vec;
    delete[] temp_prob;
}

/* ***************************************************  */
/* ***************************************************  */
/* ***************************************************  */

void update_pop_migrants(int *update_vec_migrants, int this_pop,
                         COVID19ParamStruct *Params, COVID19PopStruct AllPops,
                         COVID19MoveMatStruct MovePops, gsl_rng *rand1){

    int i, other_pop;
    int n_pop = Params->n_pop;
    int S_move_sum = 0;
    int S_remain = 0;
    int n_infect_visitors;
    double n_visitors, infect_frac;

    double event_prob = 0.0;
    int n_event = 0;

    double infect_sum  = 0.0;
    double beta_scaled;

    update_vec_migrants[0] = AllPops.S_pop[this_pop][2];
    update_vec_migrants[1] = AllPops.E_pop[this_pop][2];
    update_vec_migrants[2] = AllPops.I_asym_pop[this_pop][2];
    update_vec_migrants[3] = AllPops.I_presym_pop[this_pop][2];
    update_vec_migrants[4] = AllPops.I_sym_pop[this_pop][2]; // No Change
    update_vec_migrants[5] = AllPops.I_home_pop[this_pop][2]; // No Change
    update_vec_migrants[6] = AllPops.I_hosp_pop[this_pop][2]; // No Change
    update_vec_migrants[7] = AllPops.I_icu1_pop[this_pop][2]; // No Change
    update_vec_migrants[8] = AllPops.I_icu2_pop[this_pop][2]; // No Change
    update_vec_migrants[9] = AllPops.R_pop[this_pop][2]; // No Change
    update_vec_migrants[10] = AllPops.D_pop[this_pop][2]; // No Change

    // How many moved out of the pop?
    // These should not be able to get exposed by I_asym imports
    S_move_sum = 0;
    for(i=1;i<=n_pop;i++){
        other_pop = i;
        S_move_sum = S_move_sum + MovePops.S_move[other_pop][this_pop];
    }
    S_remain = AllPops.S_pop[this_pop][2] - S_move_sum;

    for(i=1;i<=n_pop;i++){

        other_pop = i;

        if(other_pop != this_pop){

            /* ******************** */
            // Suscept move TO other pop and get exposed...
            /* ******************** */

            //// sum of infectious:
            //// not visiting hospital...
            infect_sum =
            AllPops.I_asym_pop[other_pop][2] +
            AllPops.I_presym_pop[other_pop][2] +
            AllPops.I_sym_pop[other_pop][2];

            if(infect_sum > 0.0){

                trans_type_beta(beta_scaled,
                                other_pop, // NEEDS TO BE *OTHER* POP
                                infect_sum,
                                Params,
                                AllPops,
                                rand1);

            }else{
                beta_scaled = 0.0;
            }

            // The sucept that moved get exposed
            event_prob =
                beta_scaled * MovePops.S_move[other_pop][this_pop] *
                (AllPops.I_asym_pop[other_pop][2]*Params->frac_beta_asym +
                    AllPops.I_presym_pop[other_pop][2] + AllPops.I_sym_pop[other_pop][2]);

            n_event = gsl_ran_poisson(rand1, event_prob*Params->tau);

            // update the S and E:
            update_vec_migrants[0] = update_vec_migrants[0] - n_event;
            update_vec_migrants[1] = update_vec_migrants[1] + n_event;


            /* ******************** */
            // I move FROM other pop and CAUSE exposure in focal pop...
            /* ******************** */
            event_prob = 0;
            n_event = 0;

            //// sum of infectious from OTHER pop:
            infect_sum = MovePops.I_move[this_pop][other_pop];

            if(infect_sum > 0.0){

                trans_type_beta(beta_scaled,
                                this_pop, // NEEDS TO BE *THIS* POP
                                infect_sum,
                                Params,
                                AllPops,
                                rand1);

            }else{
                beta_scaled = 0.0;
            }

            // The remaining sucept that moved get exposed
            event_prob = beta_scaled * S_remain * infect_sum * Params->frac_beta_asym;

            n_event = gsl_ran_poisson(rand1, event_prob*Params->tau);

            // update the S and E:
            update_vec_migrants[0] = update_vec_migrants[0] - n_event;
            update_vec_migrants[1] = update_vec_migrants[1] + n_event;

        }// end other_pop != this_pop

    }//end i n_pop


    /* ******************** */
    // DO "TOURISTS" CAUSE INFECTION?
    /* ******************** */

    n_visitors = Params->imm_frac * Params->pop_N[this_pop];

    infect_frac =
        (AllPops.I_asym_pop[this_pop][2] +
        AllPops.I_presym_pop[this_pop][2] +
        AllPops.I_sym_pop[this_pop][2]) / Params->pop_N[this_pop];

    n_infect_visitors = gsl_ran_poisson(rand1, infect_frac * n_visitors);

    if(n_infect_visitors > 0){

        infect_sum = n_infect_visitors;

        trans_type_beta(beta_scaled,
                        this_pop,
                        infect_sum,
                        Params,
                        AllPops,
                        rand1);

        event_prob = beta_scaled * S_remain * infect_sum * Params->frac_beta_asym;
        n_event = gsl_ran_poisson(rand1, event_prob*Params->tau);

        // update the S and E:
        update_vec_migrants[0] = update_vec_migrants[0] - n_event;
        update_vec_migrants[1] = update_vec_migrants[1] + n_event;
    }


}


/*
 * Function: covid19_beta_calc
 *
 * Used by GSL root solver in calculateBeta.
 */
double covid19_beta_calc (double beta, void *params)
{
    double result, temp;
    struct covid_beta_calc_struct *p = (struct covid_beta_calc_struct *)params;

    result = (beta * p->Params->frac_beta_asym * p->Params->asym_rate / p->Params->recov_a);

    temp = (beta / p->Params->recov_p);
    temp += (beta / p->Params->recov_s);
    temp += (p->Params->sym_to_icu_rate * (beta * p->Params->frac_beta_hosp / p->Params->recov_icu1));
    temp += ((1 - p->Params->hosp_rate) * (beta / p->Params->recov_home));
    temp += p->Params->hosp_rate * (beta * p->Params->frac_beta_hosp / p->Params->recov_hosp);
    temp += p->Params->hosp_rate * (p->Params->icu_rate * (beta * p->Params->frac_beta_hosp / p->Params->recov_icu1));
    temp += p->Params->hosp_rate * (p->Params->icu_rate * ((1 - p->Params->death_rate) * (beta * p->Params->frac_beta_hosp / p->Params->recov_icu2)));

    result += (1 - p->Params->asym_rate) * temp;
    result -= p->r0;

    return result;
}


/*
 * Function: calculateBeta
 *
 * Calculates beta based on R0 value.
 */
double calculateBeta(float r0, COVID19ParamStruct *Params)
{
    int status;
    int iteration = 0;
    int max_iter = MAX_BRENT_ITERATIONS;

    double root = 0;
    double root_low = BETA_LOWER_LIMIT;
    double root_high = BETA_UPPER_LIMIT;

    const gsl_root_fsolver_type *root_fsolver_type;
    root_fsolver_type = gsl_root_fsolver_brent;

    gsl_root_fsolver *s;
    s = gsl_root_fsolver_alloc(root_fsolver_type);

    struct covid_beta_calc_struct params = {r0, Params};

    gsl_function F;
    F.function = &covid19_beta_calc;
    F.params = &params;

    gsl_root_fsolver_set (s, &F, root_low, root_high);

    // For debugging output in console, set the constant to 1 at the top of this file.
    if (OUTPUT_DEBUG_ON == 1)
    {
        printf ("using %s method\n", gsl_root_fsolver_name (s));

        printf ("%5s [%9s, %9s] %9s %9s\n",
                "iter", "lower", "upper", "root",
                "err(est)");
    }

    /*
     * Use GSL root solver to find the root. Number of iterations to try and
     * the upper and lower bounds for beta are set as constants at the top
     * of this file.
     */
    do
    {
        iteration++;
        status = gsl_root_fsolver_iterate (s);
        root = gsl_root_fsolver_root(s);
        root_low = gsl_root_fsolver_x_lower(s);
        root_high = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(root_low, root_high, 0, 0.001);

        // For debugging output in console, set the constant to 1 at the top of this file.
        if (OUTPUT_DEBUG_ON == 1)
        {
            if (status == GSL_SUCCESS) printf ("Converged:\n");

            printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
                    iteration, root_low, root_high,
                    root ,
                    root_high - root_low);
        }
    }
    while (status == GSL_CONTINUE && iteration < max_iter);

    // For debugging output in console, set the constant to 1 at the top of this file.
    if (OUTPUT_DEBUG_ON == 1)
    {
        printf("\tDEBUG: beta = %.4f for R0 = %.1f\n\n", root, r0);
    }

    return root;
}
