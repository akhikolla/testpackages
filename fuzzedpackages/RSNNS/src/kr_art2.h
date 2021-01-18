/************************************************************************************

   This file is part of SnnsCLib, a fork of the kernel and parts of the gui of 
   the Stuttgart Neural Network Simulator (SNNS), version 4.3.

   The file's original version is part of SNNS 4.3. It's source code can be found at

   http://www.ra.cs.uni-tuebingen.de/SNNS/

   SNNS 4.3 is under the license LGPL v2. We note that source code files of SNNS 4.3 
   state as version "4.2". Base of this fork is SNNS 4.3 with a reverse-applied 
   python patch (see http://developer.berlios.de/projects/snns-dev/).

   SnnsCLib was developed in 2010 by Christoph Bergmeir under supervision of 
   José M. Benítez, both affiliated to DiCITS Lab, Sci2s group, DECSAI, 
   University of Granada

   Changes done to the original code were performed with the objective to
   port it from C to C++ and to encapsulate all code in one class named SnnsCLib.

   Changes in header files mainly include:
   * removed all static keywords
   * moved initializations of variables to the constructor of SnnsCLib

   Changes in cpp code files mainly include:
   * changed file ending from .c to .cpp
   * removed all SNNS internal includes and only include SnnsCLib   
   * static variables within functions were turned into member variables of SnnsCLib
   * function declarations were changed to method declarations, i.e. "SnnsCLib::.."
     was added
   * calls to the function table are now "C++-style", using the "this"-pointer

   License of SnnsCLib:
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.
 
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.
 
   You should have received a copy of the GNU Library General Public License
   along with this library; see the file COPYING.LIB.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.

************************************************************************************/


/*****************************************************************************
  FILE           : $Source: /projects/higgs1/SNNS/CVS/SNNS/kernel/sources/kr_art2.ph,v $
  SHORTNAME      : kr_art2 
  SNNS VERSION   : 4.2

  PURPOSE        : SNNS Kernel Function Prototypes for ART2-Networks
  NOTES          :

  AUTHOR         : Kai-Uwe Herrmann
  DATE           : 17.05.92

  CHANGED BY     : Sven Doering
  RCS VERSION    : $Revision: 2.6 $
  LAST CHANGE    : $Date: 1998/02/25 15:26:38 $

    Copyright (c) 1990-1995  SNNS Group, IPVR, Univ. Stuttgart, FRG
    Copyright (c) 1996-1998  SNNS Group, WSI, Univ. Tuebingen, FRG

******************************************************************************/
#ifndef _KR_ART2_DEFINED_
#define  _KR_ART2_DEFINED_

/* begin global definition section */

/************* Global variables
*************/
 int              Art2_NoOfRecUnits;
 struct Unit     *Art2_cl_unit;
 struct Unit     *Art2_nc_unit;

/**************** Function Prototypes
****************/

/***************************************************************************/
/* kra2_init_propagate ()

   initializes net for propagation.
*/
 krui_err  kra2_init_propagate (
                                      void
                                     );


/***************************************************************************/
/* kra2_sort ()

   Set logical layer numbers and logical unit numbers in an ART2 network.
   Also this function checks, whether the network is an ART2 network or not.
   Returns an error code, when actual network is no ART2 architecture.
*/
 krui_err  kra2_sort (
                            void
                           );



/***************************************************************************/
/* kra2_set_params ()

   Sets the value of Parameters rho, a, b, c, d, theta, which are stored locally
   in this Module.
*/
 krui_err kra2_set_params (
                                FlintType rho,
                                FlintType param_a,
                                FlintType param_b,
                                FlintType param_c,
                                FlintType param_d,
                                FlintType theta
                               );


/***************************************************************************/
/* kra2_get_rho ()

   returns the actual value of Parameter rho.
*/
 FlintType kra2_get_rho (
                               void
                              );

/***************************************************************************/
/* kra2_get_a ()

   returns the actual value of Parameter a.
*/
 FlintType kra2_get_a (
                             void
                            );

/***************************************************************************/
/* kra2_get_b ()

   returns the actual value of Parameter b.
*/
 FlintType kra2_get_b (
                             void
                            );

/***************************************************************************/
/* kra2_get_c ()

   returns the actual value of Parameter c.
*/
 FlintType kra2_get_c (
                             void
                            );


/***************************************************************************/
/* kra2_get_d ()

   returns the actual value of Parameter d.
*/
 FlintType kra2_get_d (
                             void
                            );


/***************************************************************************/
/* kra2_get_theta ()

   returns the actual value of Parameter theta.
*/
 FlintType kra2_get_theta (
                                 void
                                );


/***************************************************************************/
/* kra2_checkReset ()

   checks if global reset has to be sent into network
*/
 void kra2_checkReset (
                             void
                            );




/***************************************************************************/
/* kra2_Reset ()

   returns TRUE if global reset is actually active, else FALSE
*/
 bool kra2_Reset (
                        void
                       );




/***************************************************************************/
/* kra2_init_pattern ()

   sets current phase to bottom up
*/
 void kra2_init_pattern (
                               void
                              );



/***************************************************************************/
/* kra2_top_dn_phase ()

   sets current phase to bottom up
*/
 bool kra2_topdn_phase (
                             void
                            );



/***************************************************************************/
/* kra2_compute_norms ()

   computes the L2 vector norms of inp, w, u, v, p, r
*/

 void kra2_compute_norms (
                                void
                               );






/***************************************************************************/
/* kra2_L2_Norm ()

   returns the L2-Norm of a vector which is determined by the number
   of the Layer in the topo_ptr_array.
*/

 FlintType kra2_L2_Norm (
                               int VectorNr
                              );



/***************************************************************************/
/* kra2_classified ()

   returns TRUE if net has classified input pattern
*/
 bool kra2_classified (
                             void
                            );



/***************************************************************************/
/* kra2_not_classifiable ()

   returns TRUE if net is not able to classify input pattern
*/
 bool kra2_not_classifiable (
                                   void
                                  );



/***************************************************************************/
/* kra2_save_for_stability_check ()

   saves informaion of relevant unit in F1-Layer for stability check
*/
 void kra2_save_for_stability_check (
                                           void
                                          );



/***************************************************************************/
/* kra2_check_f1_stability ()

   checks, if F1-Layer is stable.
*/
 void kra2_check_f1_stability (
                                     void
                                    );





/***************************************************************************/
/* kra2_f1_stable ()

   returns TRUE if F1-Layer is stable
*/
 bool kra2_f1_stable (
                            void
                           );




/***************************************************************************/
/* kra2_getClassNo ()

   Returns the number of the actually activated class J, 1 <= J <= M
*/
 int  kra2_getClassNo (void
                            );


/* end global definition section */

/* begin private definition section */


/*#################################################

GROUP: local defines

#################################################*/



#define MIN_NO_OF_DELAY_STEPS   5                /* when checking if pattern has been
                                                    classified F1 layer has to be stable
                                                    over MIN_NO_OF_DELAY_STEPS cycles
                                                 */

#define F1_STABILITY_PARAM      0.0001           /* if difference of activation in all
                                                    u units between two prop cycles are
                                                    less or equal to F1_STABILITY_PARAM
                                                    then F1 is stable
                                                 */


/*#################################################

GROUP: global variables, local to this module

#################################################*/

/* Global variable for parameter values */

FlintType      Param_rho;
FlintType      Param_a;
FlintType      Param_b;
FlintType      Param_c;
FlintType      Param_d;
FlintType      Param_theta;

/* Global variables for vector norms */

FlintType      NormInp;
FlintType      NormW;
FlintType      NormU;
FlintType      NormV;
FlintType      NormP;
FlintType      NormR;


TopoPtrArray   topo_layer[10];      /* contains pointers to first pointer
                                              to inp unit, first pointer
                                              to w unit, x unit ... in the topo
                                              pointer array
                                           */

int            NoOfDelaySteps;

bool           GlobalReset;
bool           TopDownPhase;
bool           f1_stable;   /* becomes TRUE if F1-Layer is stable
                                              (see kra2_check_f1_stability)
                                           */

/* functions that are local to this module
*/

void   kra2_set_fix_weight (

                                   struct Unit   *src_unit,
                                   struct Unit   *trgt_unit,
                                   FlintType     *weight

                                  );


FlintType  kra2_compute_l2_norm (

                                        int Layer

                                       );

int  kra2_get_NoOfRecUnits (

                                   void

                                  );


krui_err  kra2_get_InpUnits (

                                    TopoPtrArray  *topo_ptr

                                   );


krui_err  kra2_get_WUnits (

                                  TopoPtrArray *topo_ptr,
                                  int          *no_of_w_units

                                 );


krui_err  kra2_get_XUnits (

                                  TopoPtrArray *topo_ptr,
                                  int          *no_of_x_units

                                 );


krui_err  kra2_get_UUnits (

                                  TopoPtrArray *topo_ptr,
                                  int          *no_of_u_units

                                 );


krui_err  kra2_get_VUnits (

                                  TopoPtrArray *topo_ptr,
                                  int          *no_of_v_units

                                 );


krui_err  kra2_get_PUnits (

                                  TopoPtrArray *topo_ptr,
                                  int          *no_of_p_units

                                 );


krui_err  kra2_get_QUnits (

                                  TopoPtrArray *topo_ptr,
                                  int          *no_of_q_units

                                 );


krui_err  kra2_get_RUnits (

                                  TopoPtrArray *topo_ptr,
                                  int          *no_of_r_units

                                 );


krui_err  kra2_get_RecUnits (

                                    TopoPtrArray  *topo_ptr

                                   );


krui_err  kra2_get_RstUnits (

                                    TopoPtrArray  *topo_ptr,
                                    int          *no_of_rst_units

                                   );


krui_err  kra2_TopoPtrArray (

                                    void

                                   );


krui_err  kra2_LinksToInpUnits (

                                       TopoPtrArray *topo_ptr

                                      );


krui_err  kra2_LinksToWUnits (

                                     TopoPtrArray *topo_ptr

                                    );


krui_err  kra2_LinksToXUnits (

                                     TopoPtrArray *topo_ptr

                                    );


krui_err  kra2_LinksToUUnits (

                                     TopoPtrArray *topo_ptr

                                    );


krui_err  kra2_LinksToVUnits (

                                     TopoPtrArray *topo_ptr

                                    );


krui_err  kra2_LinksToPUnits (

                                     TopoPtrArray *topo_ptr

                                    );


krui_err  kra2_LinksToQUnits (

                                     TopoPtrArray *topo_ptr

                                    );


krui_err  kra2_LinksToRUnits (

                                     TopoPtrArray *topo_ptr

                                    );


krui_err  kra2_LinksToRecUnits (

                                       TopoPtrArray *topo_ptr

                                      );


krui_err  kra2_LinksToRstUnits (
                                       TopoPtrArray *topo_ptr
                                      );

krui_err kra2_init_i_act (void);
krui_err kra2_init_fix_weights (void);

/* end private definition section */

#endif 

