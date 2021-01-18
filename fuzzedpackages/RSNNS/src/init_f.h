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
  FILE           : $Source: /projects/higgs1/SNNS/CVS/SNNS/kernel/sources/init_f.ph,v $
  SHORTNAME      : init_f
  SNNS VERSION   : 4.2

  PURPOSE        : private header
  NOTES          :

  AUTHOR         : 
  DATE           : 

  CHANGED BY     : Sven Doering
  RCS VERSION    : $Revision: 2.14 $
  LAST CHANGE    : $Date: 1998/02/25 15:26:19 $

    Copyright (c) 1990-1995  SNNS Group, IPVR, Univ. Stuttgart, FRG
    Copyright (c) 1996-1998  SNNS Group, WSI, Univ. Tuebingen, FRG

******************************************************************************/
#ifndef _INIT_F_DEFINED_
#define  _INIT_F_DEFINED_

/* begin global definition section */

krui_err  INIT_randomizeWeights(float *parameterArray, int NoOfParams);
krui_err INIT_RM_randomizeWeights(float*parameterArray, int NoOfParams);
krui_err  INIT_randomizeWeights_perc(float *parameterArray, int NoOfParams);
krui_err INIT_Weights_CPNv32(float *parameterArray, int NoOfParams);
krui_err INIT_Weights_CPNv33(float *parameterArray, int NoOfParams);
krui_err INIT_Weights_CPN_Rand_Pat(float *parameterArray, int NoOfParams);
void RbfInitSetCenter(int pattern_no, int sub_pat_no, struct Unit *hidden_unit, float deviation, float bias);
void RbfInitBPCenter(struct Unit *hidden_unit);
krui_err  RbfInitNetwork(int start_pat, int end_pat, float i_bias, float i_devi, float i_f_0, float i_f_1, float i_smooth, int init_type);
#ifdef RBF_INCLUDE_KOHONEN_CONVEX

void RbfKohonenConvexInit(int start_pattern,int end_pattern,float alpha_start,
	float alpha_increment,float learn_rate,int count);
#endif
krui_err RbfKohonenInit(int start_pattern, int end_pattern, float learn_rate, int count, int shuffle);
krui_err RbfStartInit(float *parameterArray, int NoOfParams, int init_type);
krui_err INIT_RBF_Weights(float *parameterArray, int NoOfParams);
krui_err INIT_RBF_Weights_redo(float *parameterArray, int NoOfParams);
krui_err INIT_RBF_Weights_kohonen(float *parameterArray, int NoOfParams);
krui_err INIT_Weights_ART1(float *parameterArray, int NoOfParams);
krui_err INIT_Weights_ART2(float *parameterArray, int NoOfParams);
krui_err INIT_Weights_ARTMAP(float *parameterArray, int NoOfParams);
krui_err INIT_CC_Weights(float *parameterArray, int NoOfParams);
krui_err INIT_TACOMA_Weights(float *parameterArray, int NoOfParams);
krui_err INIT_SOM_Rand_Pat(float *parameterArray, int NoOfParams);
krui_err INIT_SOM_Weights_v32(float *parameterArray, int NoOfParams);
krui_err INIT_SOM_Weights_const(float *parameterArray, int NoOfParams);
krui_err INIT_JE_Weights (float *parameterArray, int NoOfParams) ;

krui_err INIT_Hebb(float *parameterArray, int NoOfParams);
krui_err INIT_ClippHebb(float *parameterArray, int NoOfParams);
krui_err INIT_HOP_FixAct(float *parameterArray, int NoOfParams);
krui_err INIT_PseudoInv(float *parameterArray, int NoOfParams);
krui_err ENZO_noinit(void);

/* end global definition section */

/* begin private definition section */

#define  INIT_PARAM1( param )   param[ 0 ]  /*    contains the 1st initialisation parameter  */
#define  INIT_PARAM2( param )   param[ 1 ]  /*    contains the 2nd initialisation parameter  */
#define  INIT_PARAM3( param )   param[ 2 ]  /*    contains the 3rd initialisation parameter  */
#define  INIT_PARAM4( param )   param[ 3 ]  /*    contains the 4th initialisation parameter  */
#define  INIT_PARAM5( param )   param[ 4 ]  /*    contains the 5th initialisation parameter  */

#define	 RBF_INIT_FULL		0
#define	 RBF_INIT_REINIT	1
#define  RBF_INIT_KOHONEN	2
#define  MY_HUGE_VAL		1E20


/* end private definition section */

krui_err   PseudoInv(RbfFloatMatrix *source, int NoOfColumns, RbfFloatMatrix *target );

#endif 
