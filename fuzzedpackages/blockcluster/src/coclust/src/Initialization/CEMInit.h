/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2015  <MODAL team @INRIA,Lille & U.M.R. C.N.R.S. 6599 Heudiasyc, UTC>

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as
 published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public
 License along with this program; if not, write to the
 Free Software Foundation, Inc.,
 59 Temple Place,
 Suite 330,
 Boston, MA 02111-1307
 USA

 Contact : parmeet.bhatia@inria.fr , bhatia.parmeet@gmail.com
 */


/** @file CEMInit.h
 *  @brief Declares CEM initialization class CEMInit derived from IInit.
 **/

#ifndef CEMINIT_H_
#define CEMINIT_H_

#include "IInit.h"

/** @brief This class provides functionalities for CEM initialization.
 *  It is  derived from IInit abstract class.
 */
class CEMInit: public IInit
{
  public:
    CEMInit(){};
    virtual bool run();
    inline virtual ~CEMInit(){};
};

inline bool CEMInit::run()
{
  return p_Model_->cemInitStep();
}

#endif /* CEMINIT_H_ */
