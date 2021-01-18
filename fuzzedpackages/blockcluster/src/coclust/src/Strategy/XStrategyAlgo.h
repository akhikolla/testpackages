/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2015  <MODAL team @INRIA,Lille & U.M.R. C.N.R.S. 6599 Heudiasyc, UTC>

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this program; if not, write to the
 Free Software Foundation, Inc.,
 59 Temple Place,
 Suite 330,
 Boston, MA 02111-1307
 USA

 Author : Parmeet Bhatia
 Contact : bhatia.parmeet@gmail.com , serge.iovleff@stkpp.org

 */


/** @file XStrategyAlgo.h
 *  @brief Declares XStrategyAlgo class.
 **/

#ifndef XSTRATEGYALGO_H_
#define XSTRATEGYALGO_H_


#include "IStrategy.h"

/** @brief This class provides functionality to run various algorithms except SEMGibbsAlgo.
 *
 */
class XStrategyAlgo: public IStrategy
{
  public:
    XStrategyAlgo(StrategyParameters & Aparam)
    :IStrategy(Aparam)
    {};
    virtual bool run();
    virtual ~XStrategyAlgo(){};
  private:
    void writeVector(int it, MatrixReal classes, std::string type);
};

#endif /* XSTRATEGYALGO_H_ */
