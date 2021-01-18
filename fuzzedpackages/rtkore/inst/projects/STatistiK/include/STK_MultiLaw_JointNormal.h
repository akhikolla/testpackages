/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

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

 Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 */

/*
 * Project:  stkpp::STatistiK::MultiLaw
 * created on: 23 janv. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_MultiLaw_JointNormal.h
 *  @brief In this file we define the joint Normal probability law.
 **/

#ifndef STK_MULTILAW_JOINTNORMAL_H
#define STK_MULTILAW_JOINTNORMAL_H

#include "STK_MultiLaw_IMultiLaw.h"
#include "STK_Law_Normal.h"

namespace STK
{

namespace MultiLaw
{
/** @ingroup Laws
 * @brief Implement the Joint Normal probability law
 **/
template<class RowVector>
class JointNormal: public JointProbability<RowVector, Law::Normal>
{
  public:
    typedef JointProbability<RowVector, Law::Normal> Base;
    /** default constructor */
    JointNormal(): Base(_T("JointNormal")) {}
    /** constructor */
    JointNormal( int const& nbVar): Base(_T("JointNormal"), nbVar) {}
    /** destructor */
    virtual ~JointNormal() {}
    /** @return the mean of the jth law */
    inline Real const& mu(int j) const { return this->law(j).mu();}
    /** set the mean of the jth law */
    inline void setMu(int j, Real const& mu)
    { this->law(j).setMu(mu);}
    /** @return the standard deviation of the jth law */
    inline Real const& sigma(int j) const { return this->law(j).sigma();}
    /** set the standard deviation of the jth law */
    inline void setSigma(int j, Real const& sigma)
    { this->law(j).setSigma(sigma);}
};

} // namespace Law

} // namespace STK

#endif /* STK_MULTILAW_JOINTNORMAL_H */
