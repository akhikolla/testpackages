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
 * Project:  stkpp::DManager
 * created on: 15 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_KernelHandler.cpp
 *  @brief In this file we implement the KernelHandler class.
 **/


#include <Clustering/include/KernelModels/STK_KernelHandler.h>

namespace STK
{
/* default constructor */
KernelHandler::KernelHandler(): Base(), nbSample_(0) {}
/* constructor with an instance of a kernel to handle*/
KernelHandler::KernelHandler( Kernel::IKernel* p_kernel, String const& idData, String const& idModel)
                            : Base()
                            , nbSample_(p_kernel->nbSample())
{}
/* destructor */
KernelHandler::~KernelHandler()
{
  for(Iterator it = v_kernel_.beginIterator(); it != v_kernel_.endIterator(); ++it)
  {
    if (it->first)
    {
      Kernel::IKernel* p = it->first;
      it->first = 0;
      if (!isHandled(p)) delete p;
    }
  }
}

/* add an instance of a kernel to the handler
 * @param p_kernel an instance of a kernel
 * @param idData can be any string given by the user for identifying data.
 * @param idModel represent the idModel of a given model (can be defined
 *  inside or outside STK++).
 * @return @c true if the pair (idData,idModel) has been successfully added
 * to the handler
 */
bool KernelHandler::addKernel(Kernel::IKernel* p_kernel, String const& idData, String const& idModel)
{
  if (!p_kernel) return false;
  if (!addInfo(idData, idModel)) return false; // there exists an idData with an other idModel

  // check if idData already exists
  Array1D< TaggedKernel >::ConstIterator it;
  for (it = v_kernel_.beginConstIterator() ; it!=v_kernel_.endConstIterator(); ++it)
  { if (it->second == idData) { break;}}

  // add kernel if it does not exist otherwise do nothing
  if (it == v_kernel_.endConstIterator())
  {
    v_kernel_.push_back(TaggedKernel(p_kernel, idData));
    nbSample_   = p_kernel->nbSample();
  }
  return true;
}
/* get an instance of a kernel from the handler
 * @param idData can be any string given by the user for identifying data.
 * @return @c 0 if the idData has not been found, a pointer to the kernel
 * otherwise
 */
Kernel::IKernel const* KernelHandler::getKernel( String const& idData) const
{
  // check if the idData already exists
  Array1D< TaggedKernel >::ConstIterator it;
  for (it = v_kernel_.beginConstIterator() ; it!=v_kernel_.endConstIterator(); ++it)
  {  if (it->second == idData) { break;}}
  // no data set found
  if (it == v_kernel_.endConstIterator()) { return 0;}
  return it->first;
}
/* remove an instance of a kernel to the handler
 *  @param idData can be any string given by the user for identifying data.
 */
void KernelHandler::removeKernel( String const& idData)
{
  // check if the idData exists
  Array1D< TaggedKernel >::Iterator it;
  for (it = v_kernel_.beginIterator() ; it!=v_kernel_.endIterator(); ++it)
  {  if (it->second == idData) { break;}
  }
  // check if kernel has been found
  if (it != v_kernel_.endIterator())
  {
    Kernel::IKernel* p = it->first;
    it->first = 0;
    if (!isHandled(p)) delete p;
    v_kernel_.erase(it.pos(), 1);
  }
}

/* utility lookup function allowing to know if some pointer on a kernel
 *  is handled.
 *  @return @c true if the pointed kernel is found, @c false otherwise
 **/
bool KernelHandler::isHandled(Kernel::IKernel* const p_kernel) const
{
  for(ConstIterator it = v_kernel_.beginConstIterator(); it != v_kernel_.endConstIterator(); ++it)
  { if (it->first == p_kernel) return true;}
  return false;
}

} // namespace STK

