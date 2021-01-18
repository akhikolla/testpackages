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

/** @file STK_KernelHandler.h
 *  @brief In this file we define the KernelHandler class.
 **/


#ifndef STK_KERNELHANDLER_H
#define STK_KERNELHANDLER_H

#include <DManager/include/STK_DataHandlerBase.h>
#include <Arrays/include/STK_CArray.h>
#include <Arrays/include/STK_Array1D.h>

#include <STatistiK/include/STK_Kernel_Util.h>
#include <STatistiK/include/STK_Kernel_Gaussian.h>
#include <STatistiK/include/STK_Kernel_Hamming.h>
#include <STatistiK/include/STK_Kernel_Laplace.h>
#include <STatistiK/include/STK_Kernel_Linear.h>
#include <STatistiK/include/STK_Kernel_Polynomial.h>
#include <STatistiK/include/STK_Kernel_RationalQuadratic.h>


namespace STK
{
// forward declaration
class KernelHandler;

namespace hidden
{
/** @ingroup hidden
 *  Specialization of the DataHandlerTraits struct for KernelHandler
 **/
template<typename Type>
struct DataHandlerTraits<KernelHandler, Type>
{
  typedef CArray<Type> Data;
};

} // namespace hidden

/** @ingroup DManager
 *  @c implementation of the DataHandlerBase class for kernel mixture models.
 *  This class stores the kernels and the idData, allow to add them and to remvoe
 *  them.
 */
class KernelHandler: public DataHandlerBase<KernelHandler>
{
  public:
    typedef DataHandlerBase<KernelHandler> Base;
    typedef DataHandlerBase<KernelHandler>::InfoMap InfoMap;
    typedef std::pair<Kernel::IKernel*, String> TaggedKernel;
    typedef Array1D< TaggedKernel >::iterator Iterator;
    typedef Array1D< TaggedKernel >::const_iterator ConstIterator;

    /** default constructor */
    KernelHandler();
    /** constructor with an instance of a kernel to handle*/
    KernelHandler(Kernel::IKernel* p_kernel, String const& idData, String const& idModel);
    /** destructor. All pointer on kernels will be deleted. */
    ~KernelHandler();
    /** @return the number of sample (the number of rows of the data) */
    inline int nbSample() const { return nbSample_;}
    /** add an instance of a kernel to the handler
     * @param p_kernel an instance of a kernel
     * @param idData can be any string given by the user for identifying data.
     * @param idModel represent the idModel of a given model (can be defined
     *  inside or outside STK++).
     * @return @c true if the pair (idData,idModel) has been successfully added
     * to the handler
     **/
    bool addKernel(Kernel::IKernel* p_kernel, String const& idData, String const& idModel);
    /** add an instance of a kernel to the handler
     * @param kernelName name of the kernel
     * @param data data set associated to the kernel
     * @param param parameters of the kernel
     * @param idData can be any string given by the user for identifying data.
     * @param idModel represent the idModel of a given model (can be defined
     *  inside or outside STK++).
     * @return @c true if the pair (idData,idModel) has been successfully added
     * to the handler
     **/
    template<class Derived, class Array>
    bool addKernel( String const& kernelName
                  , Array const& data
                  , ExprBase<Derived> const& param
                  , String const& idData
                  , String const& idModel);
    /** get an instance of a kernel from the handler
     * @param idData can be any string given by the user for identifying data.
     * @return @c 0 if the idData has not been found, a pointer to the kernel
     * otherwise
     */
    Kernel::IKernel const* getKernel( String const& idData) const;
    /** remove an instance of a kernel to the handler
     *  @param idData can be any string given by the user for identifying data.
     */
    void removeKernel( String const& idData);
    /** utility lookup function allowing to know if some pointer on a kernel
     *  is handled by the KernelHandler.
     *  @return @c true if the pointed kernel is found, @c false otherwise
     **/
    bool isHandled(Kernel::IKernel* const p_kernel) const;

  protected:
    /** Number of sample */
    int nbSample_;
    /** Array of the kernels */
    Array1D< TaggedKernel > v_kernel_;
};

template<class Derived, class Array>
bool KernelHandler::addKernel( String const& kernelName
                             , Array const& data
                             , ExprBase<Derived> const& param
                             , String const& idData
                             , String const& idModel)
{
  // check if the idData already exists
  Array1D< TaggedKernel >::Iterator it;
  for (it = v_kernel_.beginIterator() ; it!=v_kernel_.endIterator(); ++it)
  { if (it->second == idData) { break;}}
  // if idData already exists return false
  if (it != v_kernel_.endIterator()) return false;

  // create kernel
  Kernel::kernelType kt = Kernel::stringToKernelType(kernelName);
  Kernel::IKernel* p_kernel = 0;
  switch (kt)
  {
    case Kernel::laplace_:
      p_kernel = new Kernel::Laplace<Array>(data, param);
      break;
    case Kernel::gaussian_:
      p_kernel = new Kernel::Gaussian<Array>(data, param);
      break;
    case Kernel::hamming_:
      p_kernel = new Kernel::Hamming<Array>(data, param);
      break;
    case Kernel::linear_:
      p_kernel = new Kernel::Linear<Array>(data, param);
      break;
    case Kernel::polynomial_:
      p_kernel = new Kernel::Polynomial<Array>(data, param);
      break;
    case Kernel::rationalQuadratic_:
      p_kernel = new Kernel::RationalQuadratic<Array>(data, param);
      break;
    default:
      break;
  }
  // not an existing kernel
  if (!p_kernel) return false;
  // add kernel to
  v_kernel_.push_back(TaggedKernel(p_kernel, idData));
  nbSample_   = p_kernel->nbSample();
  return addInfo(idData, idModel);
}

} // namespace STK

#endif /* STK_KERNELHANDLER_H */
