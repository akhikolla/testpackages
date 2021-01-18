/*
 fstcore - R bindings to the fstlib library

 Copyright (C) 2017-present, Mark AJ Klik

 This file is part of the fstcore R package.

 The fstcore R package is free software: you can redistribute it and/or modify it
 under the terms of the GNU Affero General Public License version 3 as
 published by the Free Software Foundation.

 The fstcore R package is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License
 for more details.

 You should have received a copy of the GNU Affero General Public License along
 with the fstcore R package. If not, see <http://www.gnu.org/licenses/>.

 You can contact the author at:
 - fstcore R package source repository : https://github.com/fstpackage/fstcore
*/


#ifndef TYPE_FACTORY_H
#define TYPE_FACTORY_H


#include <interface/itypefactory.h>


class BlobContainer : public IBlobContainer
{
  unsigned char* raw_vecp;
  unsigned long long length; 

public:

  BlobContainer(unsigned long long size, const SEXP r_container)
  {
    // PROTECT raw vector by containing it in a list object
    SEXP raw_vec = PROTECT(Rf_allocVector(RAWSXP, size));
    SET_VECTOR_ELT(r_container, 0, raw_vec);
    
    this->raw_vecp = RAW(raw_vec);
    this->length = Rf_xlength(raw_vec);

    UNPROTECT(1);
  }

  ~BlobContainer()
  {
  }

  unsigned char* Data()
  {
    return raw_vecp;
  }

  unsigned long long Size()
  {
    return length;
  }
};


class TypeFactory : public ITypeFactory
{
  SEXP r_container;

public:

  TypeFactory(const SEXP r_container)
  {
    this->r_container = r_container;
  }

  ~TypeFactory() { }

  /**
   * \brief Create a BlobContainer type of size indicated.
   * \param size Size of BlobContainer to create.
   * \return Pointer to the generated BlobContainer object;
   */
  IBlobContainer* CreateBlobContainer(unsigned long long size)
  {
    return new BlobContainer(size, r_container);
  }
};


#endif  // TYPE_FACTORY_H
