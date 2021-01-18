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
 * Purpose:  Include all files of the DManager project.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file DManager.h
 *  @brief This file include all the other header files of the
 *  project DManager.
 **/

/**
 *  @defgroup DManager Data Management
 *  @brief The DManager project propose classes and functions for managing the
 *  data.
 *
 * The DManager project contains all the class related to data management.
 * It provides
 * <ul>
 *    <li> an abstract base class for statistical variables, </li>
 *    <li> a template implementation for arbitrary data, </li>
 *    <li> a Dataframe (Table) class,</li>
 *    <li> classes for read and write csv and (TODO) dbf files, </li>
 *    <li> various  classes for importing/exporting and converting data from
 *    different containers in different type, </li>
 *    <li> methods for sorting one dimensional containers and two-dimensionnal
 *    containers using one column as sorting key, </li>
 *    <li> classes for reading/writing options file. </li>
 * </ul>
 **/

/** @ingroup DManager
 *  @namespace STK::Import
 *  @brief The import namespace contains the constants and global functions used
 *  in the DManager project when importing data in STK containers.
 **/

#ifndef DMANAGER_H
#define DMANAGER_H

/* Utilities used in the DManager Project */
#include <DManager/include/STK_DManager_Util.h>
#include <DManager/include/STK_Import_Util.h>

/* Interface Variable class and Variable class */
#include <DManager/include/STK_IVariable.h>
#include <DManager/include/STK_Variable.h>

/* DataFrame. */
#include <DManager/include/STK_Cell.h>
#include <DManager/include/STK_List1D.h>
#include <DManager/include/STK_DataFrame.h>

/* Cross Validation utilities */
#include <DManager/include/STK_CvHandler.h>
#include <DManager/include/STK_PartitionHandler.h>

/* Interface class for data handlers */
#include <DManager/include/STK_DataHandlerBase.h>
#include <DManager/include/STK_IDataHandler.h>

/* Interface class for data wrappers */
#include <DManager/include/STK_IDataBridge.h>
#include <DManager/include/STK_DataBridge.h>


/* main classes for managing Csv data. */
#include <DManager/include/STK_ReadWriteCsv.h>
#include <DManager/include/STK_DataHandler.h>

/* Export or import data from one container to an other container */
#include <DManager/include/STK_ImportFromCsv.h>
#include <DManager/include/STK_ExportToCsv.h>
#include <DManager/include/STK_DataFrameToArray2D.h>
#include <DManager/include/STK_CsvToArray.h>

/* HeapSort utilities. */
#include <DManager/include/STK_HeapSort.h>

/* main classes for option files. */
#include <DManager/include/STK_Option.h>
#include <DManager/include/STK_IPage.h>
#include <DManager/include/STK_ReadWritePages.h>

#endif  /* DMANAGER_H */
