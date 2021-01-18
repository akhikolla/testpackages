/*
	xmlfile.h
	Sample implementation of InputStream and OutputStream that uses stdio.

	Copyright (c) 2000 Paul T. Miller

	LGPL DISCLAIMER
	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Library General Public
	License as published by the Free Software Foundation; either
	version 2 of the License, or (at your option) any later version.

	This library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	Library General Public License for more details.

	You should have received a copy of the GNU Library General Public
	License along with this library; if not, write to the
	Free Software Foundation, Inc., 59 Temple Place - Suite 330,
	Boston, MA  02111-1307, USA.

	http://www.gnu.org/copyleft/lgpl.html
*/

#ifndef XMLFILE_H
#define XMLFILE_H

#include "xmlstream.h"
#include <stdio.h>

#ifdef _WIN32
#include <wchar.h>
#endif

#ifdef __cplusplus
	extern "C" {
#endif

/* "C" interface */
XML_InputStream *XML_OpenFileInputStream(const char *path);
void XML_CloseFileInputStream(XML_InputStream *stream);
XML_OutputStream *XML_OpenFileOutputStream(const char *path);
void XML_CloseFileOutputStream(XML_OutputStream *stream);

#ifdef __cplusplus
	}
#endif

#ifdef __cplusplus
//
// C++ Implementation
//
XML_BEGIN_NAMESPACE

/**
	Implementation of InputStream interface that uses stdio
*/
class FileInputStream : public InputStream
{
public:
	/** Constructor - open the file for reading in text mode
		@param path		path to the file to open
		@throws			FileException if the file cannot be opened
		@see			FileException
	*/
	FileInputStream(const char *path);
#ifdef _WIN32
	FileInputStream(const wchar_t *path);
#endif
	virtual ~FileInputStream();
	virtual int read(XML_Char *buf, size_t bufLen);
private:
	FILE *mFile;
};

/**
	Implementation of OutputStream interface that uses stdio
*/
class FileOutputStream : public OutputStream
{
public:
	/** Constructor - create the file for writing in text mode.
		@param path		path to the file to create
		@throws			FileException if the file cannot be created
		@see			FileException
	*/
	FileOutputStream(const char *path);
#ifdef _WIN32
	FileOutputStream(const wchar_t *path);
#endif
	virtual ~FileOutputStream();
	virtual int write(const char *buf, size_t bufLen);
private:
	FILE *mFile;
};

/**
	Exception thrown by FileInputStream and FileOutputStream if anything goes wrong.
*/
class FileException
{
public:
	FileException(int errCode) : mErrCode(errCode) {}
	int GetErrorCode() const { return mErrCode; }
	const char *GetErrorString() const;

private:
	int mErrCode;
};

XML_END_NAMESPACE
#endif	// __cplusplus

#endif	// XMLFILE_H

