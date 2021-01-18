/*
	xmlzipfile.h
	Implementation of InputStream (only) that uses Miniz zip extract iter.

	Copyright (c) 2019 Ibrahim Umar
*/

#ifndef XMLZIPFILE_H
#define XMLZIPFILE_H

#include "xmlstream.h"
#include "../miniz/miniz.h"

#include <exception>

#ifdef __cplusplus
//
// C++ Implementation
//
XML_BEGIN_NAMESPACE

/**
	Implementation of InputStream interface that uses miniz's zip_extract_iter
*/
class ZipInputStream : public InputStream
{
public:
	/** Constructor - open the file for reading in text mode
		@param path		path to the file to open
		@throws			FileException if the file cannot be opened
		@see			FileException
	*/
	ZipInputStream(const char *path, const char *filename);
	virtual ~ZipInputStream();
	virtual int read(XML_Char *buf, size_t bufLen);
private:
	int errNo;
	mz_bool status;
	mz_zip_archive zip_archive;
	mz_zip_reader_extract_iter_state* zipstate;
};

/**
	Exception thrown by ZipInputStream if anything goes wrong.
*/
class ZipException : public std::exception
{
public:
	ZipException(int errCode) : mErrCode(errCode) {}
	int GetErrorCode() const { return mErrCode; }
	const char * what () const throw ();

private:
	int mErrCode;
};

XML_END_NAMESPACE
#endif	// __cplusplus

#endif	// XMLZIPFILE_H

