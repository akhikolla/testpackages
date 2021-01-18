/*
	xmlfile.cpp
	Implementation of InputStream (only) that uses Miniz zip extract iter.

	Copyright (c) 2019 Ibrahim Umar
*/

#include "xmlio/xmlzipfile.h"
#include "miniz/miniz.h"
#include <assert.h>

XML_BEGIN_NAMESPACE

//
// ZipInputStream
//

ZipInputStream::ZipInputStream(const char *path, const char *filename)
{
	memset(&zip_archive, 0, sizeof(zip_archive));

	status = mz_zip_reader_init_file(&zip_archive, path, 0);
	if (!status)
		throw ZipException(0);

	zipstate = mz_zip_reader_extract_file_iter_new(&zip_archive, filename, 0);

	if(!zipstate)
		throw ZipException(1);
}

ZipInputStream::~ZipInputStream()
{
	if(zipstate)
		status = mz_zip_reader_extract_iter_free(zipstate);
	mz_zip_reader_end(&zip_archive);
}

int ZipInputStream::read(XML_Char *buf, size_t bufLen)
{
	assert(buf);
	return mz_zip_reader_extract_iter_read(zipstate, buf, bufLen);
}

//
// ZipException
//

const char * ZipException::what() const throw ()
{
	static const char* const errStr[] = {
		"mz_zip_reader_init_file() failed!\n",
		"Unable to find the file inside zip archive!\n"
	};

	return errStr[mErrCode];
}

XML_END_NAMESPACE
