/*
	xmlfile.cpp
	Sample implementation of InputStream and OutputStream that uses stdio.

	Copyright (c) 1999 Paul T. Miller
	Copyright (c) 1999 Artel Software, Inc.

	Permission  is  hereby  granted,  free of charge, to any person
	obtaining  a copy of this software and associated documentation
	files  (the  "Software"),  to  deal  in  the  Software  without
	restriction,  including  without  limitation the rights to use,
	copy,  modify,  merge,  publish, distribute, sublicense, and/or
	sell  copies of the Software, and to permit persons to whom the
	Software  is  furnished  to  do  so,  subject  to the following
	conditions:

	The  above copyright notice and this permission notice shall be
	included in all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
	OF  MERCHANTABILITY,  FITNESS  FOR  A  PARTICULAR  PURPOSE  AND
	NONINFRINGEMENT.  IN  NO  EVENT  SHALL THE AUTHORS OR COPYRIGHT
	HOLDERS  BE  LIABLE  FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
	WHETHER  IN  AN  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
	FROM,  OUT  OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
	OTHER DEALINGS IN THE SOFTWARE.
*/

#include "xmlio/xmlfile.h"
#include <assert.h>
#include <errno.h>
#include <string.h>

XML_BEGIN_NAMESPACE

//
// FileInputStream
//

FileInputStream::FileInputStream(const char *path)
{
	mFile = fopen(path, "r");
	if (!mFile)
		throw FileException(errno);
}

#ifdef _WIN32
FileInputStream::FileInputStream(const wchar_t *path)
{
	mFile = _wfopen(path, L"r");
	if (!mFile)
		throw FileException(errno);
}
#endif

FileInputStream::~FileInputStream()
{
	fclose(mFile);
}

int FileInputStream::read(XML_Char *buf, size_t bufLen)
{
	assert(buf);
	return fread(buf, sizeof(XML_Char), bufLen, mFile);
}

//
// FileInputStream
//
FileOutputStream::FileOutputStream(const char *path)
{
	mFile = fopen(path, "w");
	if (!mFile)
		throw FileException(errno);
}

#ifdef _WIN32
FileOutputStream::FileOutputStream(const wchar_t *path)
{
	mFile = _wfopen(path, L"w");
	if (!mFile)
		throw FileException(errno);
}
#endif

FileOutputStream::~FileOutputStream()
{
	fclose(mFile);
}

int FileOutputStream::write(const char *buf, size_t bufLen)
{
	assert(buf);
	return fwrite(buf, sizeof(char), bufLen, mFile);
}

//
// FileException
//

const char *FileException::GetErrorString() const
{
	return strerror(mErrCode);
}

XML_END_NAMESPACE

