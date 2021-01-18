/*
	Defines interfaces for input and output streams needed by the inputter and
	outputter. These should be subclassed by application-specific IO
	classes.

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

#ifndef XMLSTREAM_H
#define XMLSTREAM_H

#include "xmlconfig.h"

#ifdef __cplusplus
	extern "C" {
#endif

/*
	C Interface
*/

struct XML_InputStream_;
struct XML_OutputStream_;

/** Read function
	@param data		buffer to receive characters
	@param len		number of characters to read
	@return			number of characters read, 0 if end of file, -1 if there is an error
*/	
typedef int (*XML_ReadProc)(struct XML_InputStream_ *str, XML_Char *data, size_t len);

/** Write function
	@param data		buffer containing characters
	@param len		number of characters to write
	@return			number of characters written, -1 if there is an error
*/	
typedef int (*XML_WriteProc)(struct XML_OutputStream_ *str, const XML_Char *data, size_t len);

/**
	Input Stream interface required by XML_Input.
*/
typedef struct XML_InputStream_
{
	XML_ReadProc readProc;
} XML_InputStream;

/**
	Output Stream interface required by XML_Output.
*/
typedef struct XML_OutputStream_
{
	XML_WriteProc writeProc;
} XML_OutputStream;

#ifdef __cplusplus
}	// extern "C"

//
// C++ Interface
//
XML_BEGIN_NAMESPACE

/**
	Input Stream interface required by XML::Input.
*/
class InputStream
{
public:
	/** read up to bufLen characters from an input source and place in buf.
		@param buf		destination buffer
		@param bufLen	maximum size of destination buffer
		@return 		the number of characters actually read, 0 if eof
	*/
	virtual int read(XML_Char *buf, size_t bufLen) = 0;
};

/**
	Output Stream interface required by XML::Output.
*/
class OutputStream
{
public:
	/** write up to bufLen characters to an output source
		@param buf		source buffer
		@param bufLen	number of characters to write
		@return 		the number of characters actually written - if this
						number is less than bufLen, there was an error
	*/
	virtual int write(const char *buf, size_t bufLen) = 0;
};

XML_END_NAMESPACE
#endif	// __cplusplus

#endif	// XMLSTREAM_H

