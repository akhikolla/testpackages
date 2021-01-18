/*
	Private implementation structures - DO NOT INCLUDE DIRECTLY!

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

#include "xmlstream.h"
#include "xmlinput.h"

#define XML_TOKEN_MAX		2048 	/* maximum length of a run */
#define XML_ELEM_NAME_MAX	80		/* maximum element name length */
#define XML_ATTR_MAX		50		/* maximum attributes in use at once */
#define XML_ATTR_NAME_MAX	32		/* maximum attribute name length */
#define XML_ATTR_VALUE_MAX	2048	/* maximum attribute value length */
#define XML_BUFFER_SIZE		2048	/* maximum number of characters to buffer at once */

#define XML_EOF				((XML_Char)-1)

typedef struct XML_Attribute_
{
	struct XML_Attribute_ *next;			/* next attribute in list */
	XML_Char name[XML_ATTR_NAME_MAX];
	XML_Char value[XML_ATTR_VALUE_MAX];
} XML_Attribute;

typedef struct XML_Element_
{
	struct XML_Input_ *input;
	XML_Char name[XML_ELEM_NAME_MAX];
	struct XML_Attribute_ *attrs;		/* attribute list */
	int level;					/* element level (0 is root) */
	int empty;					/* "empty" element */
} XML_Element;

typedef XML_Error (*ElementHandlerProc)(struct XML_Input_ *input, XML_Element *elem, const XML_Handler *handler, void *userData);
typedef XML_Error (*DataHandlerProc)(struct XML_Input_ *input, const XML_Char *data, size_t len, const XML_Handler *handler, void *userData);

typedef struct XML_Input_
{
	struct XML_InputStream_ *stream;	/* input source */
	int level;					/* current nest level */
	int column;					/* current column */
	int line;					/* current line number */
	int offset;					/* current character offset */
	void *userData;				/* user-supplied data */
	XML_Char *buffer;			/* our read buffer */
	size_t maxBufSize;			/* max buffer size */
	size_t bufSize;				/* # of bytes in buffer */
	XML_Char *bufPtr;			/* our position in the buffer */
	XML_Error error;			/* current error */
	XML_Attribute *attrPool;	/* pool of attributes */
	XML_Attribute *nextAttr;	/* next available attribute */
	int attrsUsed;				/* number of used attributes */
	ElementHandlerProc elementHandler;
 	DataHandlerProc dataHandler;
} XML_Input;

