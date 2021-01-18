/*
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

#include "xmlio/xmlinputp.h"
#include <assert.h>
#include <string.h>

enum XmlToken_
{
	XML_TOK_EOF,
	XML_TOK_INVALID,		/* invalid token */
	XML_TOK_START_TAG,		/* <Tag> */
	XML_TOK_END_TAG,		/* </Tag> */
	XML_TOK_EMPTY_TAG_END,	/* /> */
	XML_TOK_ATTRIBUTE,		/* name="value" or name='value' */
	XML_TOK_TAG_END,		/* > */
	XML_TOK_DATA,			/* element data */
	XML_TOK_LINE,			/* a linefeed (LF, CR, CRLF) */
	XML_TOK_WHITESPACE,		/* some whitespace */
	XML_TOK_PI_START,		/* <?Tag attr="value"?> */
	XML_TOK_PI_END,			/* ?> */
	XML_TOK_COMMENT_START,	/* <!-- */
	XML_TOK_COMMENT_END,	/* --> */
	XML_TOK_CDATA_START,	/* <![CDATA[ */
	XML_TOK_CDATA_END,		/* ]]> */
	XML_TOK_DOCTYPE			/* <!DOCTYPE xxx [ */
};

typedef int XmlToken;

/* return !0 if the strings match, 0 if not */
int XML_StringsMatch(const XML_Char *s1, const XML_Char *s2)
{
	assert(s1);
	assert(s2);

	while (*s1 != 0 && *s1 == *s2)
	{
		s1++;
		s2++;
	}
	return *s1 == 0 && *s2 == 0;
}

/* return !0 if the character is whitespace */
int XML_IsWhiteSpace(XML_Char c)
{
	return (c == ' ') || (c == '\t') || (c == '\r') || (c == '\n');
}

static XML_Char peekChar(XML_Input *input)
{
	if (input->bufPtr >= input->buffer + input->bufSize)
	{
		XML_InputStream *stream = input->stream;
		input->bufSize = (*stream->readProc)(stream, input->buffer, input->maxBufSize);
		if (input->bufSize <= 0)
			return XML_EOF;
		input->bufPtr = input->buffer;
	}
	return *input->bufPtr;
}

static XML_Char nextChar(XML_Input *input)
{
 	if (input->bufPtr >= input->buffer + input->bufSize)
	{
		XML_InputStream *stream = input->stream;
		input->bufSize = (*stream->readProc)(stream, input->buffer, input->maxBufSize);
		if (input->bufSize <= 0)
			return XML_EOF;
		input->bufPtr = input->buffer;
	}
	input->offset++;
	input->bufPtr++;
	return input->bufPtr[-1];
}

static void ungetChar(XML_Input *input, XML_Char c)
{
	if (c != XML_EOF)
	{
		assert(input->bufPtr > input->buffer);
		input->bufPtr--;
		assert(*input->bufPtr == c);
		input->offset--;
	}
}

static XmlToken getNextToken(XML_Input *input, XML_Char *token, size_t *tokenLen, int expectAttrs)
{
	XML_Char *tptr = token;
	/*	int attr = 0;   mwall jan01 */

	/* look at first character to see if we're in whitespace */
	XML_Char c = nextChar(input);
	if (c == XML_EOF)
		return XML_TOK_EOF;

	/* check for a run of whitespace */
	if (c == ' ' || c == '\t')
	{
		do
		{
			*tptr++ = c;
			c = nextChar(input);
		} while (c == ' ' || c == '\t');
		ungetChar(input, c);
		*tptr = 0;
		input->column += (tptr - token);
		return XML_TOK_WHITESPACE;
	}
	/* check for a run of linefeeds */
	else if (c == '\r' || c == '\n')
	{
		*tptr++ = c;
		if (c == '\r')
		{
			// check for CRLF
			if (peekChar(input) == '\n')
				*tptr++ = nextChar(input);
		}
		else if (c == '\n')
		{
			// check for LFCR
			if (peekChar(input) == '\r')
				*tptr++ = nextChar(input);
		}
		*tptr = 0;
		input->line++;
		input->column = 0;
		return XML_TOK_LINE;
	}

	*tptr++ = c;
	if (expectAttrs)
	{
		/* we should hit an attr, ?>, or /> */
		if (c == '/')		/* check for '/>' */
		{
			c = nextChar(input);
			if (c == '>')
			{
				*tptr++ = c;
				*tptr = 0;
				input->column += 2;
				return XML_TOK_EMPTY_TAG_END;
			}
			return c == XML_EOF ? XML_TOK_EOF : XML_TOK_INVALID;
		}
		else if (c == '?')	/* check for '?>' */
		{
			c = nextChar(input);
			if (c == '>')
			{
				*tptr++ = c;
				*tptr = 0;
				input->column += 2;
				return XML_TOK_PI_END;
			}
			return c == XML_EOF ? XML_TOK_EOF : XML_TOK_INVALID;
		}
		else if (c == '>')	/* end of attribute list */
		{
			*tptr = 0;
			input->column++;
			return XML_TOK_TAG_END;
		}
		else
		{
			XML_Char quote;

			/* check for 'name="value"' */
			while (c != '=')
			{
				c = nextChar(input);
				if (c == XML_EOF)
					return XML_TOK_EOF;
				*tptr++ = c;
			}
			// get first quote
			c = nextChar(input);
			if (c != '\"' && c != '\'')
				return XML_TOK_INVALID;
			quote = c;
			*tptr++ = c;
			c = 0;
			// read until next quote is grabbed
			while (c != quote)
			{
				c = nextChar(input);
				if (c == XML_EOF)
					return XML_TOK_EOF;
				*tptr++ = c;
			}
			assert((size_t)(tptr - token) < *tokenLen);
			*tptr = 0;
			input->column += tptr - token;
			return XML_TOK_ATTRIBUTE;
		}
	}
	else if (c == '<')
	{
		c = peekChar(input);

		if (c == '/')
		{
			/* looks like an end tag */
			while (1)
			{
				c = nextChar(input);
				if (c == XML_EOF)
					return XML_TOK_EOF;
				else if (XML_IsWhiteSpace(c))
					return XML_TOK_INVALID;
				*tptr++ = c;
				if (c == '>')
					break;
			}
			*tptr = 0;
			input->column += (tptr - token);
			return XML_TOK_END_TAG;
		}
		else if (c == '?')
		{
			/* looks like a processing instruction */
			while (1)
			{
				c = nextChar(input);
				if (c == XML_EOF)
					return XML_TOK_EOF;
				if (XML_IsWhiteSpace(c))
				{
					ungetChar(input, c);
					break;
				}
				*tptr++ = c;
			}
			*tptr = 0;
			input->column += (tptr - token);
			return XML_TOK_PI_START;
		}
		else if (c == '!')
		{
			/* a comment or CDATA */
			*tptr++ = nextChar(input);
			c = peekChar(input);

			if (c == '-')		/* a comment */
			{
				*tptr++ = nextChar(input);
				*tptr++ = nextChar(input);		/* hopefully a '-' */

				*tptr = 0;
				input->column += (tptr - token);
				return XML_TOK_COMMENT_START;
			}
			else if (c == '[')	/* CDATA? */
			{
				*tptr++ = nextChar(input);
				*tptr++ = nextChar(input);		/* hopefully a 'C' */
				*tptr++ = nextChar(input);		/* hopefully a 'D' */
				*tptr++ = nextChar(input);		/* hopefully a 'A' */
				*tptr++ = nextChar(input);		/* hopefully a 'T' */
				*tptr++ = nextChar(input);		/* hopefully a 'A' */
				*tptr++ = nextChar(input);		/* hopefully a '[' */

				*tptr = 0;
				input->column += (tptr - token);
				return XML_TOK_CDATA_START;
			}
			else if (c == 'D')	/* DOCTYPE? */
			{
				*tptr++ = nextChar(input);		/* D */
				*tptr++ = nextChar(input);		/* hopefully a 'O' */
				*tptr++ = nextChar(input);		/* hopefully a 'C' */
				*tptr++ = nextChar(input);		/* hopefully a 'T' */
				*tptr++ = nextChar(input);		/* hopefully a 'Y' */
				*tptr++ = nextChar(input);		/* hopefully a 'P' */
				*tptr++ = nextChar(input);		/* hopefully a 'E' */

				*tptr = 0;
				input->column += (tptr - token);
				return XML_TOK_DOCTYPE;
			}
		}
		else
		{
			/* probably a tag start */
			while (1)
			{
				c = nextChar(input);
				if (c == XML_EOF)
					return XML_TOK_EOF;
				if (XML_IsWhiteSpace(c) || c == '/' || c == '>')
				{
					ungetChar(input, c);
					break;
				}
				*tptr++ = c;
			}
			*tptr = 0;
			input->column += (tptr - token);
			return XML_TOK_START_TAG;
		}
	}
	else
	{
		/* perhaps some data? read what we can of it */
		while ((size_t)(tptr - token) < *tokenLen)
		{
			c = nextChar(input);
			if (c == XML_EOF)
				return XML_TOK_EOF;
			if (c == '<')	/* new element - end of data */
			{
				ungetChar(input, c);
				break;
			}
			*tptr++ = c;
		}
		input->column += (tptr - token);
		*tokenLen = tptr - token;
		return XML_TOK_DATA;
	}
	return XML_TOK_INVALID;
}

static XML_Error skipComment(XML_Input *input)
{
	XML_Error error = XML_Error_None;
	int done = 0;
	while (!done)
	{
		XML_Char c = nextChar(input);
		if (c == XML_EOF)
			return XML_Error_EOF;

		if (c == '-')
		{
			if (peekChar(input) == '-')
			{
				nextChar(input);
				if (peekChar(input) == '>')
				{
					nextChar(input);
					done = 1;
				}
				else
					ungetChar(input, c);
			}
		}
	}
	return error;
}

static XML_Error skipCDATA(XML_Input *input)
{
	XML_Error error = XML_Error_None;
	int done = 0;
	while (!done)
	{
		XML_Char c = nextChar(input);
		if (c == XML_EOF)
			return XML_Error_EOF;

		if (c == ']')
		{
			if (peekChar(input) == ']')
			{
				nextChar(input);
				if (peekChar(input) == '>')
				{
					nextChar(input);
					done = 1;
				}
				else
					ungetChar(input, c);
			}
		}
	}
	return error;
}

static XML_Error skipDOCTYPE(XML_Input *input)
{
	XML_Error error = XML_Error_None;
	int done = 0;
	while (!done)
	{
		XML_Char c = nextChar(input);
		if (c == XML_EOF)
			return XML_Error_EOF;

		if (c == ']')
		{
			if (peekChar(input) == '>')
			{
				nextChar(input);
				done = 1;
			}
			else
				ungetChar(input, c);
		}
	}
	return error;
}

static XML_Error skip(XML_Input *input, int level, int expectData)
{
	int inData = expectData;
	XML_Error error = XML_Error_None;
	while (input->level > level)
	{
		if (inData)
		{
			while (inData)
			{
				XML_Char c = nextChar(input);
				if (c == XML_EOF)
					return XML_Error_EOF;
				else if (c == '\r' || c == '\n')
				{
					if (c == '\r')
					{
						/* check for CRLF */
						if (peekChar(input) == '\n')
							nextChar(input);
					}
					input->line++;
					input->column = 0;
				}
				else if (c == '<')
				{
					ungetChar(input, c);
					inData = 0;
				}
				else
				{
					input->column++;
				}
			}
		}
		else
		{
			XML_Char token[XML_TOKEN_MAX];
			size_t tokenLen = XML_TOKEN_MAX;
			XmlToken tok = getNextToken(input, token, &tokenLen, 0);

			if (tok == XML_TOK_EOF)
			{
				return XML_Error_EOF;
			}
			else if (tok == XML_TOK_INVALID)
				return XML_Error_Malformed;

			if (tok == XML_TOK_START_TAG)
			{
				int attrs = 1;
				input->level++;
				/* get attributes */
				while (attrs)
				{
					tok = getNextToken(input, token, &tokenLen, 1);
					if (tok == XML_EOF)
						return XML_Error_EOF;
					else if (tok == XML_TOK_EMPTY_TAG_END)
					{
						attrs = 0;
						input->level--;
					}
					else if (tok == XML_TOK_TAG_END)
					{
						attrs = 0;
						inData = 1;
					}
				}
			}
			else if (tok == XML_TOK_END_TAG)
			{
				input->level--;
			}
			else if (tok == XML_TOK_COMMENT_START)
			{
				skipComment(input);
			}
			else if (tok == XML_TOK_CDATA_START)
			{
				skipCDATA(input);
			}
		}
	}

	return error;
}

const XML_Char *XML_InputGetErrorString(XML_Error error)
{
	static const XML_Char *errorStrings[XML_Errors] = 
	{
		"no error",
		"not well-formed",
		"invalid value",
		"premature eof",
	};
	return errorStrings[error];
}

XML_Error XML_InputGetError(const XML_Input *input)
{
	assert(input);
	return input->error;
}

static XML_Error elementHandler(XML_Input *input, XML_Element *elem, const XML_Handler *handler, void *userData)
{
	if (handler->size > 0 && handler->info.Element.userData == NULL)
	{
		/* compute user-data from offset */
		userData = (void *)((char *)userData + handler->offset);
	}
	else if (handler->info.Element.userData)
		userData = handler->info.Element.userData;

	return (*handler->info.Element.proc)(elem, userData);
}

static XML_Error dataHandler(XML_Input *input, const XML_Char *data, size_t len, const XML_Handler *handler, void *userData)
{
	return (*handler->info.Data.proc)(data, len, userData);
}

static XML_Error intHandler(XML_Input *input, XML_Element *elem, const XML_Handler *handler, void *userData)
{
	XML_Char tmp[40];
	size_t len = sizeof(tmp) / sizeof(XML_Char);
	XML_Error error = XML_ElementReadData(elem, tmp, &len);
	tmp[len] = 0;
	if (error == XML_Error_None)
	{
		int value = (int)atol(tmp);
		int *result = handler->info.Int.result ? handler->info.Int.result : (int *)((char *)userData + handler->offset);
		if (handler->info.Int.maxVal != 0 || handler->info.Int.minVal != 0)
		{
			/* do range checking */
			if (value < handler->info.Int.minVal)
				value = handler->info.Int.minVal;
			else if (value > handler->info.Int.maxVal)
				value = handler->info.Int.maxVal;
		}
		*result = value;
	}
	return error;
}

static XML_Error uintHandler(XML_Input *input, XML_Element *elem, const XML_Handler *handler, void *userData)
{
	XML_Char tmp[40];
	size_t len = sizeof(tmp) / sizeof(XML_Char);
	XML_Error error = XML_ElementReadData(elem, tmp, &len);
	tmp[len] = 0;
	if (error == XML_Error_None)
	{
		unsigned int value = (unsigned int)atol(tmp);
		unsigned int *result = handler->info.UInt.result ? handler->info.UInt.result : (unsigned int *)((char *)userData + handler->offset);
		if (handler->info.UInt.maxVal != 0 || handler->info.UInt.minVal != 0)
		{
			/* do range checking */
			if (value < handler->info.UInt.minVal)
				value = handler->info.UInt.minVal;
			else if (value > handler->info.UInt.maxVal)
				value = handler->info.UInt.maxVal;
		}
		*result = value;
	}
	return error;
}

static XML_Error floatHandler(XML_Input *input, XML_Element *elem, const XML_Handler *handler, void *userData)
{
	XML_Char tmp[40];
	size_t len = sizeof(tmp) / sizeof(XML_Char);
	XML_Error error = XML_ElementReadData(elem, tmp, &len);
	tmp[len] = 0;
	if (error == XML_Error_None)
	{
		float value = (float)atof(tmp);
		float *result = handler->info.Float.result ? handler->info.Float.result : (float *)((char *)userData + handler->offset);
		if (handler->info.Float.maxVal != 0 || handler->info.Float.minVal != 0)
		{
			/* do range checking */
			if (value < handler->info.Float.minVal)
				value = handler->info.Float.minVal;
			else if (value > handler->info.Float.maxVal)
				value = handler->info.Float.maxVal;
		}
		*result = value;
	}
	return error;
}

static XML_Error doubleHandler(XML_Input *input, XML_Element *elem, const XML_Handler *handler, void *userData)
{
	XML_Char tmp[80];
	size_t len = sizeof(tmp) / sizeof(XML_Char);
	XML_Error error = XML_ElementReadData(elem, tmp, &len);
	tmp[len] = 0;
	if (error == XML_Error_None)
	{
		double value = atof(tmp);
		double *result = handler->info.Double.result ? handler->info.Double.result : (double *)((char *)userData + handler->offset);
		/* minVal and maxVal are POINTERS to doubles */
		if (handler->info.Double.minVal != NULL)
		{
			if (value < *(handler->info.Double.minVal))
				value = (*handler->info.Double.minVal);
		}
		if (handler->info.Double.maxVal != NULL)
		{
			if (value > *(handler->info.Double.maxVal))
				value = (*handler->info.Double.maxVal);
		}
		*result = value;
	}
	return error;
}

static XML_Error stringHandler(XML_Input *input, XML_Element *elem, const XML_Handler *handler, void *userData)
{
	size_t len = handler->info.CString.maxLen / sizeof(XML_Char);
	XML_Char *str = handler->info.CString.result ? handler->info.CString.result : (XML_Char *)((char *)userData + handler->offset);
	XML_Error error = XML_ElementReadData(elem, str, &len);
	str[len] = 0;
	return error;
}

static void setValue(void *result, long value, size_t size)
{
	switch (size)
	{
		case sizeof(unsigned char):
			*((unsigned char *)result) = (unsigned char)value;
			break;

		case sizeof(unsigned short):
			*((unsigned short *)result) = (unsigned short)value;
			break;

		case sizeof(unsigned int):
			*((unsigned int *)result) = (unsigned int)value;
			break;

		default:
			/* unknown size */
			break;
	}
}

static XML_Error boolHandler(XML_Input *input, XML_Element *elem, const XML_Handler *handler, void *userData)
{
	XML_Char tmp[40];
	size_t len = sizeof(tmp) / sizeof(XML_Char);
	XML_Error error = XML_ElementReadData(elem, tmp, &len);
	tmp[len] = 0;	/* must NULL-terminate */
	if (error == XML_Error_None)
	{
		int value = XML_StringsMatch(tmp, "True") != 0 || XML_StringsMatch(tmp, "true") != 0;
		void *result = handler->info.Bool.result ? handler->info.Bool.result : (void *)((char *)userData + handler->offset);
		setValue(result, value, handler->size);
	}
	return error;
}

static XML_Error listHandler(XML_Input *input, XML_Element *elem, const XML_Handler *handler, void *userData)
{
	XML_Char tmp[80];
	size_t len = sizeof(tmp) / sizeof(XML_Char);
	XML_Error error = XML_ElementReadData(elem, tmp, &len);
	tmp[len] = 0;	/* must NULL-terminate */
	if (error == XML_Error_None)
	{
		/* loop over list entries and compare */
		void *result = handler->info.List.result ? handler->info.List.result : (void *)((char *)userData + handler->offset);
		int i;
		for (i = 0; i < handler->info.List.listSize; i++)
		{
			if (XML_StringsMatch(tmp, handler->info.List.list[i]))
			{
				setValue(result, i, handler->size);
				break;
			}
		}
	}
	return error;
}

static XML_Error handleHandler(XML_Input *input, XML_Element *elem, const XML_Handler *handler, void *userData)
{
	XML_Error error = XML_Error_None;
	int offset = input->offset;
	switch (handler->type)
	{
		case XML_Handler_Element:
			error = input->elementHandler(input, elem, handler, userData);
			break;

		case XML_Handler_Int:
			error = intHandler(input, elem, handler, userData);
			break;

		case XML_Handler_UInt:
			error = uintHandler(input, elem, handler, userData);
			break;

		case XML_Handler_Float:
			error = floatHandler(input, elem, handler, userData);
			break;

		case XML_Handler_Double:
			error = doubleHandler(input, elem, handler, userData);
			break;

		case XML_Handler_CString:
			error = stringHandler(input, elem, handler, userData);
			break;

		case XML_Handler_Bool:
			error = boolHandler(input, elem, handler, userData);
			break;

		case XML_Handler_List:
			error = listHandler(input, elem, handler, userData);
			break;

			/* mwall-begin jan01 */
	        case XML_Handler_None:
		case XML_Handler_Chain:
	        case XML_Handler_Data:
	        case XML_Handler_CDATA:
	        case XML_Handler_Comment:
	        case XML_Handler_Types:
		  error = XML_Error_None;
		  break;
			/* mwall-end jan01 */

	}
	if (error == XML_Error_None)
	{
		if (input->offset == offset && !elem->empty)
		{
			/* handler did nothing, so skip this element */
			error = skip(input, input->level - 1, !elem->empty);
		}
	}
	return error;
}

XML_Input *XML_InputCreate(XML_InputStream *stream)
{
	XML_Input *input;
	
	assert(stream);
	
	input = (XML_Input *)malloc(sizeof(XML_Input));
	if (input)
	{
		/* set up function pointers */
		input->elementHandler = elementHandler;
		input->dataHandler = dataHandler;

		input->stream = stream;
		input->level = 0;
		input->offset = 0;
		input->line = 0;
		input->column = 0;
		input->userData = NULL;
		input->error = XML_Error_None;

		/* allocate our working buffer */
		input->maxBufSize = XML_BUFFER_SIZE;
		input->bufSize = 0;
		input->buffer = (XML_Char *)malloc(input->maxBufSize);
		if (!input->buffer)
		{
			free(input);
			return NULL;
		}
		input->bufPtr = input->buffer;

		/* allocate a block of attributes */
		input->attrPool = (XML_Attribute *)malloc(sizeof(XML_Attribute) * XML_ATTR_MAX);
		if (!input->attrPool)
		{
			free(input->buffer);
			free(input);
			return NULL;
		}
		/* and link them together into a free list */
		{
			int i;
			XML_Attribute *attr = input->attrPool;
			for (i = 0; i < XML_ATTR_MAX - 1; i++, attr++)
				attr->next = attr + 1;
			attr->next = NULL;
			input->nextAttr = input->attrPool;
		}
		input->attrsUsed = 0;

		return input;
	}
	return input;
}

void XML_InputFree(XML_Input *input)
{
	if (input)
	{
		if (input->attrPool) free(input->attrPool);
		if (input->buffer) free(input->buffer);
		free(input);
	}
}

/* return an unused attribute from the free pool */
static XML_Attribute *newAttribute(XML_Input *input)
{
	XML_Attribute *attr = input->nextAttr;
	assert(input->attrsUsed < XML_ATTR_MAX);
	assert(attr);

	input->nextAttr = attr->next;
	memset(attr, 0, sizeof(XML_Attribute));
	input->attrsUsed++;
	return attr;
}

/* return an attribute back to the free pool */
static void freeAttribute(XML_Input *input, XML_Attribute *attr)
{
	assert(input);
	assert(attr);
	assert(attr->next == NULL);

	assert(input->attrsUsed > 0);

	attr->next = input->nextAttr;
	input->nextAttr = attr;
	input->attrsUsed--;
}

/* associate user-data with the Input */
void XML_InputSetUserData(XML_Input *input, void *userData)
{
	assert(input);
	input->userData = userData;
}

/* return the user-data associated with the Input */
void *XML_InputGetUserData(const XML_Input *input)
{
	assert(input);
	return input->userData;
}

/* return the current line number */
int XML_InputGetLine(const XML_Input *input)
{
	assert(input);
	return input->line + 1;
}

/* return the current column number */
int XML_InputGetColumn(const XML_Input *input)
{
	assert(input);
	return input->column + 1;
}

/* return the current character offset */
int XML_InputGetOffset(const XML_Input *input)
{
	assert(input);
	return input->offset;
}

/* find a handler for this element name */
static const XML_Handler *findHandler(const XML_Char *name, const XML_Handler handlers[], void **userData)
{
	int i = 0;
	while (handlers[i].type != XML_Handler_None)
	{
		/* deal with chains */
		if (handlers[i].type == XML_Handler_Chain)
		{
			const XML_Handler *handler = findHandler(name, handlers[i].info.Chain.handlers, userData);
			if (handler)
			{
				/* return user-data specific to the handler chain */
				if (handlers[i].size > 0)
					*userData = (void *)((char *)(*userData) + handlers[i].offset);
				else if (handlers[i].info.Chain.userData)
					*userData = handlers[i].info.Chain.userData;
				return handler;
			}
		}

		if (handlers[i].name && XML_StringsMatch(name, handlers[i].name))
			return &handlers[i];
		i++;
	}
	/* no exact match found - look for a generic (unnamed) handler */
	i = 0;
	while (handlers[i].type != XML_Handler_None)
	{
		if (!handlers[i].name && handlers[i].type == XML_Handler_Element)
			return &handlers[i];
		i++;
	}
	return NULL;
}

/* search for a data handler */
static const XML_Handler *findDataHandler(const XML_Handler handlers[], void **userData, XML_HandlerType type)
{
	int i = 0;
	while (handlers[i].type != XML_Handler_None)
	{
		/* deal with chains */
		if (handlers[i].type == XML_Handler_Chain)
		{
			const XML_Handler *handler = findDataHandler(handlers[i].info.Chain.handlers, userData, type);
			if (handler)
			{
				/* return user-data specific to the handler chain */
				if (handlers[i].size > 0)
					*userData = (void *)((char *)(*userData) + handlers[i].offset);
				else if (handlers[i].info.Chain.userData)
					*userData = handlers[i].info.Chain.userData;
				return handler;
			}
		}
		if (handlers[i].type == type)
			return &handlers[i];
		i++;
	}
	return NULL;
}


static XML_Error handlePI(XML_Input *input, XML_Char *pi, const XML_Handler handlers[], void *userData)
{
	XML_Error error = XML_Error_None;
	int gettingAttrs = 1;
	while (gettingAttrs && error == XML_Error_None)
	{
		XML_Char token[XML_TOKEN_MAX];
		size_t tokenLen = XML_TOKEN_MAX;
		XmlToken tok = getNextToken(input, token, &tokenLen, 1);
		switch (tok)
		{
			case XML_TOK_ATTRIBUTE:
				break;

			case XML_TOK_WHITESPACE:
				break;

			case XML_TOK_PI_END:
				gettingAttrs = 0;
				break;

			default:
				gettingAttrs = 0;
				error = XML_Error_Malformed;
				break;
		}
	}
	return error;
}

/* break the token of the form [name="value"] into separate components */
static void getAttribute(const XML_Char *token, XML_Char *name, XML_Char *value)
{
#if defined(DEBUG)
	const XML_Char *nameStart = name;
	const XML_Char *valueStart = value;
#endif
	XML_Char quote;

	/* copy up to '=' into name */
	while (*token && *token != '=')
		*name++ = *token++;
	*name = 0;
#if defined(DEBUG)
	assert(name - nameStart < XML_ATTR_NAME_MAX);
#endif
	assert(*token == '=');
	token++;
	quote = *token++;
	assert(quote == '\"' || quote == '\'');
	while (*token && *token != quote)
		*value++ = *token++;
	*value = 0;
#if defined(DEBUG)
	assert(value - valueStart < XML_ATTR_VALUE_MAX);
#endif
}

static XML_Error handleData(XML_Input *input, const XML_Char *token, size_t len, const XML_Handler *handler, void *userData)
{
	XML_Char data[XML_TOKEN_MAX];
	const XML_Char *end = data + len;
	int done = 0;
	XML_Error error = XML_Error_None;

	XML_Char *ptr = data;

	/* pre-initialize data with the token */
	strncpy(data, token, len);
	ptr += len;

	assert(handler->type == XML_Handler_Data);
	while (!done && error == XML_Error_None)
	{
		while (ptr < end && !done)
		{
			XML_Char c = nextChar(input);
			if (c == XML_EOF)
				return XML_Error_EOF;

			if (c == '<')	/* new element */
			{
				ungetChar(input, c);
				done = 1;
			}
			else
				*ptr++ = c;
		}
		if (!done && handler)
			error = input->dataHandler(input, data, ptr - data, handler, userData);

		/* reset pointer */
		ptr = data;
	}
	/* signal end-of-data*/
	if (handler && error == XML_Error_None)
		error = input->dataHandler(input, data, 0, handler, userData);
	return error;
}

static XML_Error handleComment(XML_Input *input, const XML_Handler *handler, void *userData)
{
	if (handler)
	{
		XML_Char data[XML_TOKEN_MAX];
		int done = 0;
		XML_Error error = XML_Error_None;

		assert(handler->type == XML_Handler_Comment);
		/* read comment data */
		while (!done && error == XML_Error_None)
		{
			XML_Char *ptr = data;
			const XML_Char *end = ptr + XML_TOKEN_MAX;

			while (ptr < end && !done)
			{
				XML_Char c = nextChar(input);
				if (c == XML_EOF)
					return XML_Error_EOF;

				if (c == '-')
				{
					if (peekChar(input) == '-')
					{
						nextChar(input);
						if (peekChar(input) == '>')
						{
							nextChar(input);
							done = 1;
						}
						else
							ungetChar(input, c);
					}
				}
				if (!done)
					*ptr++ = c;
			}
			error = input->dataHandler(input, data, ptr - data, handler, userData);
		}
		/* signal end-of-comment */
		if (error == XML_Error_None)
			error = input->dataHandler(input, data, 0, handler, userData);
		return error;
	}
	else
		return skipComment(input);
}

static XML_Error handleCDATA(XML_Input *input, const XML_Handler *handler, void *userData)
{
	if (handler)
	{
		XML_Char data[XML_TOKEN_MAX];
		int done = 0;
		XML_Error error = XML_Error_None;

		assert(handler->type == XML_Handler_CDATA);
		while (!done && error == XML_Error_None)
		{
			XML_Char *ptr = data;
			const XML_Char *end = ptr + XML_TOKEN_MAX;

			while (ptr < end && !done)
			{
				XML_Char c = nextChar(input);
				if (c == XML_EOF)
					return XML_Error_EOF;

				if (c == ']')
				{
					if (peekChar(input) == ']')
					{
						nextChar(input);
						if (peekChar(input) == '>')
						{
							nextChar(input);
							done = 1;
						}
						else
							ungetChar(input, c);
					}
				}
				if (!done)
					*ptr++ = c;
			}
			error = input->dataHandler(input, data, ptr - data, handler, userData);
		}
		/* signal end-of-CDATA */
		if (error == XML_Error_None)
			error = input->dataHandler(input, data, 0, handler, userData);
		return error;
	}
	else
		return skipCDATA(input);
}

static XML_Error handleDOCTYPE(XML_Input *input, const XML_Handler *handler, void *userData)
{
	return skipDOCTYPE(input);
}

static XML_Error handleElement(XML_Input *input, XML_Char *elemName, const XML_Handler handlers[], void *userData)
{
	XML_Error error = XML_Error_None;
	int gettingAttrs = 1;
	int emptyTag = 0;
	XML_Attribute *attrs = NULL;
	XML_Attribute *lastAttr = NULL;

	while (gettingAttrs && error == XML_Error_None)
	{
		XML_Char token[XML_TOKEN_MAX];
		size_t tokenLen = XML_TOKEN_MAX;
		XmlToken tok = getNextToken(input, token, &tokenLen, 1);
		switch (tok)
		{
			case XML_TOK_ATTRIBUTE:
			{
				/* add a new attribute to the list */
				XML_Attribute *attr = newAttribute(input);
				assert(attr);
				getAttribute(token, attr->name, attr->value);
				if (attrs)
				{
					assert(lastAttr);
					lastAttr->next = attr;
				}
				else
					attrs = attr;
				lastAttr = attr;
				break;
			}

			case XML_TOK_WHITESPACE:
				break;

			case XML_TOK_EMPTY_TAG_END:	/* /> */
				emptyTag = 1;
				/* fall-through */
			case XML_TOK_TAG_END:		/* > */
				gettingAttrs = 0;
				break;

			default:
				gettingAttrs = 0;
				error = XML_Error_Malformed;
				break;
		}
	}

	if (error == XML_Error_None)
	{
		const XML_Handler *handler = findHandler(elemName, handlers, &userData);
		if (handler)
		{
			XML_Element elem;
			input->level++;

			elem.input = input;
			strcpy(elem.name, elemName);
			elem.attrs = attrs;
			elem.empty = emptyTag;
			elem.level = input->level;
			error = handleHandler(input, &elem, handler, userData);
			if (emptyTag)
				input->level--;
		}
		else if (!emptyTag)
		{
			/* skip this entire block */
			int currLevel = input->level;
			input->level++;
			error = skip(input, currLevel, 1);
		}
	}
	/* free the attributes */
	while (attrs)
	{
		XML_Attribute *next = attrs->next;
		attrs->next = NULL;
		freeAttribute(input, attrs);
		attrs = next;
	}

	return error;
}

XML_Error XML_InputProcess(XML_Input *input, const XML_Handler handlers[], void *userData)
{
	int done = 0;
	int level = input->level;
	const XML_Handler *dataHandler = NULL;
	void *dataUserData = userData;

	assert(input);
	assert(handlers);

	/* search for a data handler in this handler set */
	dataHandler = findDataHandler(handlers, &dataUserData, XML_Handler_Data);

	while (!done && input->error == XML_Error_None)
	{
		XML_Char token[XML_TOKEN_MAX];
		size_t tokenLen = XML_TOKEN_MAX;
		XmlToken tok = getNextToken(input, token, &tokenLen, 0);
		switch (tok)
		{
			case XML_TOK_EOF:
				done = 1;
				break;

			case XML_TOK_INVALID:
				done = 1;
				input->error = XML_Error_Malformed;
				break;

			case XML_TOK_WHITESPACE:
				if (dataHandler)
					input->error = input->dataHandler(input, token, strlen(token), dataHandler, dataUserData);
				break;

			case XML_TOK_LINE:
				if (dataHandler)
					input->error = input->dataHandler(input, token, strlen(token), dataHandler, dataUserData);
				break;

			case XML_TOK_PI_START:
				input->error = handlePI(input, token + 2, handlers, userData);
				break;

			case XML_TOK_START_TAG:
			{
				input->error = handleElement(input, token + 1, handlers, userData);
				break;
			}

			case XML_TOK_END_TAG:
				input->level--;
				if (input->level == level - 1)
					done = 1;
				break;

			case XML_TOK_DATA:
				input->error = handleData(input, token, tokenLen, dataHandler, dataUserData);
				break;

			case XML_TOK_COMMENT_START:
			{
				void *data = userData;
				const XML_Handler *handler = findDataHandler(handlers, &data, XML_Handler_Comment);
				input->error = handleComment(input, handler, data);
				break;
			}

			case XML_TOK_CDATA_START:
			{
				void *data = userData;
				const XML_Handler *handler = findDataHandler(handlers, &data, XML_Handler_CDATA);
				input->error = handleCDATA(input, handler, data);
				break;
			}

			case XML_TOK_DOCTYPE:
			{
				input->error = handleDOCTYPE(input, NULL, NULL);
				break;
			}

			default:
				done = 1;
				input->error = XML_Error_Malformed;
				break;
		}
	}
	return input->error;
}

/* process the sub-element with a new set of handlers */
XML_Error XML_ElementProcess(XML_Element *elem, const XML_Handler handlers[], void *userData)
{
	assert(elem);

	if (elem->empty)
		return XML_Error_None;

	return XML_InputProcess(elem->input, handlers, userData);
}

/* read data into the buffer until we see a '<' */
XML_Error XML_ElementReadData(XML_Element *elem, XML_Char *data, size_t *readSize)
{
	XML_Char *ptr = data;
	XML_Input *input = elem->input;
	size_t maxSize = *readSize;

	assert(readSize);
	assert(elem);
	assert(data);

	while ((size_t)(ptr - data) <= maxSize)
	{
		XML_Char c = nextChar(input);
		if (c == XML_EOF)
		{
			input->error = XML_Error_EOF;
			return input->error;
		}
		else if (c == '\r')
		{
			/* check for CRLF */
			if (peekChar(input) == '\n')
				nextChar(input);

			input->line++;
			input->column = 0;
			*ptr++ = '\n';
		}
		else if (c == '<')
		{
			ungetChar(input, c);
			break;
		}
		else
		{
			*ptr++ = c;
			input->column++;
		}
	}
	/* return the actual number of characters read */
	*readSize = ptr - data;

	return XML_Error_None;
}

/* return the Element name */
const XML_Char *XML_ElementGetName(const XML_Element *elem)
{
	assert(elem);
	return elem->name;
}

/* return the Input pointer that generated this Element */
const XML_Input *XML_ElementGetInput(const XML_Element *elem)
{
	assert(elem);
	return elem->input;
}

/* return the Element level (0 is the master document level) */
int XML_ElementGetLevel(const XML_Element *elem)
{
	assert(elem);
	assert(elem->level > 0);
	return elem->level - 1;
}

XML_Error XML_ElementGetError(const XML_Element *elem)
{
	const XML_Input *input = XML_ElementGetInput(elem);
	assert(input);
	return XML_InputGetError(input);
}

int XML_ElementIsEmpty(const XML_Element *elem)
{
	assert(elem);
	return elem->empty;
}

/* return the number of attributes found for this Element */
int XML_ElementGetNumAttrs(const XML_Element *elem)
{
	int i = 0;
	const XML_Attribute *attr = elem->attrs;
	while (attr)
	{
		i++;
		attr = attr->next;
	}
	return i;
}

/* return the Attribute name of the specified attribute number */
const XML_Char *XML_ElementGetAttrName(const XML_Element *elem, int index)
{
	int i = 0;
	const XML_Attribute *attr = elem->attrs;
	while (attr)
	{
		if (i == index)
			return attr->name;
		attr = attr->next;
		i++;
	}
	return NULL;
}

/* return the Attribute value of the specified attribute number */
const XML_Char *XML_ElementGetAttrValue(const XML_Element *elem, int index)
{
	int i = 0;
	const XML_Attribute *attr = elem->attrs;
	while (attr)
	{
		if (i == index)
			return attr->value;
		attr = attr->next;
		i++;
	}
	return NULL;
}

/* look for an Attribute with the given name and return the value */
const XML_Attribute *XML_ElementFindAttr(const XML_Element *elem, const XML_Char *name)
{
  /*	int i = 0;  mwall jan01 */
	const XML_Attribute *attr = elem->attrs;
	while (attr)
	{
		if (XML_StringsMatch(attr->name, name))
			return attr;
		attr = attr->next;
	}
	return NULL;
}

/* return the first Element Attribute */
const XML_Attribute *XML_ElementGetAttrList(const XML_Element *elem)
{
	assert(elem);
	return elem->attrs;
}

/* return the next Attribute, or NULL if there are no more */
const XML_Attribute *XML_AttrGetNext(const XML_Attribute *attr)
{
	assert(attr);
	return attr->next;
}

/* return the Attribute name */
const XML_Char *XML_AttrGetName(const XML_Attribute *attr)
{
	assert(attr);
	return attr->name;
}

/* return the Attribute value */
const XML_Char *XML_AttrGetValue(const XML_Attribute *attr)
{
	return attr ? attr->value : NULL;
}

int XML_AttrGetInt(const XML_Attribute *attr, int defValue)
{
	return attr ? atoi(attr->value) : defValue;
}

unsigned int XML_AttrGetUInt(const XML_Attribute *attr, unsigned int defValue)
{
	return attr ? atoi(attr->value) : defValue;
}

float XML_AttrGetFloat(const XML_Attribute *attr, float defValue)
{
	return attr ? (float)atof(attr->value) : defValue;
}

double XML_AttrGetDouble(const XML_Attribute *attr, double defValue)
{
	return attr ? atof(attr->value) : defValue;
}

int XML_AttrGetBoolean(const XML_Attribute *attr, int defValue)
{
	if (attr)
	{
		const XML_Char *value = attr->value;
		if (XML_StringsMatch(value, "true") ||
			XML_StringsMatch(value, "True"))
			return 1;
		else if (XML_StringsMatch(value, "false") ||
			XML_StringsMatch(value, "False"))
			return 0;
	}
	return defValue;
}
