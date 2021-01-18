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

#ifndef XMLINPUT_H
#define XMLINPUT_H

#include "xmlconfig.h"
#include <stddef.h>		/* offsetof */

#ifdef __cplusplus
	extern "C" {
#endif

/*
------------------------------------------------------------------------------
BEGIN C INTERFACE
------------------------------------------------------------------------------
*/

struct XML_Input_;
struct XML_Element_;
struct XML_Attribute_;
struct XML_InputStream_;		/* declared in xmlstream.h */

/**
	Input processing errors
*/
typedef enum XML_Error_
{
	XML_Error_None,
	XML_Error_Malformed,		/* XML data was malformed/incorrect */
	XML_Error_InvalidValue,		/* an invalid value was found for a built-in handler */
	XML_Error_EOF,				/* a premature End-Of-File was found */
	XML_Errors
} XML_Error;

/**
	Element handler types
*/
typedef enum XML_HandlerType_
{
	XML_Handler_None,
	XML_Handler_Element,
	XML_Handler_Int,
	XML_Handler_UInt,
	XML_Handler_Float,
	XML_Handler_Double,
	XML_Handler_Bool,
	XML_Handler_CString,
	XML_Handler_List,
	XML_Handler_Chain,
	XML_Handler_Data,
	XML_Handler_CDATA,
	XML_Handler_Comment,
	XML_Handler_Types
} XML_HandlerType;

/**
	Prototype for a custom element handler procedure
*/
typedef XML_Error (*XML_HandlerProc)(struct XML_Element_ *elem, void *userData);
typedef XML_Error (*XML_DataProc)(const XML_Char *data, size_t len, void *userData);

/* begin specifications for handler-specific extra data */
struct XML_ElementHandler_
{
	XML_HandlerProc proc;
	void *userData;
};

struct XML_DataHandler_
{
	XML_DataProc proc;
	void *userData;
};

struct XML_ChainHandler_
{
	const struct XML_Handler_ *handlers;
	void *userData;
};

struct XML_IntHandler_
{
	int *result;
	int minVal;
	int maxVal;
};

struct XML_UIntHandler_
{
	unsigned int *result;
	unsigned int minVal;
	unsigned int maxVal;
};

struct XML_FloatHandler_
{
	float *result;
	float minVal;
	float maxVal;
};

struct XML_DoubleHandler_
{
	double *result;
	double *minVal;
	double *maxVal;
};

struct XML_BoolHandler_
{
	int *result;
};

struct XML_CStringHandler_
{
	XML_Char *result;				/* address of string */
	size_t maxLen;					/* maximum size of string in characters (includes ending 0) */
};

struct XML_ListHandler_
{
	int *result;					/* list index */
	const XML_Char *const *list;	/* NULL-terminated list values */
	int listSize;					/* entries in list */
};
/* end specifications for handler-specific extra data */

/**
	XML_Handler - specifies how to deal with specific elements. Stored in
	a (usually static) array and passed to processing functions.
*/
typedef struct XML_Handler_
{
	const char *name;		/* name of element (or NULL for a generic handler) */
	XML_HandlerType type;	/* type of handler */
	size_t offset;			/* offset of destination address */
	size_t size;			/* size of destination object */
	union {
		struct { void *dummy[3]; } Dummy;
		struct XML_ElementHandler_ Element;
		struct XML_DataHandler_ Data;
		struct XML_ChainHandler_ Chain;
		struct XML_IntHandler_ Int;
		struct XML_UIntHandler_ UInt;
		struct XML_FloatHandler_ Float;
		struct XML_DoubleHandler_ Double;
		struct XML_BoolHandler_ Bool;
		struct XML_CStringHandler_ CString;
		struct XML_ListHandler_ List;
	} info;
} XML_Handler;

#define XML_MEMBER_OFFSET(object, member)	offsetof(object, member)
#define XML_MEMBER_SIZE(object, member)		sizeof(((object *)0)->member)
#define XML_OBJECT_MEMBER(object, member)	XML_MEMBER_OFFSET(object, member), XML_MEMBER_SIZE(object, member)

#define XML_DATA_HANDLER(name, proc, data) 	{ name, XML_Handler_Data, 0, 0, { (void *)proc, (void *)data, NULL } }

#define XML_CDATA_HANDLER(name, proc, data) { name, XML_Handler_CDATA, 0, 0, { (void *)proc, (void *)data, NULL } }

#define XML_COMMENT_HANDLER(name, proc, data) { name, XML_Handler_Comment, 0, 0, { (void *)proc, (void *)data, NULL } }

#define XML_ELEMENT_HANDLER(name, proc, data) { name, XML_Handler_Element, 0, 0, { (void *)proc, (void *)data, NULL } }

#define XML_ELEMENT_HANDLER_MEMBER(name, proc, object, member) { name, XML_Handler_Element, XML_OBJECT_MEMBER(object, member), { (void *)proc, NULL, NULL} }

#define XML_INT_HANDLER(name, object, member, minVal, maxVal) { name, XML_Handler_Int, XML_OBJECT_MEMBER(object, member), { NULL, (void *)minVal, (void *)maxVal } }

#define XML_UINT_HANDLER(name, object, member, minVal, maxVal) { name, XML_Handler_UInt, XML_OBJECT_MEMBER(object, member), { NULL, (void *)minVal, (void *)maxVal } }

#define XML_FLOAT_HANDLER(name, object, member, minVal, maxVal) { name, XML_Handler_Float, XML_OBJECT_MEMBER(object, member), { NULL, (void *)minVal, (void *)maxVal } }

#define XML_DOUBLE_HANDLER(name, object, member, minVal, maxVal) { name, XML_Handler_Double, XML_OBJECT_MEMBER(object, member), { NULL, (void *)minVal, (void *)maxVal } }

#define XML_BOOL_HANDLER(name, object, member) { name, XML_Handler_Bool, XML_OBJECT_MEMBER(object, member), { NULL, NULL, NULL } }

#define XML_LIST_HANDLER(name, object, member, list, listSize) { name, XML_Handler_List, XML_OBJECT_MEMBER(object, member), { NULL, (void *)list, (void *)listSize } }

#define XML_STRING_HANDLER(name, object, member, maxLen) { name, XML_Handler_CString, XML_OBJECT_MEMBER(object, member), { NULL, (void *)maxLen, (void *)NULL } }

#define XML_CHAIN_HANDLER(handlers, userData) { NULL, XML_Handler_Chain, 0, 0, { (void *)handlers, (void *)userData, NULL } }

#define XML_CHAIN_HANDLER_MEMBER(handlers, object, member) { NULL, XML_Handler_Chain, XML_OBJECT_MEMBER(object, member), { (void *)handlers, NULL, NULL } }

#define XML_HANDLER_END			{ NULL, XML_Handler_None }

struct XML_Input_ *XML_InputCreate(struct XML_InputStream_ *stream);
void XML_InputFree(struct XML_Input_ *input);
void XML_InputSetUserData(struct XML_Input_ *input, void *userData);
void *XML_InputGetUserData(const struct XML_Input_ *input);
XML_Error XML_InputProcess(struct XML_Input_ *input, const XML_Handler handlers[], void *userData);
XML_Error XML_InputGetError(const struct XML_Input_ *input);
int XML_InputGetLine(const struct XML_Input_ *input);
int XML_InputGetColumn(const struct XML_Input_ *input);
int XML_InputGetOffset(const struct XML_Input_ *input);
const XML_Char *XML_InputGetErrorString(XML_Error error);

XML_Error XML_ElementProcess(struct XML_Element_ *elem, const XML_Handler handlers[], void *userData);
XML_Error XML_ElementReadData(struct XML_Element_ *elem, XML_Char *data, size_t *readSize);
const XML_Char *XML_ElementGetName(const struct XML_Element_ *elem);
const struct XML_Input_ *XML_ElementGetInput(const struct XML_Element_ *elem);
int XML_ElementGetLevel(const struct XML_Element_ *elem);
int XML_ElementIsEmpty(const struct XML_Element_ *elem);
XML_Error XML_ElementGetError(const struct XML_Element_ *elem);

int XML_ElementGetNumAttrs(const struct XML_Element_ *elem);
const XML_Char *XML_ElementGetAttrName(const struct XML_Element_ *elem, int index);
const XML_Char *XML_ElementGetAttrValue(const struct XML_Element_ *elem, int index);
const struct XML_Attribute_ *XML_ElementFindAttr(const struct XML_Element_ *elem, const XML_Char *name);
const struct XML_Attribute_ *XML_ElementGetAttrList(const struct XML_Element_ *elem);
const struct XML_Attribute_ *XML_AttrGetNext(const struct XML_Attribute_ *attr);
const XML_Char *XML_AttrGetName(const struct XML_Attribute_ *attr);
const XML_Char *XML_AttrGetValue(const struct XML_Attribute_ *attr);
int XML_AttrGetInt(const struct XML_Attribute_ *attr, int defValue);
unsigned int XML_AttrGetUInt(const struct XML_Attribute_ *attr, unsigned int defValue);
float XML_AttrGetFloat(const struct XML_Attribute_ *attr, float defValue);
double XML_AttrGetDouble(const struct XML_Attribute_ *attr, double defValue);
int XML_AttrGetBoolean(const struct XML_Attribute_ *attr, int defValue);
int XML_IsWhiteSpace(XML_Char c);
int XML_StringsMatch(const XML_Char *s1, const XML_Char *s2);

#ifdef __cplusplus
}	// extern "C"

// ------------------------------------------------------------------------------
// BEGIN C++ INTERFACE
// ------------------------------------------------------------------------------

XML_BEGIN_NAMESPACE

class InputStream;	// xmlstream.h
class Input;
class Handler;
class Element;

typedef ::XML_Error Error;
typedef ::XML_Char Char;

inline bool IsWhiteSpace(XML_Char c) { return XML_IsWhiteSpace(c) != 0; }
inline bool StringsMatch(const XML_Char *s1, const XML_Char *s2) { return XML_StringsMatch(s1, s2) != 0; }

/** pointer to a function that handles elements
	Applications must implement element handler callbacks using this signature
	@param element	the current element
	@param userData the user-specific data passed in the Handler data or Parse() method
*/
typedef void (*HandlerProc)(Element &element, void *userData);

/** pointer to a function that handles element data
	Applications must implement data handler callbacks using this signature
	@param data 	pointer to the data
	@param len		number of characters of data
	@param userData the user-specific data passed in the Handler data or Parse() method
*/
typedef void (*DataProc)(const XML_Char *data, size_t len, void *userData);

/**
	a Handler for specific events. Applications should pass an array of these
	(terminated with the special Handler::END handler) to the Parser::Parse()
	and Element::Parse() methods. If Handler-specific user-data is supplied, it
	is passed to the HandlerProc, otherwise the user-data passed to the Parse()
	method is passed.
*/
class Handler : public ::XML_Handler
{
friend class Input;
public:
	/** specify an Element handler
		@param elemName		the element name to handle
		@param proc			the HandlerProc to call
		@param userData		the data to pass to the HandlerProc (or NULL)
	*/
	Handler(const XML_Char *elemName, HandlerProc proc, void *userData = NULL);

	/** specify a generic handler that is called when an element is encountered
		that does not have a specific handler.
		This can be used in factory situations where elements
		are named with object types, and the handler can instantiate the
		object by element name via a factory.
		@param proc			the HandlerProc to call
		@param userData		the data to pass to the HandlerProc (or NULL)
	*/
	Handler(HandlerProc proc, void *userData = NULL);

	/** specify a handler that is called when element data is read
		@param proc			the DataProc to call
		@param userData		the data to pass to the DataProc (or NULL)
	*/
	Handler(DataProc proc, void *userData = NULL);

	/** specify a handler that is called when data of the desired type is read
		@param typw			XML_Handler_[Data | CDATA | Comment]
		@param proc			the DataProc to call
		@param userData		the data to pass to the DataProc (or NULL)
	*/
	Handler(XML_HandlerType type, DataProc proc, void *userData = NULL);

	/** specify a separate list of handlers as if the list was embedded at
		this point in the handler list. This allows inserting a handler
		list from a parent class.
		@param handlers		the list of handlers to include
		@param userData		the user-data to use for the handlers
	*/
	Handler(const Handler handlers[], void *userData = NULL);

	/** specify a handler for an element containing a single (signed) integer value.
		@param elemName		the element name to handle
		@param value		the address to place the value
		@param minVal		the minimal value to clamp to, (if minVal or maxVal != 0)
		@param maxVal		the maximum value to clamp to, (if minVal or maxVal != 0)
	*/
	Handler(const XML_Char *elemName, int *value, int minVal = 0, int maxVal = 0);
	Handler(const XML_Char *elemName, unsigned int *value, unsigned int minVal = 0, unsigned int maxVal = 0);
	Handler(const XML_Char *elemName, float *value, float minVal = 0.0f, float maxVal = 0.0f);
	Handler(const XML_Char *elemName, double *value, double *minVal = NULL, double *maxVal = NULL);
	Handler(const XML_Char *elemName, bool *value);
	Handler(const XML_Char *elemName, int *value, const XML_Char * const*list, size_t size);
	Handler(const XML_Char *elemName, XML_Char *value, size_t maxLen);
	Handler(const XML_Char *elemName, XML_HandlerType type, size_t offset, size_t size);
	Handler(const XML_Char *elemName, const XML_Char *const list[], int listSize, size_t offset, size_t size);

	/// return the name of the element the handler handles (or NULL if generic)
	const XML_Char *GetName() const;

	/// special handler to specify the end of a handler list
	static const Handler END;

private:
	Handler();
};

/**
	Represents a single attribute containing a name and value. Use GetNext()
	to get to the next attribute.
*/
class Attribute
{
friend class Element;
public:
	/// copy constructor
	Attribute(const Attribute &rhs);

	/// assignment operator
	Attribute &operator=(const Attribute &rhs);

	/// returns false if the attribute is past the end of the attribute list
	operator bool() const;

	/// return the next attribute
	Attribute GetNext() const;

	/// return the attribute name
	const XML_Char *GetName() const;

	/// return the attribute value
	const XML_Char *GetValue() const;

private:
	Attribute(const struct ::XML_Attribute_ *attr);
	const struct XML_Attribute_ *attr;
};

/**
	Represents a single element, which possibly contains zero or more
	Attributes, zero or more child Elements, or data
*/
class Element
{
friend class Input;
public:
	/// return the element name
	const XML_Char *GetName() const;

	/// return the number of attributes found
	int NumAttributes() const;

	Attribute GetAttrList() const;

	/// return the attribute value with the specified name
	void GetAttribute(const XML_Char *name, XML_Char *str) const;

	/** find the specified attribute and return it, using the default
		value if the attribute is not found.
	*/
	const XML_Char *GetAttribute(const XML_Char *name, const XML_Char *def = NULL) const;

	/** find the specified number attribute and return it, using the default
		value if the argument is not found.
	*/
	void GetAttribute(const XML_Char *name, int &value, int defValue = 0) const;
	void GetAttribute(const XML_Char *name, unsigned int &value, unsigned int defValue = 0) const;
	void GetAttribute(const XML_Char *name, double &value, double defValue = 0.0) const;
	void GetAttribute(const XML_Char *name, float &value, float defValue = 0.0f) const;

	/** find the specified boolean attribute and return it, using the default
		value if the argument is not found. If the attribute value is "true"
		or "True", true is returned, otherwise false
	*/
	void GetAttribute(const XML_Char *name, bool &value, bool defValue = false) const;

	/** read up to len characters of data, or until an element tag is found (<)
		@param buf	buffer to read characters into
		@param len  maximum number of characters to read
		@return 	the number of characters actually read
	*/
	size_t ReadData(XML_Char *buf, size_t len);

	/// parse nested child elements based on new set of handlers
	void Process(const Handler handlers[], void *userData);

	/// return the Input object associated with this element
	const Input &GetInput() const;

	bool IsEmpty() const;

private:
	Element(struct XML_Element_ *element);

	struct XML_Element_ *element;
};

/**
	Parse an XML input stream using per-element handler callbacks
*/
class Input
{
friend class Element;
friend class ParseException;
public:
	/// constructor
	Input(InputStream &stream);

	/// destructor
	~Input();

	/// parse elements based on set of handlers
	void Process(const Handler handlers[], void *userData);

	/// return the current line number
	int GetLine() const;

	/// return the current column number
	int GetColumn() const;

	/// return the current byte offset
	int GetOffset() const;

	/// set user-specific data that can be queried with GetUserData
	void SetUserData(void *userData);

	/// return the user-specific data set with SetUserData
	void *GetUserData() const;

	/// return the last error code
	XML_Error GetError() const;

private:
	struct InputImpl *impl;	// implementation
	static XML_Error elementHandler(struct ::XML_Input_ *input, struct ::XML_Element_ *elem, const ::XML_Handler *handler, void *userData);
	static XML_Error dataHandler(struct ::XML_Input_ *input, const ::XML_Char *data, size_t len, const ::XML_Handler *handler, void *userData);
};

/**
	All parsing-related exceptions are derived from ParseException.
*/
class ParseException
{
friend class Input;
public:
	ParseException(const Input &input);

	/// return the line the exception occured at
	int GetLine() const { return line; }
	/// return the column the exception occured at
	int GetColumn() const { return column; }
	/// return the byte offset the exception occured at
	int GetOffset() const { return offset; }

	XML_Error GetError() const { return error; }

	/// return a string describing the error
	virtual const XML_Char *What() const;

protected:
	XML_Error error;
	int line;
	int column;
	int offset;
};

class InvalidValue : public ParseException
{
public:
	InvalidValue(const Input &input);
	InvalidValue(const Input &input, int line, int column);
};

XML_END_NAMESPACE
#endif	// __cplusplus

#endif	// XMLINPUT_H


