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
#include "xmlio/xmlinput.h"
#include "xmlio/xmlstream.h"
#include <assert.h>	// assert()
#include <stdio.h>	// sprintf()
#include <string.h>	// memset()

XML_BEGIN_NAMESPACE

struct InputImpl
{
	::XML_InputStream adapter;
	InputStream *stream;
	::XML_Input *input;
	void *userData;
};

// convert XML_Element to XML::Element and handle exceptions
::XML_Error Input::elementHandler(::XML_Input *input, ::XML_Element *elem, const ::XML_Handler *handler, void *userData)
{
	Element e(elem);
	HandlerProc proc = (HandlerProc)handler->info.Element.proc;
	try {
		(*proc)(e, handler->info.Element.userData ? handler->info.Element.userData : userData);
	}
	catch (const ParseException &e)
	{
		return e.GetError();
	}
	return ::XML_Error_None;
}

::XML_Error Input::dataHandler(::XML_Input *input, const ::XML_Char *data, size_t len, const ::XML_Handler *handler, void *userData)
{
	DataProc proc = (DataProc)handler->info.Data.proc;
	(*proc)(data, len, userData);
	return ::XML_Error_None;
}

static int sInputReadProc(::XML_InputStream *str, ::XML_Char *buf, size_t bufLen)
{
	InputImpl *impl = (InputImpl *)str;
	assert(impl);
	assert(buf);
	return impl->stream->read(buf, bufLen);
}

Input::Input(InputStream &stream) : impl(new InputImpl)
{
	impl->stream = &stream;
	impl->adapter.readProc = sInputReadProc;
	impl->input = ::XML_InputCreate(&impl->adapter);
	impl->input->elementHandler = elementHandler;
	impl->input->dataHandler = dataHandler;
	impl->input->userData = this;
	impl->userData = NULL;
}

Input::~Input()
{
	::XML_InputFree(impl->input);
	delete impl;
}

void Input::Process(const Handler handlers[], void *userData)
{
	::XML_Error error = ::XML_InputProcess(impl->input, (const ::XML_Handler *)handlers, userData);
	if (error != ::XML_Error_None)
		throw ParseException(*this);
}

int Input::GetLine() const
{
	return ::XML_InputGetLine(impl->input);
}

int Input::GetColumn() const
{
	return ::XML_InputGetColumn(impl->input);
}

int Input::GetOffset() const
{
	return ::XML_InputGetOffset(impl->input);
}

void Input::SetUserData(void *userData)
{
	impl->userData = userData;
}

void *Input::GetUserData() const
{
	return impl->userData;
}

::XML_Error Input::GetError() const
{
	return ::XML_InputGetError(impl->input);
}

//
// Handler
//

const Handler Handler::END;

Handler::Handler()
{
	name = NULL;
	type = ::XML_Handler_None;
}

Handler::Handler(const ::XML_Char *elemName, HandlerProc proc, void *userData)
{
	assert(elemName);
	assert(proc);

	name = elemName;
	type = ::XML_Handler_Element;
	offset = size = 0;
	info.Element.proc = (::XML_HandlerProc)proc;
	info.Element.userData = userData;
}

Handler::Handler(HandlerProc proc, void *userData)
{
	assert(proc);

	name = NULL;
	type = ::XML_Handler_Element;
	offset = size = 0;
	info.Element.proc = (::XML_HandlerProc)proc;
	info.Element.userData = userData;
}

Handler::Handler(DataProc proc, void *userData)
{
	assert(proc);

	name = NULL;
	type = ::XML_Handler_Data;
	offset = size = 0;
	info.Data.proc = (::XML_DataProc)proc;
	info.Data.userData = userData;
}

Handler::Handler(XML_HandlerType type, DataProc proc, void *userData)
{
	assert(proc);
	assert(type == XML_Handler_Data || type == XML_Handler_CDATA || type == XML_Handler_Comment);

	name = NULL;
	this->type = type;
	offset = size = 0;
	info.Data.proc = (::XML_DataProc)proc;
	info.Data.userData = userData;
}

Handler::Handler(const Handler handlers[], void *userData)
{
	assert(handlers);

	name = NULL;
	type = ::XML_Handler_Chain;
	offset = size = 0;
	info.Chain.handlers = handlers;
	info.Chain.userData = userData;
}

Handler::Handler(const ::XML_Char *elemName, int *value, int minVal, int maxVal)
{
	assert(elemName);
	assert(value);

	name = elemName;
	type = ::XML_Handler_Int;
	offset = 0;
	size = sizeof(int);
	info.Int.result = value;
	info.Int.minVal = minVal;
	info.Int.maxVal = maxVal;
}

Handler::Handler(const ::XML_Char *elemName, unsigned int *value, unsigned int minVal, unsigned int maxVal)
{
	assert(elemName);
	assert(value);

	name = elemName;
	type = ::XML_Handler_UInt;
	offset = 0;
	size = sizeof(unsigned int);
	info.UInt.result = value;
	info.UInt.minVal = minVal;
	info.UInt.maxVal = maxVal;
}

Handler::Handler(const ::XML_Char *elemName, float *value, float minVal, float maxVal)
{
	assert(elemName);
	assert(value);

	name = elemName;
	type = ::XML_Handler_Float;
	offset = 0;
	size = sizeof(float);
	info.Float.result = value;
	info.Float.minVal = minVal;
	info.Float.maxVal = maxVal;
}

Handler::Handler(const ::XML_Char *elemName, double *value, double *minVal, double *maxVal)
{
	assert(elemName);
	assert(value);

	name = elemName;
	type = ::XML_Handler_Double;
	offset = 0;
	size = sizeof(double);
	info.Double.result = value;
	info.Double.minVal = minVal;
	info.Double.maxVal = maxVal;
}

Handler::Handler(const ::XML_Char *elemName, bool *value)
{
	assert(elemName);
	assert(value);

	name = elemName;
	type = ::XML_Handler_Bool;
	offset = 0;
	size = sizeof(bool);
	info.Bool.result = (int *)value;
}

Handler::Handler(const ::XML_Char *elemName, int *value, const ::XML_Char * const*list, size_t size)
{
	assert(elemName);
	assert(value);

	name = elemName;
	type = ::XML_Handler_List;
	offset = 0;
	size = sizeof(int);
	info.List.result = value;
	info.List.list = list;
	info.List.listSize = size;
}

Handler::Handler(const ::XML_Char *elemName, ::XML_Char *value, size_t maxLen)
{
	assert(elemName);
	assert(value);

	name = elemName;
	type = ::XML_Handler_CString;
	offset = 0;
	size = sizeof(::XML_Char);
	info.CString.result = value;
	info.CString.maxLen = maxLen;
}

Handler::Handler(const ::XML_Char *elemName, ::XML_HandlerType type, size_t offset, size_t size)
{
	assert(elemName);

	::memset(this, 0, sizeof(Handler));

	name = elemName;
	this->type = type;
	this->offset = offset;
	this->size = size;
}

Handler::Handler(const ::XML_Char *elemName, const ::XML_Char *const list[], int listSize, size_t offset, size_t size)
{
	assert(elemName);
	assert(list);

	name = elemName;
	this->type = ::XML_Handler_List;
	this->offset = offset;
	this->size = size;
	this->info.List.result = NULL;
	this->info.List.list = list;
	this->info.List.listSize = listSize;
}

const ::XML_Char *Handler::GetName() const
{
	return name;
}

Attribute::Attribute(const ::XML_Attribute *attr) : attr(attr)
{
}

Attribute::Attribute(const Attribute &rhs) : attr(rhs.attr)
{
}

Attribute &Attribute::operator=(const Attribute &rhs)
{
	attr = rhs.attr;
	return *this;
}

Attribute Attribute::GetNext() const
{
	assert(attr);
	return Attribute(::XML_AttrGetNext(attr));
}

const XML_Char *Attribute::GetName() const
{
	assert(attr);
	return ::XML_AttrGetName(attr);
}

const XML_Char *Attribute::GetValue() const
{
	assert(attr);
	return ::XML_AttrGetValue(attr);
}

Attribute::operator bool() const
{
	return attr != NULL;
}


//
// Element
//

Element::Element(::XML_Element *elem)
{
	assert(elem);
	element = elem;
}

const ::XML_Char *Element::GetName() const
{
	return element->name;
}

int Element::NumAttributes() const
{
	return ::XML_ElementGetNumAttrs(element);
}

Attribute Element::GetAttrList() const
{
	return Attribute(::XML_ElementGetAttrList(element));
}

const ::XML_Char *Element::GetAttribute(const ::XML_Char *name, const ::XML_Char *def) const
{
	const ::XML_Attribute *attr = ::XML_ElementFindAttr(element, name);
	return attr ? attr->value : def;
}

void Element::GetAttribute(const ::XML_Char *name, int &value, int defValue) const
{
	value = ::XML_AttrGetInt(::XML_ElementFindAttr(element, name), defValue);
}

void Element::GetAttribute(const ::XML_Char *name, unsigned int &value, unsigned int defValue) const
{
	value = ::XML_AttrGetUInt(::XML_ElementFindAttr(element, name), defValue);
}

void Element::GetAttribute(const ::XML_Char *name, double &value, double defValue) const
{
	value = ::XML_AttrGetDouble(::XML_ElementFindAttr(element, name), defValue);
}

void Element::GetAttribute(const ::XML_Char *name, float &value, float defValue) const
{
	value = ::XML_AttrGetFloat(::XML_ElementFindAttr(element, name), defValue);
}

void Element::GetAttribute(const ::XML_Char *name, bool &value, bool defValue) const
{
	value = ::XML_AttrGetBoolean(::XML_ElementFindAttr(element, name), defValue) != 0;
}

size_t Element::ReadData(::XML_Char *buf, size_t len)
{
	::XML_Error error = ::XML_ElementReadData(element, buf, &len);
	if (error != ::XML_Error_None)
		throw ParseException(GetInput());
	return len;
}

void Element::Process(const Handler handlers[], void *userData)
{
	::XML_Error error = ::XML_ElementProcess(element, handlers, userData);
	if (error != ::XML_Error_None)
		throw ParseException(GetInput());
}

const Input &Element::GetInput() const
{
	Input *input = (Input *)element->input->userData;
	assert(input);
	return *input;
}

bool Element::IsEmpty() const
{
	return element->empty != 0;
}

//
// Exceptions
//

ParseException::ParseException(const Input &input)
{
	error = input.GetError();
	line = input.GetLine();
	column = input.GetColumn();
	offset = input.GetOffset();
}

const ::XML_Char *ParseException::What() const
{
	return ::XML_InputGetErrorString(error);
}

InvalidValue::InvalidValue(const Input &input) : ParseException(input)
{
	// possibly thrown by user code so assign error code here
	this->error = ::XML_Error_InvalidValue;
}

InvalidValue::InvalidValue(const Input &input, int line, int column) : ParseException(input)
{
	// possibly thrown by user code so assign error code here
	this->error = ::XML_Error_InvalidValue;
	this->line = line;
	this->column = column;
}

XML_END_NAMESPACE

