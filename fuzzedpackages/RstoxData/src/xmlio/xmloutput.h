/*
	Utility class for outputting XML data files.

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

#ifndef XMLOUT_H
#define XMLOUT_H

#include "xmlstream.h"

#ifdef __cplusplus
#include <vector>
#include <string>

XML_BEGIN_NAMESPACE

/**
	Utility class for saving XML data files to an OutputStream.
	Keeps track of element nesting level and assists with writing
	attribute lists.
*/
class Output
{
public:
	/// Constructor - initialize and Output object from an OutputStream
	Output(OutputStream &stream);

	/** Write the beginning ?xml processing instruction - call before
		any other method
	*/
	void BeginDocument(const char *version = "1.0", 
					const char *encoding = "US-ASCII", 
					bool standalone = true);

	/** Does nothing except performs sanity checking on nesting level -
		call after finished writing data
	*/
	void EndDocument();

	enum Mode
	{
		indent,	/// normal nested indenting
		terse	/// do not add newlines or indenting (for compact elements)
	};

	/**
		Begin a new Element, indenting to the next nesting level if mode is 'indent'
		@param name		the element name
		@param mode		if terse, do not write a newline after
	*/
	void BeginElement(const char *name, Mode mode = indent);
	/**
		Begin a new Element, leaving the element tag open for attributes
		@param name		the element name
	*/
	void BeginElementAttrs(const char *name);
	/**
		End an element tag after writing attributes
		@param mode		if terse, do not write a newline after
	*/
	void EndAttrs(Mode mode = indent);

	/// write the std::string attribute
	void WriteAttr(const char *name, const std::string &value);
	/// write a "C" string attribute
	void WriteAttr(const char *name, const char *value);
	/// write an integer attribute
	void WriteAttr(const char *name, int value);
	/// write a floating-point attribute
	void WriteAttr(const char *name, double value);
	/// write a boolean attribute ("true" or "false")
	void WriteAttr(const char *name, bool value);
	/** finish writing the current element by writing the end-element tag
		@param mode		if terse, do not write a newline or indent before
	*/
	void EndElement(Mode mode = indent);
	/// explicitely indent to the current nest level (when writing element data)
	void Indent();

	/// write out a "terse" element with the specified data
	void WriteElement(const char *name, const std::string &value);
	/// write out a "terse" element with the specified data
	void WriteElement(const char *name, const char *value);
	/// write out a "terse" element with the specified data
	void WriteElement(const char *name, int value);
	/// write out a "terse" element with the specified data
	void WriteElement(const char *name, unsigned int value);
	/// write out a "terse" element with the specified data
	void WriteElement(const char *name, double value);
	/// write out a "terse" element with the specified data
	void WriteElement(const char *name, bool value);

	Output &operator<<(const std::string &str);
	Output &operator<<(const char *str);
	Output &operator<<(int value);
	Output &operator<<(unsigned int value);
	Output &operator<<(double value);
	Output &operator<<(bool value);

	void write(const char *str, size_t len);
	void writeString(const char *str);
	void writeLine(const char *str);

private:
	OutputStream &mStream;
	int mLevel;					// nesting level
	typedef std::vector<const char *> ElementStack;
	ElementStack mElements;		// needed to write EndElement tag name
	bool mAttributes;			// used for sanity-checking
};

XML_END_NAMESPACE

#endif	/* __cplusplus */

#endif	// XMLOUT_H
